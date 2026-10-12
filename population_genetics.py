"""
Population Genetics Statistical Analysis for Diploid SNP Genotype Data.
Implements: Ho/He/PIC, an exact conditional HWE test for biallelic
diploid SNPs, Weir & Cockerham (1984) ANOVA variance-component F-statistics,
sample-aware pairwise-difference AMOVA (Excoffier et al. 1992, Arlequin-compatible), and Bonferroni + BH-FDR
multiple-testing correction. Multiallelic HWE testing is not supported in v1.0.0.

F-statistics reference:
  Weir BS, Cockerham CC. Estimating F-statistics for the analysis of
  population structure. Evolution. 1984;38(6):1358-1370.
"""

import random as _random_mod
from typing import Dict, List, Tuple, Optional
from collections import Counter, defaultdict
from dataclasses import dataclass

# ---------- Dataclasses ----------


@dataclass
class LocusStats:
    locus_name: str
    population: str
    n_individuals: int
    n_alleles_observed: int
    allele_frequencies: Dict[float, float]
    ho: float
    he: float
    he_unbiased: float
    pic: float
    hwe_p_value: Optional[float]
    hwe_method: str
    fis: Optional[float]
    null_allele_freq: Optional[float]


@dataclass
class LocusSummary:
    locus_name: str
    populations: List[str]
    overall_ho: float
    overall_he: float
    overall_pic: float
    overall_n_alleles: int
    per_population: Dict[str, LocusStats]
    fst: Optional[float]
    fit: Optional[float]
    fst_p_value: Optional[float]
    fst_method: str = ""
    permutation_count: int = 0
    random_seed: int = 0


@dataclass
class AMOVAResult:
    source: str
    df: int
    ss: float
    ms: float
    variance_component: float
    percentage: float
    phi_stat: Optional[float] = None
    p_value: Optional[float] = None
    method: str = "pairwise_difference"
    n_individuals: int = 0
    n_loci: int = 0
    n_permutations: int = 0
    random_seed: int = 42
    warning: Optional[str] = None


@dataclass
class GenotypeFreqTestResult:
    locus_name: str
    chi2_stat: float
    degrees_of_freedom: int
    raw_p_value: Optional[float]
    p_value: Optional[float]
    bonferroni_p_value: Optional[float]
    fdr_q_value: Optional[float]
    significant_raw: bool
    significant_bonferroni: bool
    significant_fdr: bool
    significant: bool
    genotype_labels: List[str]
    populations: List[str]
    freq_table: List[List[float]]
    original_contingency_table: Optional[List[List[int]]] = None
    merged_contingency_table: Optional[List[List[int]]] = None
    merged_genotype_map: Optional[Dict[str, str]] = None
    test_method: str = "chi2"
    warning: Optional[str] = None


# ---------- Core Engine ----------


class PopulationGenetics:
    """Population genetics for diploid SNP markers."""

    @staticmethod
    def compute_allele_frequencies(
        genotypes: List[Tuple[float, ...]],
    ) -> Tuple[Dict[float, float], int]:
        c = Counter()
        total = 0
        for gt in genotypes:
            for a in gt:
                c[a] += 1
                total += 1
        if total == 0:
            return {}, 0
        return {k: v / total for k, v in c.items()}, total

    @staticmethod
    def compute_heterozygosity(
        genotypes: List[Tuple[float, ...]],
    ) -> Tuple[float, float, int]:
        n = len(genotypes)
        if n == 0:
            return 0.0, 0.0, 0
        ho = sum(1 for gt in genotypes if len(gt) == 2 and gt[0] != gt[1]) / n
        freqs, _ = PopulationGenetics.compute_allele_frequencies(genotypes)
        he = 1.0 - sum(f**2 for f in freqs.values())
        if n > 1:
            he = (2.0 * n) / (2.0 * n - 1.0) * he
        return ho, he, n

    @staticmethod
    def compute_pic(freqs: Dict[float, float]) -> float:
        p = list(freqs.values())
        k = len(p)
        if k <= 1:
            return 0.0
        ss = sum(x**2 for x in p)
        ct = sum(2.0 * p[i] ** 2 * p[j] ** 2 for i in range(k) for j in range(i + 1, k))
        return max(0.0, 1.0 - ss - ct)

    # ---------- HWE ----------

    @staticmethod
    def _valid_hwe_genotypes(genotypes):
        """Return complete, non-missing diploid genotypes for HWE testing.

        HWE analysis in v1.0.0 is restricted to complete diploid calls. Missing
        and malformed calls are ignored. Allele labels may be strings or
        numbers, but they must be hashable.
        """
        missing_strings = {"", "-", "--", "N", "NN", "N/A", "NA", "NAN"}
        valid = []

        for genotype in genotypes or []:
            if genotype is None:
                continue

            try:
                if len(genotype) != 2:
                    continue
                allele1, allele2 = genotype
            except (TypeError, ValueError):
                continue

            if allele1 is None or allele2 is None:
                continue

            if isinstance(allele1, str) and allele1.strip().upper() in missing_strings:
                continue
            if isinstance(allele2, str) and allele2.strip().upper() in missing_strings:
                continue

            try:
                hash(allele1)
                hash(allele2)
            except TypeError:
                continue

            valid.append((allele1, allele2))

        return valid

    @staticmethod
    def _hwe_exact_biallelic(obs_hom1: int, obs_hets: int, obs_hom2: int) -> float:
        """Compute the two-sided exact HWE p-value for a biallelic locus.

        This uses the heterozygote-probability recurrence described by
        Wigginton et al. (2005). The test conditions on the observed allele
        counts and sums probabilities no greater than that of the observed
        genotype configuration.
        """
        if min(obs_hom1, obs_hets, obs_hom2) < 0:
            raise ValueError("Observed genotype counts must be non-negative")

        n_genotypes = obs_hom1 + obs_hets + obs_hom2
        if n_genotypes == 0:
            return 1.0

        obs_hom_rare = min(obs_hom1, obs_hom2)
        rare_copies = 2 * obs_hom_rare + obs_hets

        probabilities = [0.0] * (rare_copies + 1)

        midpoint = int(
            rare_copies * (2 * n_genotypes - rare_copies) / (2 * n_genotypes)
        )
        if midpoint % 2 != rare_copies % 2:
            midpoint += 1

        current_hets = midpoint
        current_hom_rare = (rare_copies - midpoint) // 2
        current_hom_common = n_genotypes - current_hets - current_hom_rare

        probabilities[midpoint] = 1.0
        probability_sum = 1.0

        # Recur toward fewer heterozygotes.
        while current_hets >= 2:
            next_probability = (
                probabilities[current_hets]
                * current_hets
                * (current_hets - 1)
                / (4.0 * (current_hom_rare + 1) * (current_hom_common + 1))
            )
            probabilities[current_hets - 2] = next_probability
            probability_sum += next_probability

            current_hets -= 2
            current_hom_rare += 1
            current_hom_common += 1

        # Reset and recur toward more heterozygotes.
        current_hets = midpoint
        current_hom_rare = (rare_copies - midpoint) // 2
        current_hom_common = n_genotypes - current_hets - current_hom_rare

        while current_hets <= rare_copies - 2:
            next_probability = (
                probabilities[current_hets]
                * 4.0
                * current_hom_rare
                * current_hom_common
                / ((current_hets + 2) * (current_hets + 1))
            )
            probabilities[current_hets + 2] = next_probability
            probability_sum += next_probability

            current_hets += 2
            current_hom_rare -= 1
            current_hom_common -= 1

        if probability_sum <= 0.0:
            return 1.0

        probabilities = [value / probability_sum for value in probabilities]
        observed_probability = probabilities[obs_hets]
        tolerance = 1e-12
        p_value = sum(
            value
            for value in probabilities
            if value <= observed_probability + tolerance
        )
        return max(0.0, min(1.0, p_value))

    @staticmethod
    def hwe_exact_test(
        genotypes: List[Tuple[float, ...]], max_perm: int = 100000
    ) -> Tuple[Optional[float], str]:
        """Perform an exact conditional HWE test for a biallelic diploid SNP.

        ``max_perm`` is retained only for backward API compatibility; this
        deterministic exact test does not use permutations or MCMC.

        Returns
        -------
        (p_value, method)
            ``empty`` when no complete diploid calls remain;
            ``monomorphic`` for a one-allele locus;
            ``exact-biallelic`` for a two-allele locus; and
            ``unsupported-multiallelic`` when more than two alleles occur.
        """
        del max_perm  # Backward-compatible argument; intentionally unused.

        valid_genotypes = PopulationGenetics._valid_hwe_genotypes(genotypes)
        if not valid_genotypes:
            return None, "empty"

        alleles = []
        for allele1, allele2 in valid_genotypes:
            if allele1 not in alleles:
                alleles.append(allele1)
            if allele2 not in alleles:
                alleles.append(allele2)

        if len(alleles) == 1:
            return 1.0, "monomorphic"
        if len(alleles) > 2:
            return None, "unsupported-multiallelic"

        allele1, allele2 = alleles
        obs_hom1 = 0
        obs_hets = 0
        obs_hom2 = 0

        for first, second in valid_genotypes:
            if first == allele1 and second == allele1:
                obs_hom1 += 1
            elif first == allele2 and second == allele2:
                obs_hom2 += 1
            else:
                obs_hets += 1

        p_value = PopulationGenetics._hwe_exact_biallelic(obs_hom1, obs_hets, obs_hom2)
        return p_value, "exact-biallelic"

    # ---------- Weir & Cockerham (1984) F-statistics ----------
    # Ref: Weir BS, Cockerham CC. Evolution 38(6):1358-1370, 1984.
    #
    # WC84 estimates Fst (theta), Fis (f), and Fit (F) via method-of-moments
    # ANOVA variance components:
    #   a = variance among populations
    #   b = variance among individuals within populations
    #   w = variance between gametes within individuals (= c in WC84)
    #
    # For multiallelic loci, variance components are summed over all alleles
    # and the ratios are computed from the totals.

    @staticmethod
    def _wc84_one_allele(
        allele_tag: float,
        genotypes_by_pop: Dict[str, List[Tuple[float, ...]]],
    ) -> Tuple[float, float, float, int, int]:
        """WC84 variance components for a single binary-encoded allele.

        Each genotype tuple is re-encoded as the count (0, 1, or 2) of
        *allele_tag* copies.  The three ANOVA sums of squares (among
        populations, among individuals, within individuals) are computed and
        converted to variance components a, b, w.

        Returns (a, b, w, r, N) where r = number of populations, N = total
        sample size (individuals).
        """
        n_pop = []   # sample sizes per population
        A_pop = []   # total copies of allele_tag per population
        p_pop = []   # frequency of allele_tag per population

        for genotypes in genotypes_by_pop.values():
            # Binary-encode: count how many copies of allele_tag each diploid has
            counts = []
            for gt in genotypes:
                if len(gt) != 2:
                    continue
                c = (1 if gt[0] == allele_tag else 0) + (1 if gt[1] == allele_tag else 0)
                counts.append(c)
            ni = len(counts)
            if ni < 2:
                continue
            Ai = sum(counts)
            pi = Ai / (2.0 * ni)
            n_pop.append(ni)
            A_pop.append(Ai)
            p_pop.append(pi)

        r = len(n_pop)
        if r == 0:
            return 0.0, 0.0, 0.0, 0, 0

        N = sum(n_pop)
        if N == 0:
            return 0.0, 0.0, 0.0, 0, 0

        # Weighted mean frequency
        p_bar = sum(A_pop) / (2.0 * N)

        # --- Within-individual (between-gamete) SS ---
        SSW = 0.0
        for genotypes in genotypes_by_pop.values():
            for gt in genotypes:
                if len(gt) != 2:
                    continue
                c = (1 if gt[0] == allele_tag else 0) + (1 if gt[1] == allele_tag else 0)
                SSW += c * (2 - c) / 2.0

        # --- Total SS ---
        R = sum(A_pop)
        SST = R * (1.0 - p_bar)

        if r < 2:
            # Single population: no among-pop component; only b and w
            SSI = max(0.0, SST - SSW)
            MSI = SSI / float(N - 1) if N > 1 else 0.0
            MSW = SSW / float(N) if N > 0 else 0.0
            b_comp = max(0.0, (MSI - MSW) / 2.0)
            w_comp = max(0.0, MSW)
            return 0.0, b_comp, w_comp, r, N

        # --- Multi-population: full ANOVA ---
        # SSP: among-population sum of squares
        SSP = sum(2.0 * ni * (pi - p_bar) ** 2
                  for ni, pi in zip(n_pop, p_pop))

        # WC84 correction for unequal sample sizes (Eq. 13 in WC84)
        n_C = (N - sum(ni * ni for ni in n_pop) / float(N)) / float(r - 1)
        if n_C <= 0.0:
            return 0.0, 0.0, 0.0, r, N

        SSI = max(0.0, SST - SSP - SSW)   # among individuals within pops

        # Mean squares
        MSG = SSP / float(r - 1)
        MSI = SSI / float(N - r)
        MSW = SSW / float(N)

        # Variance components (WC84)
        a_comp = max(0.0, (MSG - MSI) / (2.0 * n_C))
        b_comp = max(0.0, (MSI - MSW) / 2.0)
        w_comp = max(0.0, MSW)

        return a_comp, b_comp, w_comp, r, N

    @staticmethod
    def compute_f_statistics(
        genotypes_by_pop: Dict[str, List[Tuple[float, ...]]],
    ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """Weir & Cockerham (1984) F-statistics — ANOVA variance components.

        Reference
        ---------
        Weir BS, Cockerham CC.  Estimating F-statistics for the analysis of
        population structure.  *Evolution* 38(6):1358–1370, 1984.

        Method
        ------
        Decomposes total allele-frequency variance into three components
        (a = among populations, b = among individuals within populations,
        w = within individuals / between gametes) and computes:

        .. math::
           F_{ST} = \\theta = a / (a + b + w) \\\\
           F_{IS} = f = b / (b + w) \\\\
           F_{IT} = F = (a + b) / (a + b + w)

        For multiallelic loci, the components are summed over all alleles
        and the ratios computed from the totals (WC84 multiallelic extension).

        Notes
        -----
        Negative variance-component estimates (possible in small samples) are
        clamped to zero.  F-statistics are clipped to [0, 1].

        Returns
        -------
        (Fis, Fst, Fit) or (None, None, None) if the data are insufficient.
        For single-population analyses, Fst and Fit are set to None.
        """
        if len(genotypes_by_pop) == 0:
            return None, None, None

        # Collect all alleles present
        all_genotypes = []
        for gt_list in genotypes_by_pop.values():
            all_genotypes.extend(gt_list)
        allele_freqs, _ = PopulationGenetics.compute_allele_frequencies(all_genotypes)
        alleles = sorted(allele_freqs.keys())

        if len(alleles) < 2:
            return None, None, None   # monomorphic

        # Sum variance components across all alleles
        total_a = 0.0
        total_b = 0.0
        total_w = 0.0
        n_valid_pops = 0
        for allele in alleles:
            a, b, w, r, _ = PopulationGenetics._wc84_one_allele(
                allele, genotypes_by_pop
            )
            total_a += a
            total_b += b
            total_w += w
            n_valid_pops = max(n_valid_pops, r)

        sigma2_total = total_a + total_b + total_w
        if sigma2_total <= 0.0:
            return None, None, None

        Fst = total_a / sigma2_total
        Fis = total_b / (total_b + total_w) if (total_b + total_w) > 0.0 else 0.0
        Fit = (total_a + total_b) / sigma2_total

        # Clamp to [0, 1]
        Fst = max(0.0, min(1.0, Fst))
        Fis = max(0.0, min(1.0, Fis))
        Fit = max(0.0, min(1.0, Fit))

        # Single-population: Fst and Fit not meaningful
        if n_valid_pops <= 1:
            return Fis, None, None

        return Fis, Fst, Fit

    @staticmethod
    def fst_permutation_test(
        genotypes_by_pop: Dict[str, List[Tuple[float, ...]]],
        n_permutations: int = 1000,
        random_seed: int = 42,
    ) -> Tuple[float, float]:
        """Permutation test for WC84 Fst with +1 correction: p=(count+1)/(n_valid+1).

        Population labels are randomly shuffled among individuals, WC84 Fst is
        re-computed on each permuted dataset, and the p-value is the proportion
        of permuted Fst values >= the observed Fst (with +1 pseudo-count).
        """
        ag = []
        pl = []
        for pop, gt in genotypes_by_pop.items():
            for g in gt:
                ag.append(g)
                pl.append(pop)
        if len(set(pl)) <= 1:
            return 0.0, 1.0
        _, ofst, _ = PopulationGenetics.compute_f_statistics(genotypes_by_pop)
        if ofst is None:
            return 0.0, 1.0
        rng = _random_mod.Random(random_seed)
        ce = 0
        nv = 0
        for _ in range(n_permutations):
            sl = list(pl)
            rng.shuffle(sl)
            pp: Dict[str, List[Tuple[float, ...]]] = defaultdict(list)
            for g, p in zip(ag, sl):
                pp[p].append(g)
            _, pf, _ = PopulationGenetics.compute_f_statistics(pp)
            if pf is not None:
                nv += 1
                ce += 1 if pf >= ofst else 0
        pv = (ce + 1) / (nv + 1) if nv > 0 else 1.0
        return ofst, max(0.0, min(1.0, pv))

    # ---------- Sample-aware multi-locus AMOVA ----------

    @staticmethod
    def _valid_diploid_genotype(genotype):
        """Return a complete diploid genotype tuple, or None for missing data."""
        if genotype is None:
            return None
        try:
            if len(genotype) != 2:
                return None
            allele1, allele2 = genotype
        except (TypeError, ValueError):
            return None

        missing_strings = {"", "-", "--", "N", "NN", "N/A", "NA", "NAN"}
        if allele1 is None or allele2 is None:
            return None
        if isinstance(allele1, str) and allele1.strip().upper() in missing_strings:
            return None
        if isinstance(allele2, str) and allele2.strip().upper() in missing_strings:
            return None
        try:
            hash(allele1)
            hash(allele2)
        except TypeError:
            return None
        return allele1, allele2

    @staticmethod
    def _diploid_allele_sharing_distance(genotype1, genotype2) -> float:
        """Return dosage-aware allele-sharing distance for two diploid calls.

        Distance is ``1 - shared_allele_copies / 2``. Counter intersection is
        used so that homozygous allele dosage is retained.
        """
        first = PopulationGenetics._valid_diploid_genotype(genotype1)
        second = PopulationGenetics._valid_diploid_genotype(genotype2)
        if first is None or second is None:
            raise ValueError("Both genotypes must be complete diploid calls")
        shared_copies = sum((Counter(first) & Counter(second)).values())
        return 1.0 - shared_copies / 2.0

    @staticmethod
    def multilocus_distance_matrix(
        sample_genotypes,
        sample_populations,
        min_common_loci: int = 1,
    ):
        """Build a sample-level multi-locus allele-sharing distance matrix.

        Parameters
        ----------
        sample_genotypes
            ``{sample_id: {locus_name: diploid_genotype_or_missing}}``.
        sample_populations
            ``{sample_id: population_name}``.
        min_common_loci
            Minimum number of loci called in both individuals for every
            pairwise distance.

        Returns
        -------
        (sample_ids, population_labels, distance_matrix, common_locus_counts)
        """
        if min_common_loci < 1:
            raise ValueError("min_common_loci must be at least 1")

        retained = []
        normalized = {}
        for sample_id in sorted(sample_genotypes, key=str):
            if sample_id not in sample_populations:
                continue
            calls = {}
            for locus, genotype in sample_genotypes[sample_id].items():
                valid = PopulationGenetics._valid_diploid_genotype(genotype)
                if valid is not None:
                    calls[locus] = valid
            if calls:
                retained.append(sample_id)
                normalized[sample_id] = calls

        if len(retained) < 2:
            raise ValueError("At least two genotyped individuals are required")

        labels = [sample_populations[sample_id] for sample_id in retained]
        if len(set(labels)) < 2:
            raise ValueError("At least two populations are required")

        n_individuals = len(retained)
        distances = [[0.0] * n_individuals for _ in range(n_individuals)]
        common_counts = [[0] * n_individuals for _ in range(n_individuals)]

        for i in range(n_individuals):
            common_counts[i][i] = len(normalized[retained[i]])
            for j in range(i + 1, n_individuals):
                first_calls = normalized[retained[i]]
                second_calls = normalized[retained[j]]
                common_loci = sorted(set(first_calls) & set(second_calls), key=str)
                common_counts[i][j] = len(common_loci)
                common_counts[j][i] = len(common_loci)

                if len(common_loci) < min_common_loci:
                    raise ValueError(
                        f"Samples {retained[i]!r} and {retained[j]!r} share "
                        f"{len(common_loci)} called loci; at least "
                        f"{min_common_loci} are required"
                    )

                locus_distances = [
                    PopulationGenetics._diploid_allele_sharing_distance(
                        first_calls[locus], second_calls[locus]
                    )
                    for locus in common_loci
                ]
                distance = sum(locus_distances) / len(locus_distances)
                distances[i][j] = distance
                distances[j][i] = distance

        return retained, labels, distances, common_counts

    @staticmethod
    def _amova_components_from_distance(distance_matrix, population_labels):
        """Smouse & Peakall (1999) distance-based AMOVA from N-individual matrix."""
        n_individuals = len(population_labels)
        if n_individuals != len(distance_matrix):
            raise ValueError("Distance matrix and population labels differ in size")

        populations = list(dict.fromkeys(population_labels))
        n_populations = len(populations)
        if n_populations < 2:
            raise ValueError("At least two populations are required")
        if n_individuals <= n_populations:
            raise ValueError(
                "Within-population degrees of freedom must be positive; "
                "at least one population needs more than one individual"
            )

        ss_total = (
            sum(
                distance_matrix[i][j] ** 2
                for i in range(n_individuals)
                for j in range(i + 1, n_individuals)
            )
            / n_individuals
        )

        population_indices = {
            population: [
                index
                for index, label in enumerate(population_labels)
                if label == population
            ]
            for population in populations
        }
        population_sizes = {
            population: len(indices)
            for population, indices in population_indices.items()
        }

        ss_within = 0.0
        for indices in population_indices.values():
            if len(indices) > 1:
                ss_within += sum(
                    distance_matrix[indices[a]][indices[b]] ** 2
                    for a in range(len(indices))
                    for b in range(a + 1, len(indices))
                ) / len(indices)

        ss_among = ss_total - ss_within
        df_among = n_populations - 1
        df_within = n_individuals - n_populations
        df_total = n_individuals - 1
        ms_among = ss_among / df_among
        ms_within = ss_within / df_within

        n_correction = (
            n_individuals
            - sum(size**2 for size in population_sizes.values()) / n_individuals
        ) / df_among
        if n_correction <= 0.0:
            raise ValueError("AMOVA sample-size correction is not positive")

        variance_within = ms_within
        variance_among = (ms_among - ms_within) / n_correction
        variance_total = variance_among + variance_within

        if abs(variance_total) <= 1e-15:
            phi_st = 0.0
            percentage_among = 0.0
            percentage_within = 100.0
        else:
            phi_st = variance_among / variance_total
            percentage_among = variance_among / variance_total * 100.0
            percentage_within = variance_within / variance_total * 100.0

        return {
            "ss_total": ss_total,
            "ss_within": ss_within,
            "ss_among": ss_among,
            "df_total": df_total,
            "df_within": df_within,
            "df_among": df_among,
            "ms_within": ms_within,
            "ms_among": ms_among,
            "variance_within": variance_within,
            "variance_among": variance_among,
            "variance_total": variance_total,
            "percentage_within": percentage_within,
            "percentage_among": percentage_among,
            "phi_st": phi_st,
        }

    @staticmethod
    def _pairwise_difference_haplotypes(sample_genotypes, sample_populations):
        """Build 2N-haplotype pairwise-difference matrix (Excoffier 1992).

        Each diploid individual contributes two haplotypes (allele sequences).
        Distance between haplotypes is the Hamming distance — the number of
        loci at which their alleles differ.

        Returns ``(haplotype_ids, haplotype_pop_labels, distance_matrix)``.
        """
        # Collect all loci that appear in at least one call
        loci = set()
        for genotype_dict in sample_genotypes.values():
            loci.update(genotype_dict.keys())
        loci = sorted(loci, key=str)
        if not loci:
            raise ValueError("No loci found in genotype data")

        haplotypes = []  # list of (id, population, allele_vector)
        for sample_id in sorted(sample_genotypes, key=str):
            if sample_id not in sample_populations:
                continue
            pop = sample_populations[sample_id]
            genotype_dict = sample_genotypes[sample_id]

            hap1 = []
            hap2 = []
            for locus in loci:
                gt = genotype_dict.get(locus)
                valid = PopulationGenetics._valid_diploid_genotype(gt)
                if valid is None:
                    hap1.append(None)
                    hap2.append(None)
                else:
                    a1, a2 = valid
                    hap1.append(a1)
                    hap2.append(a2)

            # Skip individuals with no valid loci at all
            if all(a is None for a in hap1):
                continue

            haplotypes.append((f"{sample_id}_1", pop, hap1))
            haplotypes.append((f"{sample_id}_2", pop, hap2))

        n_hap = len(haplotypes)
        if n_hap < 4:
            raise ValueError(
                "At least two complete diploid individuals from two "
                "populations are required"
            )

        distances = [[0.0] * n_hap for _ in range(n_hap)]
        for i in range(n_hap):
            vec_i = haplotypes[i][2]
            for j in range(i + 1, n_hap):
                d = sum(
                    1
                    for a, b in zip(vec_i, haplotypes[j][2])
                    if a is not None and b is not None and a != b
                )
                distances[i][j] = float(d)
                distances[j][i] = float(d)

        hap_ids = [h[0] for h in haplotypes]
        hap_pops = [h[1] for h in haplotypes]

        if len(set(hap_pops)) < 2:
            raise ValueError("At least two populations are required")

        return hap_ids, hap_pops, distances

    @staticmethod
    def _excoffier_amova_components(hap_distances, hap_populations):
        """Excoffier (1992) AMOVA from 2N-haplotype pairwise-difference matrix.

        Returns a dictionary with SS, df, MS, variance components, and F_ST
        for *among*- and *within*-population levels plus the total.
        """
        n_hap = len(hap_populations)
        populations = list(dict.fromkeys(hap_populations))
        n_pop = len(populations)

        pop_indices = {
            p: [i for i, lbl in enumerate(hap_populations) if lbl == p]
            for p in populations
        }
        pop_sizes = {p: len(idx) for p, idx in pop_indices.items()}

        # --- SS(T) = (1 / 2N)  Σ_{i<j}  d_ij^2  (Excoffier et al. 1992 Eq. 1) ---
        ss_total = (
            sum(
                hap_distances[i][j] ** 2
                for i in range(n_hap)
                for j in range(i + 1, n_hap)
            )
            / n_hap
        )

        # --- SS(W) summed over populations ---
        ss_within = 0.0
        for indices in pop_indices.values():
            nk = len(indices)
            if nk > 1:
                ss_within += (
                    sum(
                        hap_distances[indices[a]][indices[b]] ** 2
                        for a in range(nk)
                        for b in range(a + 1, nk)
                    )
                    / nk
                )

        ss_among = ss_total - ss_within

        df_among = n_pop - 1
        df_within = n_hap - n_pop
        df_total = n_hap - 1

        ms_among = ss_among / df_among if df_among > 0 else 0.0
        ms_within = ss_within / df_within if df_within > 0 else 0.0

        n_correction = (
            n_hap - sum(s ** 2 for s in pop_sizes.values()) / n_hap
        ) / df_among if df_among > 0 else 0.0

        if n_correction <= 0.0:
            raise ValueError("AMOVA sample-size correction is not positive")

        var_within = ms_within
        var_among = max(0.0, (ms_among - ms_within) / n_correction)
        var_total = var_among + var_within

        if abs(var_total) <= 1e-15:
            fst = 0.0
            pct_among = 0.0
            pct_within = 100.0
        else:
            fst = var_among / var_total
            pct_among = var_among / var_total * 100.0
            pct_within = var_within / var_total * 100.0

        return {
            "ss_total": ss_total,
            "ss_within": ss_within,
            "ss_among": ss_among,
            "df_total": df_total,
            "df_within": df_within,
            "df_among": df_among,
            "ms_within": ms_within,
            "ms_among": ms_among,
            "variance_within": var_within,
            "variance_among": var_among,
            "variance_total": var_total,
            "percentage_within": pct_within,
            "percentage_among": pct_among,
            "fst": fst,
            "n_haplotypes": n_hap,
        }

    @staticmethod
    def amova(
        sample_genotypes,
        sample_populations=None,
        n_permutations: int = 999,
        random_seed: int = 42,
        min_common_loci: int = 1,
        method: str = "pairwise_difference",
    ):
        """Perform sample-aware one-level AMOVA.

        Parameters
        ----------
        method : str
            ``"pairwise_difference"`` — Excoffier et al. (1992) AMOVA with
            2N-haplotype Hamming distance (compatible with Arlequin 3.5).

            ``"distance_based_multilocus"`` — Smouse & Peakall (1999)
            distance-based AMOVA with N-individual allele-sharing distance.
        """
        if sample_populations is None:
            raise ValueError(
                "AMOVA now requires sample-level genotypes and a separate "
                "sample-to-population mapping; pooled locus observations are "
                "not a valid AMOVA input"
            )
        if n_permutations < 0:
            raise ValueError("n_permutations must be non-negative")
        if method not in ("pairwise_difference", "distance_based_multilocus"):
            raise ValueError(
                f"Unknown AMOVA method: {method!r}; "
                f"use 'pairwise_difference' or 'distance_based_multilocus'"
            )

        # ---- Build distance matrix and run core decomposition ----
        if method == "pairwise_difference":
            # Excoffier (1992): 2N haplotypes × Hamming distance
            hap_ids, hap_labels, hap_dists = (
                PopulationGenetics._pairwise_difference_haplotypes(
                    sample_genotypes, sample_populations
                )
            )
            observed = PopulationGenetics._excoffier_amova_components(
                hap_dists, hap_labels
            )

            # Permutation: shuffle INDIVIDUALS among populations.
            # The two haplotypes of one individual must stay together —
            # shuffling the 2N haplotype labels independently (v1.0.0)
            # breaks the within-individual structure and makes the null
            # distribution too narrow, i.e. p-values come out too small.
            # Individuals occupy consecutive pairs (2k, 2k+1) by construction.
            n_individuals = len(hap_labels) // 2
            if 2 * n_individuals != len(hap_labels):
                raise ValueError("Haplotype list is not paired by individual")
            individual_labels = [hap_labels[2 * k] for k in range(n_individuals)]
            if any(hap_labels[2 * k + 1] != individual_labels[k]
                   for k in range(n_individuals)):
                raise ValueError(
                    "Haplotype pairing broken: an individual's two haplotypes "
                    "carry different population labels")
            rng = _random_mod.Random(random_seed)
            count_extreme = 0
            for _ in range(n_permutations):
                shuffled = list(individual_labels)
                rng.shuffle(shuffled)
                perm_labels: List[str] = []
                for label in shuffled:
                    perm_labels.extend([label, label])
                permuted = PopulationGenetics._excoffier_amova_components(
                    hap_dists, perm_labels
                )
                if permuted["fst"] >= observed["fst"]:
                    count_extreme += 1

            p_value = (count_extreme + 1) / (n_permutations + 1)

            all_loci = set()
            for genotypes in sample_genotypes.values():
                all_loci.update(
                    l for l, g in genotypes.items()
                    if PopulationGenetics._valid_diploid_genotype(g) is not None
                )
            n_individuals = len([
                sid for sid in sample_genotypes
                if sid in sample_populations
            ])
            metadata = {
                "method": method,
                "n_individuals": n_individuals,
                "n_loci": len(all_loci),
                "n_permutations": n_permutations,
                "random_seed": random_seed,
            }

            return [
                AMOVAResult(
                    source="Among populations",
                    df=observed["df_among"],
                    ss=observed["ss_among"],
                    ms=observed["ms_among"],
                    variance_component=observed["variance_among"],
                    percentage=observed["percentage_among"],
                    phi_stat=observed["fst"],
                    p_value=p_value,
                    **metadata,
                ),
                AMOVAResult(
                    source="Within populations",
                    df=observed["df_within"],
                    ss=observed["ss_within"],
                    ms=observed["ms_within"],
                    variance_component=observed["variance_within"],
                    percentage=observed["percentage_within"],
                    **metadata,
                ),
                AMOVAResult(
                    source="Total",
                    df=observed["df_total"],
                    ss=observed["ss_total"],
                    ms=0.0,
                    variance_component=observed["variance_total"],
                    percentage=100.0,
                    **metadata,
                ),
            ]

        else:
            # Smouse & Peakall (1999): N individuals × allele-sharing distance
            sample_ids, labels, distances, _ = (
                PopulationGenetics.multilocus_distance_matrix(
                    sample_genotypes,
                    sample_populations,
                    min_common_loci=min_common_loci,
                )
            )
            observed = PopulationGenetics._amova_components_from_distance(
                distances, labels
            )

            rng = _random_mod.Random(random_seed)
            count_extreme = 0
            for _ in range(n_permutations):
                perm_labels = list(labels)
                rng.shuffle(perm_labels)
                permuted = PopulationGenetics._amova_components_from_distance(
                    distances, perm_labels
                )
                if permuted["phi_st"] >= observed["phi_st"]:
                    count_extreme += 1

            p_value = (count_extreme + 1) / (n_permutations + 1)
            all_loci = {
                locus
                for sample_id in sample_ids
                for locus in sample_genotypes[sample_id]
                if PopulationGenetics._valid_diploid_genotype(
                    sample_genotypes[sample_id][locus]
                )
                is not None
            }
            metadata = {
                "method": method,
                "n_individuals": len(sample_ids),
                "n_loci": len(all_loci),
                "n_permutations": n_permutations,
                "random_seed": random_seed,
            }

            return [
                AMOVAResult(
                    source="Among populations",
                    df=observed["df_among"],
                    ss=observed["ss_among"],
                    ms=observed["ms_among"],
                    variance_component=observed["variance_among"],
                    percentage=observed["percentage_among"],
                    phi_stat=observed["phi_st"],
                    p_value=p_value,
                    **metadata,
                ),
                AMOVAResult(
                    source="Within populations",
                    df=observed["df_within"],
                    ss=observed["ss_within"],
                    ms=observed["ms_within"],
                    variance_component=observed["variance_within"],
                    percentage=observed["percentage_within"],
                    **metadata,
                ),
                AMOVAResult(
                    source="Total",
                    df=observed["df_total"],
                    ss=observed["ss_total"],
                    ms=0.0,
                    variance_component=observed["variance_total"],
                    percentage=100.0,
                    **metadata,
                ),
            ]

    @staticmethod
    def estimate_null_allele_frequency(ho, he):
        if he <= 0:
            return None
        return max(0.0, (he - ho) / (he + ho))

    # ---------- Multiple Testing Correction ----------

    @staticmethod
    def bonferroni_correction(p_values: List[Optional[float]]) -> List[Optional[float]]:
        vi = [i for i, p in enumerate(p_values) if p is not None]
        m = len(vi)
        if m == 0:
            return [None] * len(p_values)
        r: List[Optional[float]] = [None] * len(p_values)
        for i in vi:
            r[i] = min(1.0, p_values[i] * m)  # type: ignore
        return r

    @staticmethod
    def benjamini_hochberg_fdr(
        p_values: List[Optional[float]],
    ) -> List[Optional[float]]:
        valid = [(i, p) for i, p in enumerate(p_values) if p is not None]
        if not valid:
            return [None] * len(p_values)
        m = len(valid)
        sv = sorted(valid, key=lambda x: x[1])
        ranks = {idx: r + 1 for r, (idx, _) in enumerate(sv)}
        qr: Dict[int, float] = {idx: min(1.0, p * m / ranks[idx]) for idx, p in valid}
        pq = 1.0
        for idx, _ in reversed(sv):
            qr[idx] = min(qr[idx], pq)
            pq = qr[idx]
        result: List[Optional[float]] = [None] * len(p_values)
        for idx, q in qr.items():
            result[idx] = q
        return result

    # ---------- Per-locus analysis ----------

    @staticmethod
    def analyze_locus(
        genotypes_by_pop: Dict[str, List[Tuple[float, ...]]],
        locus_name: str = "",
        n_permutations: int = 999,
        random_seed: int = 42,
    ) -> LocusSummary:
        per_pop = {}
        for pop, gt in genotypes_by_pop.items():
            ho, he, nv = PopulationGenetics.compute_heterozygosity(gt)
            if nv == 0:
                continue
            fr, tot = PopulationGenetics.compute_allele_frequencies(gt)
            # compute_heterozygosity already returns Nei's (1978) unbiased
            # gene diversity, He = (2n / (2n - 1)) * (1 - sum(p_i^2)).  The
            # correction must NOT be applied a second time here (v1.0.0 did,
            # which inflated the rarely-used he_unbiased field).
            hu = he
            pic = PopulationGenetics.compute_pic(fr)
            hp, hm = PopulationGenetics.hwe_exact_test(gt)
            fis = 1.0 - ho / he if he > 0 else None
            nf = PopulationGenetics.estimate_null_allele_frequency(ho, he)
            per_pop[pop] = LocusStats(
                locus_name, pop, nv, len(fr), fr, ho, he, hu, pic, hp, hm, fis, nf
            )
        allg = []
        [allg.extend(gt) for gt in genotypes_by_pop.values()]
        oho, ohe, _ = PopulationGenetics.compute_heterozygosity(allg)
        ofr, _ = PopulationGenetics.compute_allele_frequencies(allg)
        opic = PopulationGenetics.compute_pic(ofr)
        ona = len(ofr)
        _, fst, fit = PopulationGenetics.compute_f_statistics(genotypes_by_pop)
        fst_p = None
        if fst is not None and len(genotypes_by_pop) > 1:
            fst, fst_p = PopulationGenetics.fst_permutation_test(
                genotypes_by_pop, n_permutations, random_seed
            )
        return LocusSummary(
            locus_name,
            list(genotypes_by_pop.keys()),
            oho,
            ohe,
            opic,
            ona,
            per_pop,
            fst,
            fit,
            fst_p,
            "WC84_variance_component",
            n_permutations,
            random_seed,
        )

    # ---------- Genotype freq chi-square ----------

    @staticmethod
    def _merge_low_frequency_genotypes(obs, gls):
        """Merge genotype rows with expected count < 5 under independence.

        Returns (merged_obs, merged_gls, merge_map) where merge_map maps
        original genotype labels to the merged label they were collapsed into.
        When a 2x2 table has expected counts < 5 after merging, returns the
        original table unchanged with a flag for Fisher's exact test.
        """
        obs = [list(row) for row in obs]
        gls = list(gls)
        n_pops = len(obs[0]) if obs else 0
        if n_pops == 0:
            return obs, gls, {}

        row_totals = [sum(row) for row in obs]
        col_totals = [sum(obs[r][c] for r in range(len(obs))) for c in range(n_pops)]
        grand_total = sum(row_totals)
        if grand_total == 0:
            return obs, gls, {}

        # Build expected counts
        expected = [
            [row_totals[r] * col_totals[c] / grand_total for c in range(n_pops)]
            for r in range(len(obs))
        ]

        merge_map: Dict[str, str] = {}
        merged_obs = []
        merged_gls = []

        # Collect rows needing merge (expected < 5 in any cell)
        low_rows = []
        low_names = []
        for i in range(len(obs)):
            if any(expected[i][c] < 5 for c in range(n_pops)):
                low_rows.append(obs[i])
                low_names.append(gls[i])
            else:
                merged_obs.append(obs[i])
                merged_gls.append(gls[i])

        if low_rows:
            # Merge all low-frequency rows into an "Others" category
            merged_row = [sum(low_rows[r][c] for r in range(len(low_rows))) for c in range(n_pops)]
            merged_label = "Others (merged: " + ", ".join(low_names) + ")"
            for name in low_names:
                merge_map[name] = merged_label
            merged_obs.append(merged_row)
            merged_gls.append(merged_label)

        # If only 1 row left after merging, test is meaningless
        if len(merged_obs) < 2:
            return obs, gls, merge_map

        return merged_obs, merged_gls, merge_map

    @staticmethod
    def _chi2_statistic(counts) -> float:
        """Pearson chi-square statistic of a contingency table (no Yates)."""
        import numpy as np
        table = np.asarray(counts, dtype=float)
        total = table.sum()
        if total <= 0:
            return 0.0
        row_totals = table.sum(axis=1, keepdims=True)
        col_totals = table.sum(axis=0, keepdims=True)
        expected = row_totals * col_totals / total
        mask = expected > 0
        residuals = (table - expected)[mask]
        return float((residuals ** 2 / expected[mask]).sum())

    @staticmethod
    def _permutation_chi2_p(labels: List[str], populations: List[str],
                            n_permutations: int = 999,
                            random_seed: int = 42) -> Optional[float]:
        """Permutation p-value for genotype-count differentiation.

        Population labels are shuffled among individuals (genotype classes stay
        attached to their individual), the contingency table is rebuilt and the
        chi-square statistic recomputed; p = (count + 1) / (n + 1).  This gives
        a valid test for R x C tables whose expected counts are too small for
        the asymptotic chi-square approximation.
        """
        if len(labels) != len(populations) or not labels:
            return None
        classes = sorted(set(labels))
        pops = sorted(set(populations))
        class_index = {c: i for i, c in enumerate(classes)}
        pop_index = {p: i for i, p in enumerate(pops)}

        def build(label_list, pop_list):
            table = [[0] * len(pops) for _ in classes]
            for label, pop in zip(label_list, pop_list):
                table[class_index[label]][pop_index[pop]] += 1
            return table

        observed = PopulationGenetics._chi2_statistic(build(labels, populations))
        rng = _random_mod.Random(random_seed)
        count_extreme = 0
        for _ in range(n_permutations):
            shuffled = list(populations)
            rng.shuffle(shuffled)
            statistic = PopulationGenetics._chi2_statistic(
                build(labels, shuffled))
            if statistic >= observed - 1e-12:
                count_extreme += 1
        return (count_extreme + 1) / (n_permutations + 1)

    @staticmethod
    def test_locus_genotype_freq(
        genotypes_by_pop: Dict[str, List[Tuple[float, ...]]], locus_name: str = "",
        n_permutations: int = 999,
        random_seed: int = 42,
    ) -> Optional[GenotypeFreqTestResult]:
        pops = list(genotypes_by_pop.keys())
        if len(pops) < 2:
            return None

        def gl(gt):
            if len(gt) == 0:
                return "missing"
            if len(gt) == 1:
                return str(gt[0])
            # Sort the two alleles so that "T/A" and "A/T" land in the same
            # class regardless of the order in the source matrix.
            return "/".join(sorted([str(gt[0]), str(gt[1])]))

        als = set()
        [
            als.add(gl(gt))
            for gt_lst in genotypes_by_pop.values()
            for gt in gt_lst
            if len(gt) > 0
        ]
        gls = sorted(als)
        if len(gls) < 2:
            return GenotypeFreqTestResult(
                locus_name, 0.0, 0, 1.0, 1.0,
                None, None, False, False, False, False,
                gls, pops, [],
                warning="Only one genotype class",
            )

        # Build original contingency table
        obs_orig = []
        for glb in gls:
            row = [
                sum(1 for gt in genotypes_by_pop.get(p, []) if len(gt) > 0 and gl(gt) == glb)
                for p in pops
            ]
            obs_orig.append(row)

        # Compute frequency table
        ft = []
        for row in obs_orig:
            tr = sum(row)
            ft.append([v / tr for v in row] if tr > 0 else [0.0] * len(row))

        # Merge low-frequency genotypes
        obs_merged, gls_merged, merge_map = PopulationGenetics._merge_low_frequency_genotypes(
            obs_orig, gls
        )

        warnings = []
        if merge_map:
            warnings.append(f"Merged {len(merge_map)} low-frequency genotype(s): {'; '.join(f'{k}->{v}' for k, v in merge_map.items())}")

        try:
            from scipy.stats import chi2_contingency
            from scipy.stats import fisher_exact as _fisher_exact

            # Use merged table for chi-square if merging occurred
            test_table = obs_merged if merge_map else obs_orig
            n_rows = len(test_table)
            n_cols = len(test_table[0]) if test_table else 0

            # Check if we should use Fisher's exact test (2x2 with small expected counts)
            use_fisher = False
            if n_rows == 2 and n_cols == 2:
                # Compute expected counts for the test table
                row_sums = [sum(r) for r in test_table]
                col_sums = [sum(test_table[r][c] for r in range(n_rows)) for c in range(n_cols)]
                grand = sum(row_sums)
                if grand > 0:
                    exp_table = [
                        [row_sums[r] * col_sums[c] / grand for c in range(n_cols)]
                        for r in range(n_rows)
                    ]
                    if any(exp_table[r][c] < 5 for r in range(n_rows) for c in range(n_cols)):
                        use_fisher = True

            if use_fisher:
                # Fisher's exact test for 2x2 table
                odds_ratio, pv = _fisher_exact(test_table)
                dof = 1
                chi2 = 0.0
                test_method = "fisher_exact"
                warnings.append("Fisher's exact test used (expected counts < 5 in 2x2 table)")
            else:
                chi2, pv, dof, exp = chi2_contingency(test_table)
                test_method = "chi2"
                # Check remaining low expected counts
                ef = exp.flatten()
                nl = sum(1 for e in ef if e < 5)
                if nl > 0 and not use_fisher:
                    # R x C table with small expected counts: the chi-square
                    # approximation is unreliable with the small per-population
                    # samples typical of Sanger validation panels.  Replace the
                    # asymptotic p-value with a permutation p-value (population
                    # labels shuffled among individuals, genotype classes kept).
                    merged_labels = []
                    for pop, gt_list in genotypes_by_pop.items():
                        for gt in gt_list:
                            if len(gt) == 0:
                                continue
                            label = gl(gt)
                            merged_labels.append(merge_map.get(label, label))
                    perm_p = PopulationGenetics._permutation_chi2_p(
                        merged_labels,
                        [p for p, gt_list in genotypes_by_pop.items()
                         for gt in gt_list if len(gt) > 0],
                        n_permutations=n_permutations,
                        random_seed=random_seed,
                    )
                    if perm_p is not None:
                        pv = perm_p
                        test_method = "chi2_permutation"
                        warnings.append(
                            f"{nl} cell(s) have expected frequency < 5 after "
                            f"merging — asymptotic chi-square p replaced by a "
                            f"permutation p-value ({n_permutations} shuffles)")
                    else:
                        warnings.append(
                            f"{nl} cell(s) still have expected frequency < 5 "
                            "after merging — interpret with caution")
                elif nl > 0:
                    warnings.append(
                        f"{nl} cell(s) still have expected frequency < 5 after "
                        "merging — interpret with caution")
        except (ImportError, ValueError) as exc:
            return GenotypeFreqTestResult(
                locus_name, 0.0, 0, None, None,
                None, None, False, False, False, False,
                gls, pops, ft,
                warning=f"unavailable: {exc}",
            )

        return GenotypeFreqTestResult(
            locus_name, chi2, dof, pv, pv,
            None, None,
            bool(pv is not None and pv < 0.05), False, False, False,
            gls, pops, ft,
            original_contingency_table=obs_orig,
            merged_contingency_table=obs_merged if merge_map else None,
            merged_genotype_map=merge_map if merge_map else None,
            test_method=test_method,
            warning="; ".join(warnings) if warnings else None,
        )

    @staticmethod
    def test_all_loci_genotype_freq(
        all_data: Dict[str, Dict[str, List[Tuple[float, ...]]]], alpha: float = 0.05
    ) -> Dict[str, Optional[GenotypeFreqTestResult]]:
        results = {}
        for ln, gbp in all_data.items():
            results[ln] = PopulationGenetics.test_locus_genotype_freq(gbp, ln)
        lns = list(all_data.keys())
        rp = [
            results[n].raw_p_value if results.get(n) is not None else None for n in lns
        ]
        bonf = PopulationGenetics.bonferroni_correction(rp)
        fdrv = PopulationGenetics.benjamini_hochberg_fdr(rp)
        for i, n in enumerate(lns):
            r = results.get(n)
            if r is None:
                continue
            r.bonferroni_p_value = bonf[i]
            r.fdr_q_value = fdrv[i]
            # Explicit bool() to avoid numpy.bool_ issues
            r.significant_raw = bool(
                r.raw_p_value is not None and r.raw_p_value < alpha
            )
            r.significant_bonferroni = bool(bonf[i] is not None and bonf[i] < alpha)
            r.significant_fdr = bool(fdrv[i] is not None and fdrv[i] < alpha)
            r.significant = r.significant_fdr
            r.p_value = r.raw_p_value
        return results
