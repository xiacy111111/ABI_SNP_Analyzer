# -*- coding: utf-8 -*-
"""Linkage disequilibrium analysis and independent-locus selection.

SNPs called from the same Sanger amplicon are physically linked and usually in
strong linkage disequilibrium (LD).  Treating them as independent loci inflates
the degrees of freedom of multi-locus tests, underestimates the variance and
violates the independent-loci assumption of software such as STRUCTURE
(Excoffier et al. 1992; Falush et al. 2003).  This module provides the two
remedies requested for ABI-SNP Analyzer:

1. **LD report** — pairwise r² and D' for every SNP pair, computed *within*
   each population (pooling diverged populations would inflate LD through the
   Wahlburg effect), plus a Monte-Carlo permutation p-value.
2. **Independent-locus selection** — SNPs are clustered into linkage blocks
   (single linkage on r² >= threshold) and one representative SNP per block is
   returned, so that the population-genetics module can be run in a
   "one SNP per amplicon block" mode alongside the "all SNPs" mode.

Estimators
----------
r²
   Squared Pearson correlation of the minor-allele dosages (0/1/2) over
   individuals genotyped at both loci — the standard unphased estimator used by
   PLINK ``--r2`` (Hill & Weir 1994; Purcell et al. 2007).
D'
   Lewontin's normalized D (Lewontin 1964) computed from two-locus haplotype
   frequencies estimated with the expectation-maximization algorithm for
   unphased genotype data (Excoffier & Slatkin 1995).
p
   Permutation test on the dosage r² (locus B genotypes shuffled among
   individuals; the observed r² is compared against the permutation
   distribution, ``p = (1 + #{r2_perm >= r2_obs}) / (n_perm + 1)``).

References
----------
Lewontin, R.C. (1964) The interaction of selection and linkage. I. General
considerations; heterotic models. *Genetics* 49: 49-67.
Excoffier, L. & Slatkin, M. (1995) Maximum-likelihood estimation of molecular
haplotype frequencies in a diploid population. *Mol Biol Evol* 12: 921-927.
Hill, W.G. & Weir, B.S. (1994) Maximum-likelihood estimation of gene location
by linkage disequilibrium. *Theor Appl Genet* 89: 840-844.
Purcell, S. et al. (2007) PLINK: a tool set for whole-genome association and
population-based linkage analyses. *Am J Hum Genet* 81: 559-575.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "LdPair",
    "LinkageBlock",
    "LdSummary",
    "compute_pairwise_ld",
    "summarise_ld_pairs",
    "cluster_linkage_blocks",
    "choose_representatives",
    "select_independent_loci",
    "write_linkage_report",
    "DEFAULT_R2_THRESHOLD",
    "DEFAULT_N_PERMUTATIONS",
]

DEFAULT_R2_THRESHOLD = 0.80
DEFAULT_N_PERMUTATIONS = 999
DEFAULT_RANDOM_SEED = 42
#: Minimum individuals genotyped at both loci for a pairwise LD estimate.
MIN_PAIR_N = 5
EM_MAX_ITER = 300
EM_TOLERANCE = 1e-12


@dataclass
class LdPair:
    """Pairwise LD between two SNPs, estimated within one population."""

    locus1: str
    locus2: str
    population: str
    n_individuals: int
    r2: Optional[float] = None          # unphased dosage estimator
    d_prime: Optional[float] = None     # EM haplotype estimator
    p_value: Optional[float] = None
    note: str = ""


@dataclass
class LdSummary:
    """Cross-population summary for one SNP pair."""

    locus1: str
    locus2: str
    n_populations: int
    n_individuals: int
    r2_max: Optional[float] = None
    r2_min: Optional[float] = None
    d_prime_max: Optional[float] = None
    p_value_min: Optional[float] = None
    populations_compared: List[str] = field(default_factory=list)
    note: str = ""


@dataclass
class LinkageBlock:
    """A set of SNPs that cannot be treated as independent loci."""

    block_id: int
    loci: List[str]
    representative: Optional[str] = None
    representative_reason: str = ""


# ---------------------------------------------------------------------------
# genotype encoding helpers
# ---------------------------------------------------------------------------

def _locus_codes(sample_genotypes: Dict[str, Dict[str, object]],
                 sample_order: Sequence[str], locus: str
                 ) -> Tuple[np.ndarray, Dict[str, int]]:
    """Encode one locus as 0 = homo major, 1 = het, 2 = homo minor, -1 missing.

    Returns ``(codes, minor_allele)``; ``minor_allele`` maps ``{"major": x,
    "minor": y}`` onto the base letters actually used.
    """
    from collections import Counter

    allele_counter: Counter = Counter()
    raw: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for sample in sample_order:
        call = sample_genotypes.get(sample, {}).get(locus)
        alleles: Tuple[Optional[str], Optional[str]] = (None, None)
        if call is not None:
            if isinstance(call, (tuple, list)):
                a1, a2 = call[0], call[1]
            else:
                text = str(call).strip().upper()
                a1, a2 = (text[0], text[1]) if len(text) >= 2 else (None, None)
            a1 = a1 if a1 in ("A", "C", "G", "T") else None
            a2 = a2 if a2 in ("A", "C", "G", "T") else None
            alleles = (a1, a2)
            for allele in alleles:
                if allele:
                    allele_counter[allele] += 1
        raw[sample] = alleles

    if not allele_counter:
        codes = np.full(len(sample_order), -1, dtype=np.int8)
        return codes, {"major": "", "minor": ""}

    ordered = [a for a, _ in allele_counter.most_common()]
    if len(ordered) == 1:
        # Monomorphic locus: every called individual is homozygous for the
        # only allele.  Coding it as 0 (rather than as missing) lets the
        # caller report the pair with an explicit "no variation" note instead
        # of silently dropping it.
        codes = np.full(len(sample_order), -1, dtype=np.int8)
        for index, sample in enumerate(sample_order):
            if raw[sample][0] is not None:
                codes[index] = 0
        return codes, {"major": ordered[0], "minor": ""}

    major, minor = ordered[0], ordered[1]
    codes = np.empty(len(sample_order), dtype=np.int8)
    for index, sample in enumerate(sample_order):
        a1, a2 = raw[sample]
        if a1 is None or a2 is None:
            codes[index] = -1
            continue
        if (a1, a2) not in ((major, major), (major, minor), (minor, major),
                            (minor, minor)):
            codes[index] = -1        # third allele -> treat as missing
            continue
        minor_count = sum(1 for a in (a1, a2) if a == minor)
        codes[index] = minor_count   # 0, 1 or 2
    return codes, {"major": major, "minor": minor}


# ---------------------------------------------------------------------------
# pairwise LD
# ---------------------------------------------------------------------------

def _em_haplotype_frequencies(codes_a: np.ndarray, codes_b: np.ndarray
                              ) -> Optional[Tuple[np.ndarray, int]]:
    """EM estimate of the four two-locus haplotype frequencies.

    ``codes_a``/``codes_b`` are 0/1/2 minor-allele counts (missing already
    removed).  Haplotypes are ordered ``(A0B0, A0B1, A1B0, A1B1)`` where
    ``A1``/``B1`` are the minor alleles.
    """
    n = codes_a.size
    if n < MIN_PAIR_N:
        return None
    p_a1 = float(codes_a.mean()) / 2.0
    p_b1 = float(codes_b.mean()) / 2.0
    if not (0.0 < p_a1 < 1.0) or not (0.0 < p_b1 < 1.0):
        return None

    x = np.array([(1 - p_a1) * (1 - p_b1), (1 - p_a1) * p_b1,
                  p_a1 * (1 - p_b1), p_a1 * p_b1], dtype=float)
    double_het = (codes_a == 1) & (codes_b == 1)
    n_double_het = int(double_het.sum())

    determined_a = codes_a[~double_het]
    determined_b = codes_b[~double_het]

    for _ in range(EM_MAX_ITER):
        counts = np.zeros(4, dtype=float)
        if determined_a.size:
            # per individual the two haplotypes are (major/minor split):
            a1_idx = np.where(determined_a == 2, 1, 0)
            a2_idx = np.where(determined_a >= 1, 1, 0)
            b1_idx = np.where(determined_b == 2, 1, 0)
            b2_idx = np.where(determined_b >= 1, 1, 0)
            h1 = a1_idx * 2 + b1_idx
            h2 = a2_idx * 2 + b2_idx
            counts += np.bincount(h1, minlength=4)
            counts += np.bincount(h2, minlength=4)
        if n_double_het:
            w = (x[0] + x[3]) / max(x[0] + x[3] + x[1] + x[2], 1e-300)
            counts[0] += n_double_het * w
            counts[3] += n_double_het * w
            counts[1] += n_double_het * (1.0 - w)
            counts[2] += n_double_het * (1.0 - w)
        total = counts.sum()
        if total <= 0:
            return None
        x_new = counts / total
        if np.max(np.abs(x_new - x)) < EM_TOLERANCE:
            x = x_new
            break
        x = x_new
    return x, n


def _ld_from_haplotypes(x: np.ndarray) -> Tuple[float, float]:
    """Return ``(r2, d_prime)`` from haplotype frequencies."""
    p_a1 = x[2] + x[3]
    p_b1 = x[1] + x[3]
    d = x[0] * x[3] - x[1] * x[2]
    if d >= 0:
        d_max = min((1 - p_a1) * (1 - p_b1), p_a1 * p_b1)
    else:
        d_max = -min((1 - p_a1) * p_b1, p_a1 * (1 - p_b1))
    d_prime = d / d_max if d_max else None
    denom = p_a1 * (1 - p_a1) * p_b1 * (1 - p_b1)
    r2 = (d * d) / denom if denom > 0 else None
    return (r2 if r2 is not None else float("nan"),
            d_prime if d_prime is not None else float("nan"))


def _dosage_r2(dosage_a: np.ndarray, dosage_b: np.ndarray) -> float:
    """Squared Pearson correlation of minor-allele dosages."""
    if dosage_a.size < MIN_PAIR_N:
        return float("nan")
    sd_a = dosage_a.std()
    sd_b = dosage_b.std()
    if sd_a <= 0 or sd_b <= 0:
        return float("nan")
    cov = float(((dosage_a - dosage_a.mean()) *
                 (dosage_b - dosage_b.mean())).mean())
    r = cov / (sd_a * sd_b)
    return float(max(-1.0, min(1.0, r)) ** 2)


def _permutation_p(dosage_a: np.ndarray, dosage_b: np.ndarray,
                   observed_r2: float, n_permutations: int,
                   rng: np.random.Generator) -> Optional[float]:
    """Monte-Carlo p-value for the dosage r² (shuffling locus B genotypes)."""
    if not np.isfinite(observed_r2) or n_permutations <= 0:
        return None
    m = dosage_a.size
    if m < MIN_PAIR_N:
        return None
    if n_permutations == 1:
        return None
    indices = np.tile(np.arange(m), (n_permutations, 1))
    rng.permuted(indices, axis=1, out=indices)
    permuted = dosage_b[indices]                       # (P, m)
    za = (dosage_a - dosage_a.mean()) / dosage_a.std()
    zb = (permuted - permuted.mean(axis=1, keepdims=True)) / \
        np.maximum(permuted.std(axis=1, keepdims=True), 1e-300)
    r_perm = (zb @ za) / m
    exceed = int(np.sum(r_perm ** 2 >= observed_r2 - 1e-12))
    return (1.0 + exceed) / (1.0 + n_permutations)


def compute_pairwise_ld(sample_genotypes: Dict[str, Dict[str, object]],
                        loci: Sequence[str],
                        sample_to_population: Dict[str, str],
                        n_permutations: int = DEFAULT_N_PERMUTATIONS,
                        random_seed: int = DEFAULT_RANDOM_SEED,
                        progress_cb: Optional[Callable[[float, str], None]] = None
                        ) -> List[LdPair]:
    """Pairwise LD for every SNP pair, estimated within each population.

    ``sample_genotypes`` maps ``sample -> {locus: (allele1, allele2) | None}``
    (the layout produced by :func:`build_population_analysis_inputs`).  Pooling
    populations is deliberately avoided: divergence between populations
    inflates LD (Wahlburg effect), so every pair is estimated per population
    and summarised across populations by :func:`summarise_ld_pairs`.
    """
    sample_order = sorted(sample_genotypes)
    populations: Dict[str, List[str]] = {}
    for sample in sample_order:
        populations.setdefault(sample_to_population.get(sample, "Unassigned"),
                               []).append(sample)

    codes = {locus: _locus_codes(sample_genotypes, sample_order, locus)
             for locus in loci}

    pairs: List[LdPair] = []
    total_pairs = max(len(loci) * (len(loci) - 1) // 2, 1)
    completed = 0
    rng = np.random.default_rng(random_seed)

    for i in range(len(loci)):
        for j in range(i + 1, len(loci)):
            locus1, locus2 = loci[i], loci[j]
            codes1, _alleles1 = codes[locus1]
            codes2, _alleles2 = codes[locus2]

            per_population: List[LdPair] = []
            for pop_name in sorted(populations):
                members = populations[pop_name]
                idx = [sample_order.index(s) for s in members]
                ca = codes1[idx]
                cb = codes2[idx]
                valid = (ca >= 0) & (cb >= 0)
                n_valid = int(valid.sum())
                if n_valid < MIN_PAIR_N:
                    continue
                da = ca[valid].astype(float)
                db = cb[valid].astype(float)

                r2 = _dosage_r2(da, db)
                note = ""
                em = _em_haplotype_frequencies(ca[valid], cb[valid])
                if em is None:
                    d_prime = float("nan")
                    note = "monomorphic or too few informative individuals"
                else:
                    _r2_hap, d_prime = _ld_from_haplotypes(em[0])
                p_value = _permutation_p(da, db, r2, n_permutations, rng) \
                    if np.isfinite(r2) else None
                per_population.append(LdPair(
                    locus1=locus1, locus2=locus2, population=pop_name,
                    n_individuals=n_valid, r2=r2, d_prime=d_prime,
                    p_value=p_value, note=note))

            pairs.extend(per_population)
            completed += 1
            if progress_cb and completed % 10 == 0:
                progress_cb(100.0 * completed / total_pairs,
                            f"LD: {completed}/{total_pairs} SNP pairs")
    return pairs


def summarise_ld_pairs(pairs: Sequence[LdPair]) -> List[LdSummary]:
    """Collapse per-population pairs into one row per SNP pair."""
    grouped: Dict[Tuple[str, str], List[LdPair]] = {}
    for pair in pairs:
        key = (pair.locus1, pair.locus2)
        grouped.setdefault(key, []).append(pair)

    summaries: List[LdSummary] = []
    for (locus1, locus2), members in grouped.items():
        r2_values = [p.r2 for p in members if p.r2 is not None and
                     np.isfinite(p.r2)]
        dp_values = [p.d_prime for p in members if p.d_prime is not None and
                     np.isfinite(p.d_prime)]
        p_values = [p.p_value for p in members if p.p_value is not None]
        summaries.append(LdSummary(
            locus1=locus1, locus2=locus2,
            n_populations=len(members),
            n_individuals=sum(p.n_individuals for p in members),
            r2_max=max(r2_values) if r2_values else None,
            r2_min=min(r2_values) if r2_values else None,
            d_prime_max=max(dp_values) if dp_values else None,
            p_value_min=min(p_values) if p_values else None,
            populations_compared=[p.population for p in members],
            note=members[0].note if len(members) == 1 else "",
        ))
    return summaries


# ---------------------------------------------------------------------------
# linkage blocks and representative selection
# ---------------------------------------------------------------------------

def cluster_linkage_blocks(loci: Sequence[str],
                           summaries: Sequence[LdSummary],
                           r2_threshold: float = DEFAULT_R2_THRESHOLD
                           ) -> List[LinkageBlock]:
    """Single-linkage clustering of SNPs on ``r2_max >= r2_threshold``."""
    parent: Dict[str, str] = {locus: locus for locus in loci}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for summary in summaries:
        if summary.r2_max is None or not np.isfinite(summary.r2_max):
            continue
        if summary.r2_max >= r2_threshold:
            union(summary.locus1, summary.locus2)

    groups: Dict[str, List[str]] = {}
    for locus in loci:                      # keep the input order
        groups.setdefault(find(locus), []).append(locus)

    blocks: List[LinkageBlock] = []
    multi = [members for members in groups.values() if len(members) > 1]
    singles = [members for members in groups.values() if len(members) == 1]
    block_id = 1
    for members in multi + singles:
        blocks.append(LinkageBlock(block_id=block_id, loci=list(members)))
        block_id += 1
    return blocks


def choose_representatives(sample_genotypes: Dict[str, Dict[str, object]],
                           sample_order: Sequence[str],
                           blocks: Sequence[LinkageBlock]
                           ) -> None:
    """Pick one representative SNP per block, in place.

    Ranking: fewest missing calls first, then highest minor-allele frequency,
    then the position in the original locus order.  Missing rate and MAF are
    what drive the reliability of the downstream statistics, so both beat the
    arbitrary "first locus" rule.
    """
    stats: Dict[str, Tuple[float, float]] = {}
    for locus in {locus for block in blocks for locus in block.loci}:
        codes, _alleles = _locus_codes(sample_genotypes, sample_order, locus)
        valid = codes >= 0
        n_valid = int(valid.sum())
        missing_rate = 1.0 - (n_valid / len(sample_order)) if sample_order else 1.0
        maf = float(codes[valid].mean()) / 2.0 if n_valid else 0.0
        maf = min(maf, 1.0 - maf)
        stats[locus] = (missing_rate, maf)

    for block in blocks:
        if len(block.loci) == 1:
            block.representative = block.loci[0]
            block.representative_reason = "only member of the block"
            continue
        best = min(
            block.loci,
            key=lambda locus: (round(stats[locus][0], 6),
                               -round(stats[locus][1], 6),
                               block.loci.index(locus)))
        block.representative = best
        block.representative_reason = (
            f"lowest missing rate ({stats[best][0]:.1%}), "
            f"highest MAF ({stats[best][1]:.3f})")


def select_independent_loci(
        sample_genotypes: Dict[str, Dict[str, object]],
        loci: Sequence[str],
        sample_to_population: Dict[str, str],
        r2_threshold: float = DEFAULT_R2_THRESHOLD,
        n_permutations: int = DEFAULT_N_PERMUTATIONS,
        random_seed: int = DEFAULT_RANDOM_SEED,
        progress_cb: Optional[Callable[[float, str], None]] = None,
) -> Tuple[List[str], List[LinkageBlock], List[LdPair], List[LdSummary]]:
    """Full pipeline: LD -> blocks -> one representative SNP per block.

    Returns ``(kept_loci, blocks, pairs, summaries)``.  ``kept_loci`` preserves
    the original locus order and is what the population-genetics module should
    be fed in "one SNP per linkage block" mode.
    """
    pairs = compute_pairwise_ld(
        sample_genotypes, loci, sample_to_population,
        n_permutations=n_permutations, random_seed=random_seed,
        progress_cb=progress_cb)
    summaries = summarise_ld_pairs(pairs)
    blocks = cluster_linkage_blocks(loci, summaries, r2_threshold)
    sample_order = sorted(sample_genotypes)
    choose_representatives(sample_genotypes, sample_order, blocks)
    kept = [block.representative for block in blocks
            if block.representative]
    kept = [locus for locus in loci if locus in kept]
    return kept, blocks, pairs, summaries


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

_METHOD_NOTE = (
    "r2: squared Pearson correlation of minor-allele dosages (0/1/2) over "
    "individuals genotyped at both loci — the unphased estimator used by "
    "PLINK --r2 (Hill & Weir 1994; Purcell et al. 2007). "
    "D': Lewontin's normalised D (Lewontin 1964) from two-locus haplotype "
    "frequencies estimated by EM for unphased genotypes (Excoffier & Slatkin "
    "1995). p: permutation test on r2 (locus-B genotypes shuffled among "
    "individuals; p = (1 + #{r2_perm >= r2_obs}) / (n_perm + 1)). "
    "All estimates are computed WITHIN each population: pooling diverged "
    "populations inflates LD through the Wahlburg effect. "
    "D' is strongly upward-biased in small samples, so r2 is the statistic to "
    "quote for block definitions.")


def write_linkage_report(path: str, pairs: Sequence[LdPair],
                         summaries: Sequence[LdSummary],
                         blocks: Sequence[LinkageBlock],
                         kept_loci: Sequence[str],
                         r2_threshold: float = DEFAULT_R2_THRESHOLD,
                         n_permutations: int = DEFAULT_N_PERMUTATIONS,
                         random_seed: int = DEFAULT_RANDOM_SEED) -> str:
    """Write ``linkage_report.xlsx`` (Pairwise_LD, Linkage_Blocks, Method)."""
    from openpyxl import Workbook
    from openpyxl.styles import Font

    workbook = Workbook()

    sheet = workbook.active
    sheet.title = "Pairwise_LD"
    headers = ["Locus1", "Locus2", "Population", "N", "r2", "D'", "p_perm",
               "Note"]
    sheet.append(headers)
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    for pair in pairs:
        sheet.append([
            pair.locus1, pair.locus2, pair.population, pair.n_individuals,
            _round_or_blank(pair.r2), _round_or_blank(pair.d_prime),
            _round_or_blank(pair.p_value), pair.note])
    _autosize(sheet)

    sheet = workbook.create_sheet("LD_Summary")
    sheet.append(["Locus1", "Locus2", "N_populations", "N_individuals",
                  "r2_max", "r2_min", "D'_max", "p_min", "Populations"])
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    for summary in summaries:
        sheet.append([
            summary.locus1, summary.locus2, summary.n_populations,
            summary.n_individuals, _round_or_blank(summary.r2_max),
            _round_or_blank(summary.r2_min), _round_or_blank(summary.d_prime_max),
            _round_or_blank(summary.p_value_min),
            ", ".join(summary.populations_compared)])
    _autosize(sheet)

    sheet = workbook.create_sheet("Linkage_Blocks")
    sheet.append(["Block", "Loci", "N_loci", "Representative", "Reason",
                  "Kept_for_analysis"])
    for cell in sheet[1]:
        cell.font = Font(bold=True)
    for block in blocks:
        sheet.append([
            f"Block{block.block_id}", ", ".join(block.loci), len(block.loci),
            block.representative or "", block.representative_reason,
            "yes" if block.representative in kept_loci else "no"])
    _autosize(sheet)

    sheet = workbook.create_sheet("Method")
    sheet.column_dimensions["A"].width = 120
    sheet.append(["Linkage disequilibrium report — method notes"])
    sheet["A1"].font = Font(bold=True)
    sheet.append([f"r2 threshold for block definition: {r2_threshold}"])
    sheet.append([f"Permutations per pair: {n_permutations}"])
    sheet.append([f"Random seed: {random_seed}"])
    sheet.append([f"Minimum individuals per pair/population: {MIN_PAIR_N}"])
    sheet.append(["Kept (independent) loci:", ", ".join(kept_loci)])
    sheet.append([_METHOD_NOTE])
    sheet.append([""])
    sheet.append(["References"])
    for reference in (
            "Lewontin, R.C. (1964) Genetics 49:49-67.",
            "Excoffier, L. & Slatkin, M. (1995) Mol Biol Evol 12:921-927.",
            "Hill, W.G. & Weir, B.S. (1994) Theor Appl Genet 89:840-844.",
            "Purcell, S. et al. (2007) Am J Hum Genet 81:559-575 (PLINK)."):
        sheet.append([reference])

    workbook.save(path)
    return path


def _round_or_blank(value: Optional[float], digits: int = 4):
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return round(float(value), digits)


def _autosize(sheet) -> None:
    for column_cells in sheet.columns:
        length = max(len(str(cell.value or "")) for cell in column_cells)
        letter = column_cells[0].column_letter
        sheet.column_dimensions[letter].width = min(max(length + 2, 10), 60)
