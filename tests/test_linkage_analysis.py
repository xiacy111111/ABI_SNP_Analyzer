# -*- coding: utf-8 -*-
"""Tests for linkage_analysis.py (LD estimation + block selection)."""

import os
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import linkage_analysis as la  # noqa: E402

LOCI = ["SNP_A", "SNP_B", "SNP_C"]
N = 80
MAJOR_A, MINOR_A = "A", "T"     # SNP_A
MAJOR_B, MINOR_B = "C", "G"     # SNP_B
MAJOR_C, MINOR_C = "G", "T"     # SNP_C


def _genotypes_from_dosage(major, minor, dosage, rng):
    """Turn a minor-allele dosage (0/1/2) into a diploid call."""
    calls = {0: (major, major), 1: (major, minor), 2: (minor, minor)}
    if dosage is None:
        return None
    return calls[int(dosage)]


def _make_population(population, ab_haplotypes, c_minor_freq, seed,
                     c_minor_freq_b=None, n=N):
    """Build one population.

    ``ab_haplotypes`` fixes the SNP_A/SNP_B haplotypes ([(a, b), ...]);
    SNP_C is drawn independently with the requested minor-allele frequency
    (``c_minor_freq_b`` for population PopB when two populations are pooled).
    """
    rng = np.random.default_rng(seed)
    samples = {}
    for i in range(n):
        h1, h2 = ab_haplotypes[int(rng.integers(0, len(ab_haplotypes)))], \
            ab_haplotypes[int(rng.integers(0, len(ab_haplotypes)))]
        p_c = c_minor_freq_b if (c_minor_freq_b is not None
                                 and population == "PopB") else c_minor_freq
        dosage_c = int(rng.choice([0, 1, 2], p=[(1 - p_c) ** 2,
                                                2 * p_c * (1 - p_c),
                                                p_c ** 2]))
        samples[f"{population}_Ind{i:03d}"] = {
            "SNP_A": (h1[0], h2[0]),
            "SNP_B": (h1[1], h2[1]),
            "SNP_C": _genotypes_from_dosage(MAJOR_C, MINOR_C, dosage_c, rng),
        }
    return samples


def _assemble(populations):
    genotypes, mapping = {}, {}
    for name, samples in populations.items():
        for sample, calls in samples.items():
            genotypes[sample] = calls
            mapping[sample] = name
    return genotypes, mapping


# perfect LD between A and B (haplotypes A-C / T-G), C independent
_link_haps = [("A", "C"), ("T", "G")]
P1 = _make_population("Pop1", _link_haps, c_minor_freq=0.5, seed=11)
genotypes, mapping = _assemble({"Pop1": P1})

_PAIRS = la.compute_pairwise_ld(genotypes, LOCI, mapping,
                                n_permutations=999, random_seed=1)


def _key(locus1, locus2, population=None):
    for pair in _PAIRS:
        if {pair.locus1, pair.locus2} == {locus1, locus2} and (
                population is None or pair.population == population):
            return pair
    raise AssertionError(f"pair {locus1}/{locus2} not found")


def test_linked_pair_has_maximal_ld():
    pair = _key("SNP_A", "SNP_B")
    assert pair.n_individuals == N
    assert pair.r2 is not None and pair.r2 > 0.95, pair.r2
    assert pair.d_prime is not None and abs(pair.d_prime) > 0.99
    assert pair.p_value is not None and pair.p_value < 0.05


def test_independent_pair_is_not_flagged():
    pair = _key("SNP_A", "SNP_C")
    assert pair.r2 is not None and pair.r2 < 0.25, pair.r2
    assert pair.p_value is None or pair.p_value > 0.05, pair.p_value


def test_blocks_split_linked_from_independent():
    summaries = la.summarise_ld_pairs(_PAIRS)
    blocks = la.cluster_linkage_blocks(LOCI, summaries, r2_threshold=0.8)
    la.choose_representatives(genotypes, sorted(genotypes), blocks)
    multi = [b for b in blocks if len(b.loci) > 1]
    assert len(multi) == 1
    assert set(multi[0].loci) == {"SNP_A", "SNP_B"}
    assert multi[0].representative in ("SNP_A", "SNP_B")
    kept = [b.representative for b in blocks]
    assert set(kept) == {"SNP_A", "SNP_C"} or set(kept) == {"SNP_B", "SNP_C"}


def test_select_independent_loci_end_to_end():
    kept, blocks, pairs, summaries = la.select_independent_loci(
        genotypes, LOCI, mapping, r2_threshold=0.8)
    assert len(kept) == 2
    assert "SNP_C" in kept
    assert pairs and summaries
    assert all(b.representative for b in blocks)


# --------------------------------------------------------------------------
# 2. Wahlburg effect: pooling diverged populations manufactures LD
# --------------------------------------------------------------------------

# each population is in linkage equilibrium for A-C, but the minor-allele
# frequency of SNP_A/SNP_C differs between populations (0.1 vs 0.9), so
# pooling them manufactures a spurious correlation.
P_A = _make_population("PopA", [("A", "C")] * 9 + [("T", "C")],
                       c_minor_freq=0.1, seed=21)
P_B = _make_population("PopB", [("A", "C")] + [("T", "C")] * 9,
                       c_minor_freq=0.1, c_minor_freq_b=0.9, seed=22)

genotypes_split, mapping_split = _assemble({"PopA": P_A, "PopB": P_B})
pairs_within = la.compute_pairwise_ld(genotypes_split, LOCI, mapping_split,
                                      n_permutations=999, random_seed=1)
genotypes_pooled = {s: dict(c) for s, c in genotypes_split.items()}
mapping_pooled = {s: "Pooled" for s in genotypes_pooled}
pairs_pooled = la.compute_pairwise_ld(genotypes_pooled, LOCI, mapping_pooled,
                                      n_permutations=999, random_seed=1)


def _max_r2(pairs, locus1, locus2):
    values = [p.r2 for p in pairs
              if {p.locus1, p.locus2} == {locus1, locus2}
              and p.r2 is not None and np.isfinite(p.r2)]
    return max(values) if values else float("nan")


def test_within_population_ld_is_not_inflated():
    # independent loci inside each population must stay near zero
    assert _max_r2(pairs_within, "SNP_A", "SNP_C") < 0.3


def test_pooled_ld_is_inflated_by_structure():
    # the same genotypes pooled into one "population" show spurious LD
    assert _max_r2(pairs_pooled, "SNP_A", "SNP_C") > 0.3


# --------------------------------------------------------------------------
# 3. degenerate inputs
# --------------------------------------------------------------------------

def test_monomorphic_locus_is_skipped():
    genotypes_mono = {s: {**c, "SNP_D": ("A", "A")} for s, c in genotypes.items()}
    pairs = la.compute_pairwise_ld(genotypes_mono, ["SNP_A", "SNP_D"], mapping,
                                   n_permutations=99, random_seed=1)
    pair = next(p for p in pairs if {p.locus1, p.locus2} == {"SNP_A", "SNP_D"})
    assert pair.note or pair.r2 is None or not np.isfinite(pair.r2)


def test_missing_genotypes_reduce_n():
    broken = {s: dict(c) for s, c in genotypes.items()}
    for sample in list(broken)[:20]:
        broken[sample]["SNP_B"] = None
    pairs = la.compute_pairwise_ld(broken, ["SNP_A", "SNP_B"], mapping,
                                   n_permutations=99, random_seed=1)
    pair = pairs[0]
    assert pair.n_individuals == N - 20


def test_blocks_with_threshold_above_max_r2_are_singletons():
    summaries = la.summarise_ld_pairs(_PAIRS)
    observed_max = max(s.r2_max for s in summaries
                       if s.r2_max is not None and np.isfinite(s.r2_max))
    blocks = la.cluster_linkage_blocks(LOCI, summaries,
                                       r2_threshold=observed_max + 0.01)
    assert all(len(b.loci) == 1 for b in blocks)


# --------------------------------------------------------------------------
# 4. Excel report
# --------------------------------------------------------------------------

def test_report_is_written():
    summaries = la.summarise_ld_pairs(_PAIRS)
    blocks = la.cluster_linkage_blocks(LOCI, summaries, 0.8)
    la.choose_representatives(genotypes, sorted(genotypes), blocks)
    kept = [b.representative for b in blocks]
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "linkage_report.xlsx")
        la.write_linkage_report(path, _PAIRS, summaries, blocks, kept)
        assert os.path.getsize(path) > 0
        from openpyxl import load_workbook
        wb = load_workbook(path)
        assert wb.sheetnames == ["Pairwise_LD", "LD_Summary",
                                 "Linkage_Blocks", "Method"]
        assert wb["Pairwise_LD"]["A1"].value == "Locus1"
