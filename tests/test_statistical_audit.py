# -*- coding: utf-8 -*-
"""Statistical audit tests (SoftwareX Major 4v follow-up).

Covers the corrections applied after the v1.0.0 statistical review:
  * He must not receive the Nei unbiasedness correction twice
  * AMOVA permutations must move individuals (both haplotypes together)
  * genotype labels must not depend on the allele order in the source matrix
  * R x C genotype-frequency tables with small expected counts must use a
    permutation p-value instead of the asymptotic chi-square
  * WC84 components must stay in [0, 1] and Fst must separate structured
    populations
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from population_genetics import PopulationGenetics  # noqa: E402


def _structured_populations(n_per_pop=25, seed=1):
    """Two populations with clearly different allele frequencies."""
    import random
    rng = random.Random(seed)

    def sample(minor_freq):
        calls = []
        for _ in range(n_per_pop):
            u = rng.random()
            if u < (1 - minor_freq) ** 2:
                calls.append(("A", "A"))
            elif u < 1 - minor_freq ** 2:
                calls.append(("A", "T"))
            else:
                calls.append(("T", "T"))
        return calls

    return {"PopA": sample(0.05), "PopB": sample(0.60)}


class TestHeterozygosityCorrection(unittest.TestCase):
    """Nei's unbiased correction must be applied exactly once."""

    def test_he_matches_hand_computed_value(self):
        genotypes = [("A", "A")] * 6 + [("A", "T")] * 4
        ho, he, n = PopulationGenetics.compute_heterozygosity(genotypes)
        # p_A = 16/20 = 0.8, p_T = 0.2 -> MLE He = 1 - (0.64 + 0.04) = 0.32
        expected = (2 * n / (2 * n - 1)) * 0.32
        self.assertAlmostEqual(he, expected, places=12)
        self.assertAlmostEqual(ho, 0.4, places=12)

    def test_analyze_locus_does_not_double_correct(self):
        by_pop = _structured_populations()
        summary = PopulationGenetics.analyze_locus(by_pop, "L1")
        for stats in summary.per_population.values():
            self.assertAlmostEqual(stats.he, stats.he_unbiased, places=12,
                                   msg="He was unbiased-corrected twice")


class TestAmovaPermutation(unittest.TestCase):
    """AMOVA permutations must preserve the individual (2-haplotype) unit."""

    def test_haplotypes_are_paired_by_individual(self):
        sample_genotypes = {
            "S1": {"L1": ("A", "T"), "L2": ("C", "G")},
            "S2": {"L1": ("A", "A"), "L2": ("G", "G")},
            "S3": {"L1": ("T", "T"), "L2": ("C", "C")},
            "S4": {"L1": ("A", "T"), "L2": ("G", "C")},
        }
        populations = {"S1": "P1", "S2": "P1", "S3": "P2", "S4": "P2"}
        _ids, labels, _dist = PopulationGenetics._pairwise_difference_haplotypes(
            sample_genotypes, populations)
        self.assertEqual(len(labels), 8)
        for k in range(4):
            self.assertEqual(labels[2 * k], labels[2 * k + 1],
                             "an individual's two haplotypes must share a label")

    def test_amova_runs_and_p_value_is_valid(self):
        sample_genotypes = {}
        populations = {}
        index = 0
        for pop, minor_freq in (("PopA", 0.05), ("PopB", 0.60)):
            for _ in range(20):
                sample = f"S{index:03d}"
                index += 1
                import random
                u = random.Random(index).random()
                dosage = 2 if u < minor_freq ** 2 else (
                    1 if u < 1 - (1 - minor_freq) ** 2 else 0)
                sample_genotypes[sample] = {
                    "L1": [("A", "A"), ("A", "T"), ("T", "T")][dosage],
                    "L2": [("C", "C"), ("C", "G"), ("G", "G")][dosage],
                }
                populations[sample] = pop
        results = PopulationGenetics.amova(
            sample_genotypes, populations,
            n_permutations=199, random_seed=7)
        phi = results[0].phi_stat
        p_value = results[0].p_value
        self.assertIsNotNone(phi)
        self.assertGreaterEqual(phi, 0.0)
        self.assertLessEqual(phi, 1.0)
        self.assertGreater(p_value, 0.0)
        self.assertLessEqual(p_value, 1.0)
        # clearly differentiated populations must be detected by the
        # individual-level permutation scheme
        self.assertLessEqual(p_value, 0.05)


class TestGenotypeLabelCanonicalisation(unittest.TestCase):
    """Allele order in the source matrix must not split genotype classes."""

    def test_swapped_alleles_land_in_one_class(self):
        by_pop = {
            "PopA": [("T", "A")] * 3 + [("A", "A")] * 9,
            "PopB": [("A", "T")] * 3 + [("C", "C")] * 9,
        }
        result = PopulationGenetics.test_locus_genotype_freq(by_pop, "L1")
        self.assertIsNotNone(result)
        labels = set(result.genotype_labels)
        self.assertIn("A/T", labels)
        self.assertNotIn("T/A", labels)


class TestGenotypeFreqSmallSamples(unittest.TestCase):
    """Small expected counts must switch to a permutation p-value."""

    def test_permutation_used_for_sparse_rxc_table(self):
        # 3 genotype classes x 3 populations; the TT row merges into "Others",
        # leaving a 3 x 3 table whose "Others" column still has expected < 5,
        # so the asymptotic chi-square must be replaced by a permutation p.
        by_pop = {
            "PopA": [("A", "A")] * 10 + [("A", "T")] * 9 + [("T", "T")] * 1,
            "PopB": [("A", "A")] * 9 + [("A", "T")] * 8 + [("T", "T")] * 3,
            "PopC": [("A", "A")] * 10 + [("A", "T")] * 9 + [("T", "T")] * 1,
        }
        result = PopulationGenetics.test_locus_genotype_freq(
            by_pop, "L1", n_permutations=299, random_seed=3)
        self.assertIsNotNone(result)
        self.assertEqual(result.test_method, "chi2_permutation")
        self.assertIsNotNone(result.raw_p_value)
        self.assertGreater(result.raw_p_value, 0.0)
        self.assertLessEqual(result.raw_p_value, 1.0)

    def test_asymptotic_chi2_kept_for_large_counts(self):
        by_pop = {
            "PopA": [("A", "A")] * 30 + [("A", "T")] * 10,
            "PopB": [("A", "A")] * 10 + [("A", "T")] * 30,
        }
        result = PopulationGenetics.test_locus_genotype_freq(
            by_pop, "L1", n_permutations=99, random_seed=3)
        self.assertIsNotNone(result)
        self.assertEqual(result.test_method, "chi2")
        self.assertLess(result.raw_p_value, 0.05)


class TestWC84Sanity(unittest.TestCase):
    """WC84 components must be valid and separate structured populations."""

    def test_fst_bounds_and_structure(self):
        structured = _structured_populations()
        fis, fst, fit = PopulationGenetics.compute_f_statistics(structured)
        self.assertIsNotNone(fst)
        for value in (fis, fst, fit):
            self.assertIsNotNone(value)
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)
        self.assertGreater(fst, 0.05)
        # panmictic-like populations (same frequencies, random genotypes)
        import random
        rng = random.Random(11)

        def sample():
            calls = []
            for _ in range(30):
                u = rng.random()
                calls.append(("A", "A") if u < 0.64 else
                             ("A", "T") if u < 0.96 else ("T", "T"))
            return calls
        _, fst_same, _ = PopulationGenetics.compute_f_statistics(
            {"P1": sample(), "P2": sample()})
        self.assertIsNotNone(fst_same)
        self.assertLess(fst_same, fst)


if __name__ == "__main__":
    unittest.main()
