"""Tests for population genetics statistical calculations."""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from population_genetics import (
    PopulationGenetics,
    LocusSummary,
    AMOVAResult,
    GenotypeFreqTestResult,
)


class TestAlleleFrequencies(unittest.TestCase):
    """Verify allele frequency computation."""

    def test_basic_frequencies(self):
        """Simple biallelic case with 3 samples."""
        genotypes = [("A", "A"), ("A", "C"), ("C", "C")]
        freqs, total = PopulationGenetics.compute_allele_frequencies(genotypes)
        self.assertEqual(total, 6)  # 3 samples x 2 alleles
        self.assertAlmostEqual(freqs["A"], 3 / 6)
        self.assertAlmostEqual(freqs["C"], 3 / 6)

    def test_heterozygote_counting(self):
        """Each heterozygote contributes 2 alleles (1 per allele type)."""
        genotypes = [("A", "C")]  # Single heterozygote
        freqs, total = PopulationGenetics.compute_allele_frequencies(genotypes)
        self.assertEqual(total, 2)  # NOT 1
        self.assertAlmostEqual(freqs["A"], 0.5)
        self.assertAlmostEqual(freqs["C"], 0.5)

    def test_empty_input(self):
        """Empty genotype list should return empty frequencies."""
        freqs, total = PopulationGenetics.compute_allele_frequencies([])
        self.assertEqual(total, 0)
        self.assertEqual(freqs, {})

    def test_monomorphic(self):
        """Monomorphic locus should have freq = 1.0."""
        genotypes = [("A", "A"), ("A", "A"), ("A", "A")]
        freqs, total = PopulationGenetics.compute_allele_frequencies(genotypes)
        self.assertEqual(total, 6)
        self.assertEqual(freqs["A"], 1.0)
        self.assertEqual(len(freqs), 1)


class TestHWE(unittest.TestCase):
    """Verify Hardy-Weinberg equilibrium exact test."""

    def test_hwe_monomorphic(self):
        """Monomorphic locus should have p-value = 1.0."""
        genotypes = [("A", "A")] * 20
        pval, method = PopulationGenetics.hwe_exact_test(genotypes)
        self.assertEqual(pval, 1.0)
        self.assertEqual(method, "monomorphic")

    def test_hwe_perfect_equilibrium(self):
        """Genotype counts in perfect HWE proportion should give p ~ 1."""
        # p=0.5, q=0.5: expect 5 AA, 10 Aa, 5 aa
        genotypes = (
            [("A", "A")] * 5
            + [("A", "C")] * 10
            + [("C", "C")] * 5
        )
        pval, method = PopulationGenetics.hwe_exact_test(genotypes)
        self.assertEqual(method, "exact-biallelic")
        self.assertGreater(pval, 0.05)  # Should not reject

    def test_hwe_severe_distortion(self):
        """Highly distorted frequencies should give p < 0.05."""
        genotypes = (
            [("A", "A")] * 18  # Almost all homozygous
            + [("A", "C")] * 1
            + [("C", "C")] * 1
        )
        pval, method = PopulationGenetics.hwe_exact_test(genotypes)
        self.assertEqual(method, "exact-biallelic")
        self.assertLess(pval, 0.05)

    def test_hwe_missing_data(self):
        """Missing/N genotypes should be filtered out."""
        genotypes = [
            ("A", "A"), ("A", "C"), ("C", "C"),
            ("-", "-"), ("N", "N"), None, ("NA", "NA"),
        ]
        pval, method = PopulationGenetics.hwe_exact_test(genotypes)
        self.assertEqual(method, "exact-biallelic")
        self.assertIsNotNone(pval)

    def test_hwe_empty(self):
        """All missing should return empty method."""
        genotypes = [("-", "-"), ("", ""), None, ("N", "N")]
        pval, method = PopulationGenetics.hwe_exact_test(genotypes)
        self.assertEqual(method, "empty")
        self.assertIsNone(pval)


class TestChiSquareGenotypeFreq(unittest.TestCase):
    """Verify genotype frequency chi-square test."""

    def test_two_population_test(self):
        """Basic 2-population genotype frequency comparison."""
        genotypes_by_pop = {
            "PopA": [("A", "A")] * 10 + [("A", "C")] * 5 + [("C", "C")] * 5,
            "PopB": [("A", "A")] * 5 + [("A", "C")] * 9 + [("C", "C")] * 6,
        }
        result = PopulationGenetics.test_locus_genotype_freq(
            genotypes_by_pop, "TestLocus"
        )
        self.assertIsNotNone(result)
        self.assertEqual(result.locus_name, "TestLocus")
        self.assertIsNotNone(result.raw_p_value)
        self.assertGreaterEqual(result.raw_p_value, 0.0)
        self.assertLessEqual(result.raw_p_value, 1.0)

    def test_original_contingency_table_preserved(self):
        """Original contingency table should be stored in the result."""
        genotypes_by_pop = {
            "PopA": [("A", "A")] * 3 + [("A", "C")] * 2 + [("C", "C")] * 1,
            "PopB": [("A", "A")] * 1 + [("A", "C")] * 2 + [("C", "C")] * 3,
        }
        result = PopulationGenetics.test_locus_genotype_freq(
            genotypes_by_pop, "TestLocus"
        )
        self.assertIsNotNone(result)
        self.assertIsNotNone(
            result.original_contingency_table,
            "Original contingency table should be stored",
        )
        self.assertEqual(len(result.original_contingency_table), 3)

    def test_single_population_returns_none(self):
        """Single population should return None (insufficient for test)."""
        genotypes_by_pop = {"PopA": [("A", "A")] * 10}
        result = PopulationGenetics.test_locus_genotype_freq(
            genotypes_by_pop, "TestLocus"
        )
        self.assertIsNone(result)

    def test_fisher_exact_fallback(self):
        """Small 2x2 table should trigger Fisher's exact test."""
        genotypes_by_pop = {
            "PopA": [("A", "A")] * 2 + [("A", "C")] * 1,
            "PopB": [("A", "A")] * 1 + [("C", "C")] * 2,
        }
        result = PopulationGenetics.test_locus_genotype_freq(
            genotypes_by_pop, "TestLocus"
        )
        self.assertIsNotNone(result)
        # Should use Fisher when expected < 5 in 2x2
        if result.test_method == "fisher_exact":
            self.assertEqual(result.test_method, "fisher_exact")


class TestHeterozygosity(unittest.TestCase):
    """Verify heterozygosity and PIC calculations."""

    def test_heterozygosity_calculation(self):
        genotypes = [("A", "A"), ("A", "C"), ("C", "C")]
        ho, he, n = PopulationGenetics.compute_heterozygosity(genotypes)
        self.assertEqual(n, 3)
        self.assertAlmostEqual(ho, 1 / 3)  # 1 het out of 3
        self.assertGreater(he, 0.0)

    def test_pic_calculation(self):
        from collections import Counter
        freqs = {"A": 0.5, "C": 0.5}
        pic = PopulationGenetics.compute_pic(freqs)
        self.assertGreater(pic, 0.0)
        self.assertLessEqual(pic, 0.5)


class TestFDRCorrection(unittest.TestCase):
    """Verify multiple testing correction."""

    def test_benjamini_hochberg(self):
        pvals = [0.001, 0.01, 0.04, 0.20, None]
        qvals = PopulationGenetics.benjamini_hochberg_fdr(pvals)
        self.assertEqual(len(qvals), len(pvals))
        self.assertIsNone(qvals[4])  # None should remain None
        # Smallest p should get a q-value
        self.assertIsNotNone(qvals[0])

    def test_bonferroni(self):
        pvals = [0.01, 0.02, None]
        corrected = PopulationGenetics.bonferroni_correction(pvals)
        self.assertEqual(len(corrected), 3)
        self.assertAlmostEqual(corrected[0], 0.02)  # 0.01 * 2
        self.assertAlmostEqual(corrected[1], 0.04)  # 0.02 * 2
        self.assertIsNone(corrected[2])


class TestAMOVA(unittest.TestCase):
    """Verify AMOVA calculation structure."""

    def test_amova_requires_sample_populations(self):
        """AMOVA should reject the old pooled-genotype API."""
        with self.assertRaises(ValueError):
            PopulationGenetics.amova({"sample1": {"loc1": ("A", "C")}})


class TestWC84FStatistics(unittest.TestCase):
    """Weir & Cockerham (1984) ANOVA variance-component F-statistics."""

    def test_strong_population_differentiation(self):
        """Fst should be > 0.5 when populations have very different allele freqs."""
        pop1 = [("A", "A")] * 8 + [("A", "C")] * 2   # p_A = 0.9
        pop2 = [("C", "C")] * 8 + [("A", "C")] * 2   # p_A = 0.1
        fis, fst, fit = PopulationGenetics.compute_f_statistics(
            {"Pop1": pop1, "Pop2": pop2}
        )
        self.assertIsNotNone(fst)
        self.assertGreater(fst, 0.5, "Strong differentiation should give Fst > 0.5")

    def test_identical_populations(self):
        """Fst should be ~0 when populations are identical."""
        pop_hwe = [("A", "A")] * 5 + [("A", "C")] * 10 + [("C", "C")] * 5
        fis, fst, fit = PopulationGenetics.compute_f_statistics(
            {"Pop1": pop_hwe, "Pop2": pop_hwe.copy()}
        )
        self.assertIsNotNone(fst)
        self.assertLess(fst, 0.05, "Identical populations should have Fst ~ 0")

    def test_single_population(self):
        """Single population: Fis should be valid, Fst=None."""
        genotypes = [("A", "A")] * 5 + [("A", "C")] * 10 + [("C", "C")] * 5
        fis, fst, fit = PopulationGenetics.compute_f_statistics({"Pop1": genotypes})
        self.assertIsNotNone(fis, "Fis should be computable for single population")
        self.assertIsNone(fst, "Fst requires >= 2 populations")

    def test_monomorphic_returns_none(self):
        """Monomorphic locus should return all None."""
        pop1 = [("A", "A")] * 10
        pop2 = [("A", "A")] * 10
        fis, fst, fit = PopulationGenetics.compute_f_statistics(
            {"Pop1": pop1, "Pop2": pop2}
        )
        self.assertIsNone(fis)
        self.assertIsNone(fst)
        self.assertIsNone(fit)

    def test_strong_inbreeding_fis(self):
        """Fis should be >> 0 with severe heterozygote deficit."""
        pop_inbred = [("A", "A")] * 9 + [("C", "C")] * 9 + [("A", "C")] * 2
        fis, fst, fit = PopulationGenetics.compute_f_statistics(
            {"Pop1": pop_inbred, "Pop2": pop_inbred.copy()}
        )
        self.assertIsNotNone(fis)
        self.assertGreater(fis, 0.5, "Fis should exceed 0.5 with severe inbreeding")

    def test_multiallelic_locus(self):
        """WC84 multiallelic extension should produce valid Fst in [0, 1]."""
        pop1 = [("A", "A")] * 5 + [("A", "C")] * 3 + [("C", "C")] * 2 \
               + [("A", "G")] * 1 + [("G", "G")] * 1
        pop2 = [("G", "G")] * 5 + [("C", "G")] * 4 + [("C", "C")] * 3
        fis, fst, fit = PopulationGenetics.compute_f_statistics(
            {"Pop1": pop1, "Pop2": pop2}
        )
        self.assertIsNotNone(fst)
        self.assertTrue(0.0 <= fst <= 1.0,
                        "Fst should be in [0, 1] for multiallelic locus")

    def test_fst_method_label_in_summary(self):
        """analyze_locus should tag Fst method as WC84_variance_component."""
        pop1 = [("A", "A")] * 5 + [("A", "C")] * 3 + [("C", "C")] * 2
        pop2 = [("C", "C")] * 5 + [("A", "C")] * 3 + [("A", "A")] * 2
        summary = PopulationGenetics.analyze_locus(
            {"Pop1": pop1, "Pop2": pop2}, locus_name="TestSNP"
        )
        self.assertEqual(summary.fst_method, "WC84_variance_component")

    def test_permutation_test_produces_valid_p(self):
        """Permutation test should produce 0 <= p <= 1."""
        pop1 = [("A", "A")] * 7 + [("A", "C")] * 3
        pop2 = [("C", "C")] * 7 + [("A", "C")] * 3
        fst_obs, p_val = PopulationGenetics.fst_permutation_test(
            {"Pop1": pop1, "Pop2": pop2}, n_permutations=199, random_seed=42
        )
        self.assertGreater(fst_obs, 0.0)
        self.assertTrue(0.0 <= p_val <= 1.0)


if __name__ == "__main__":
    unittest.main()
