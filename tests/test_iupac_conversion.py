"""Tests for IUPAC base conversion logic."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestIupacConversion(unittest.TestCase):
    """Verify IUPAC ambiguity code expansion and conversion."""

    def _get_iupac_map(self):
        """Load the IUPAC_TO_BASES dict from the main module."""
        import abi_to_snp_analyzer_v1_1 as app
        if not hasattr(app, "_IUPAC_TO_BASES"):
            self.skipTest("_IUPAC_TO_BASES not found in main module")
        return app._IUPAC_TO_BASES

    def test_homozygous_iupac(self):
        """Homozygous calls should produce single IUPAC code A/C/G/T."""
        iupac = self._get_iupac_map()
        self.assertIn("A", iupac)
        self.assertIn("C", iupac)
        self.assertIn("G", iupac)
        self.assertIn("T", iupac)
        # Pure A should expand to {A}
        self.assertEqual(iupac["A"], {"A"})

    def test_heterozygous_iupac(self):
        """Heterozygous IUPAC codes should expand to two bases."""
        iupac = self._get_iupac_map()
        self.assertEqual(iupac.get("R"), {"A", "G"})  # A or G
        self.assertEqual(iupac.get("Y"), {"C", "T"})  # C or T
        self.assertEqual(iupac.get("S"), {"G", "C"})  # G or C
        self.assertEqual(iupac.get("W"), {"A", "T"})  # A or T
        self.assertEqual(iupac.get("K"), {"G", "T"})  # G or T
        self.assertEqual(iupac.get("M"), {"A", "C"})  # A or C

    def test_diploid_conversion(self):
        """IUPAC -> diploid genotype conversion."""
        import abi_to_snp_analyzer_v1_1 as app
        if hasattr(app, "convert_iupac_to_diploid"):
            result = app.convert_iupac_to_diploid("R")
            self.assertIsInstance(result, tuple)
            self.assertEqual(len(result), 2)
            self.assertIn("A", result)
            self.assertIn("G", result)

    def test_ambiguous_iupac(self):
        """Ambiguous codes (N, ?) should either be missing or handled gracefully."""
        import abi_to_snp_analyzer_v1_1 as app
        if hasattr(app, "_is_reliable_diploid_iupac"):
            self.assertFalse(app._is_reliable_diploid_iupac("N"))
            self.assertFalse(app._is_reliable_diploid_iupac("-"))
            self.assertFalse(app._is_reliable_diploid_iupac(""))


if __name__ == "__main__":
    unittest.main()
