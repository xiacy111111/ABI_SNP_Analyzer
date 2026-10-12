"""Tests for IUPAC base conversion logic."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _app_loader import load_app



class TestIupacConversion(unittest.TestCase):
    """Verify IUPAC ambiguity code expansion and conversion."""

    def _get_iupac_map(self):
        """Load the IUPAC_TO_BASES dict from the main module."""
        app = load_app()
        if not hasattr(app, "_IUPAC_TO_BASES"):
            self.skipTest("_IUPAC_TO_BASES not found in main module")
        return app._IUPAC_TO_BASES

    def test_homozygous_iupac(self):
        """Homozygous calls should map to a single-base list."""
        iupac = self._get_iupac_map()
        self.assertIn("A", iupac)
        self.assertIn("C", iupac)
        self.assertIn("G", iupac)
        self.assertIn("T", iupac)
        # Pure A expands to the single base A
        self.assertEqual(iupac["A"], ["A"])

    def test_heterozygous_iupac(self):
        """Heterozygous IUPAC codes should expand to two bases."""
        iupac = self._get_iupac_map()
        self.assertEqual(iupac.get("R"), ["A", "G"])  # A or G
        self.assertEqual(iupac.get("Y"), ["C", "T"])  # C or T
        self.assertEqual(iupac.get("S"), ["C", "G"])  # G or C
        self.assertEqual(iupac.get("W"), ["A", "T"])  # A or T
        self.assertEqual(iupac.get("K"), ["G", "T"])  # G or T
        self.assertEqual(iupac.get("M"), ["A", "C"])  # A or C

    def test_diploid_conversion(self):
        """IUPAC -> diploid genotype conversion."""
        app = load_app()
        if hasattr(app, "convert_iupac_to_diploid"):
            result = app.convert_iupac_to_diploid("R")
            # the implementation returns the two alleles as a joined string
            self.assertIn(str(result), ("AG", "GA", ("A", "G"), ["A", "G"]))
            self.assertEqual(len(result), 2)

    def test_ambiguous_iupac(self):
        """Ambiguous codes (N, ?) should either be missing or handled gracefully."""
        app = load_app()
        if hasattr(app, "_is_reliable_diploid_iupac"):
            self.assertFalse(app._is_reliable_diploid_iupac("N"))
            self.assertFalse(app._is_reliable_diploid_iupac("-"))
            self.assertFalse(app._is_reliable_diploid_iupac(""))


if __name__ == "__main__":
    unittest.main()
