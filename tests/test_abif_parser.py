"""Tests for ABIF parser: reading chromatogram binary format."""

import os
import sys
import unittest

# Make the parent package importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _app_loader import load_app



class TestAbifParser(unittest.TestCase):
    """Verify ABIF binary parsing correctness."""

    @unittest.skipIf(
        not os.path.isdir(os.path.join(os.path.dirname(__file__), "..", "example")),
        "Example directory not found",
    )
    def test_parse_example_ab1(self):
        """Parse the first available .ab1 file and check basic structure."""
        example_dir = os.path.join(os.path.dirname(__file__), "..", "example")
        ab1_files = [
            f for f in os.listdir(example_dir)
            if f.lower().endswith(".ab1")
        ]
        if not ab1_files:
            self.skipTest("No .ab1 files in example/")
        # Use the smallest available file
        ab1_path = os.path.join(example_dir, sorted(ab1_files)[0])

        AbifParser = load_app().AbifParser
        parser = AbifParser()
        trace = parser.read_abi(ab1_path)

        self.assertIsNotNone(trace, "read_abi should return a Trace")
        self.assertIsNotNone(trace.base_calls, "Base calls should not be None")
        self.assertGreater(len(trace.base_calls), 0,
                           "Should have at least one base call")

    def test_invalid_files_are_rejected(self):
        """Empty and truncated files must be rejected gracefully (no Trace)."""
        import tempfile
        AbifParser = load_app().AbifParser
        parser = AbifParser()

        fd, tmp = tempfile.mkstemp(suffix=".ab1")
        os.close(fd)
        try:
            self.assertIsNone(parser.read_abi(tmp),
                              "an empty file must not yield a Trace")
        finally:
            os.remove(tmp)

        fd, tmp = tempfile.mkstemp(suffix=".ab1")
        os.write(fd, b"ABIF" + bytes(range(1, 5)))
        os.close(fd)
        try:
            self.assertIsNone(parser.read_abi(tmp),
                              "a truncated file must not yield a Trace")
        finally:
            os.remove(tmp)


if __name__ == "__main__":
    unittest.main()
