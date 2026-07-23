"""Tests for ABIF parser: reading chromatogram binary format."""

import os
import sys
import unittest

# Make the parent package importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


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

        from abi_to_snp_analyzer_v1_1 import AbifParser
        parser = AbifParser(ab1_path)
        parser.parse()

        # Should have base calls
        self.assertIsNotNone(parser.base_calls, "Base calls should not be None")
        self.assertGreater(len(parser.base_calls), 0, "Should have at least one base call")

    def test_truncated_file_handling(self):
        """Truncated/corrupt files should raise a clear error."""
        import tempfile
        fd, tmp = tempfile.mkstemp(suffix=".ab1")
        os.write(fd, b"ABIF\x00\x01\x02\x03")
        os.close(fd)

        from abi_to_snp_analyzer_v1_1 import AbifParser
        parser = AbifParser(tmp)
        with self.assertRaises((ValueError, IndexError, KeyError, struct.error)):
            parser.parse()
        os.unlink(tmp)

    def test_empty_file_handling(self):
        """Empty file should raise an error immediately."""
        import tempfile
        fd, tmp = tempfile.mkstemp(suffix=".ab1")
        os.close(fd)

        from abi_to_snp_analyzer_v1_1 import AbifParser
        parser = AbifParser(tmp)
        with self.assertRaises(Exception):
            parser.parse()
        os.unlink(tmp)


if __name__ == "__main__":
    unittest.main()
