# -*- coding: utf-8 -*-
"""Round-trip tests for the SNP genotype import/export module.

Every writer must produce a file that the matching reader turns back into the
identical genotype table (calls, locus names, sample order and population
grouping), so that a file exported by ABI-SNP Analyzer can be re-imported with
"Load Existing Data" and fed straight into the population-genetics module.
"""

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import snp_format_io as snpio  # noqa: E402


# --------------------------------------------------------------------------
# Fixture data: 12 samples x 5 SNP loci in 3 populations, with missing calls
# --------------------------------------------------------------------------

LOCI = ["SNP01_A>C", "SNP02", "SNP03_G>T", "SNP04", "SNP05"]

SAMPLES = [
    "TL01", "TL02", "TL03", "TL04",
    "GH01", "GH02", "GH03", "GH04",
    "HZ01", "HZ02", "HZ03", "HZ04",
]

CALLS = {
    "SNP01_A>C": {"TL01": "AA", "TL02": "AC", "TL03": "CC", "TL04": "AA",
                  "GH01": "AA", "GH02": "AC", "GH03": "CC", "GH04": "--",
                  "HZ01": "CC", "HZ02": "CC", "HZ03": "AC", "HZ04": "AA"},
    "SNP02":     {"TL01": "TT", "TL02": "AT", "TL03": "AA", "TL04": "TT",
                  "GH01": "AT", "GH02": "TT", "GH03": "AA", "GH04": "AT",
                  "HZ01": "AA", "HZ02": "AT", "HZ03": "TT", "HZ04": "AA"},
    "SNP03_G>T": {"TL01": "GG", "TL02": "GT", "TL03": "TT", "TL04": "GG",
                  "GH01": "GG", "GH02": "GG", "GH03": "GT", "GH04": "TT",
                  "HZ01": "TT", "HZ02": "GT", "HZ03": "GG", "HZ04": "GG"},
    "SNP04":     {"TL01": "CC", "TL02": "CC", "TL03": "CG", "TL04": "GG",
                  "GH01": "CG", "GH02": "CC", "GH03": "GG", "GH04": "CC",
                  "HZ01": "GG", "HZ02": "GG", "HZ03": "CG", "HZ04": "CC"},
    "SNP05":     {"TL01": "AT", "TL02": "AA", "TL03": "TT", "TL04": "--",
                  "GH01": "TT", "GH02": "AT", "GH03": "AA", "GH04": "TT",
                  "HZ01": "AT", "HZ02": "TT", "HZ03": "--", "HZ04": "AT"},
}

POPULATIONS = {
    "Tilapia_A": ["TL01", "TL02", "TL03", "TL04"],
    "Tilapia_B": ["GH01", "GH02", "GH03", "GH04"],
    "Tilapia_C": ["HZ01", "HZ02", "HZ03", "HZ04"],
}

# Formats that keep individual IDs -> sample names must survive the round trip.
ID_PRESERVING = {"genepop", "arlequin", "genalex"}

# STRUCTURE has no header at all, so locus names cannot be carried.
LOCUS_NAME_LOSSY = {"structure"}

# FSTAT / hierfstat / POPGENE / STRUCTURE encode populations as bare numbers
# (FSTAT/hierfstat also carry no individual labels), so names cannot survive
# there by design.
POP_NAME_LOSSY = {"fstat", "hierfstat", "popgene", "structure"}


def _expected_order(fmt):
    _groups, ordered = snpio._resolve_populations(SAMPLES, POPULATIONS)
    return ordered


def _locus_pairs(fmt, data):
    """Pair fixture locus names with the names recovered from the file."""
    if fmt in LOCUS_NAME_LOSSY:
        return list(zip(LOCI, data.loci))
    return [(locus, locus) for locus in LOCI]


def _compare(fmt, data):
    """Compare a re-imported table against the fixture."""
    assert data.n_loci == len(LOCI), f"{fmt}: locus count changed"
    assert data.n_samples == len(SAMPLES), f"{fmt}: sample count changed"

    expected_order = _expected_order(fmt)
    assert len(expected_order) == len(data.samples)

    if fmt in ID_PRESERVING:
        assert set(data.samples) == set(SAMPLES), f"{fmt}: sample IDs changed"
    if fmt not in LOCUS_NAME_LOSSY:
        assert data.loci == LOCI, f"{fmt}: locus names changed: {data.loci}"

    for row, (expected_sample, actual_sample) in enumerate(
            zip(expected_order, data.samples)):
        for locus, recovered_locus in _locus_pairs(fmt, data):
            expected = CALLS[locus].get(expected_sample, "--")
            actual = data.calls[recovered_locus].get(actual_sample, "--")
            assert expected == actual, (
                f"{fmt}: row {row} locus {locus}: "
                f"expected {expected!r}, got {actual!r}")

    if fmt not in POP_NAME_LOSSY:
        assert sorted(data.populations) == sorted(POPULATIONS), \
            f"{fmt}: population names changed: {sorted(data.populations)}"
    assert len(data.populations) == len(POPULATIONS), f"{fmt}: population count changed"
    assert sorted(len(v) for v in data.populations.values()) == [4, 4, 4], \
        f"{fmt}: population sizes changed"


@pytest.mark.parametrize("fmt,writer,filename", [
    ("genepop", snpio.write_genepop, "roundtrip.txt"),
    ("arlequin", snpio.write_arlequin, "roundtrip.arp"),
    ("hierfstat", snpio.write_hierfstat, "roundtrip.hierfstat.txt"),
    ("genalex", snpio.write_genalex, "roundtrip.genalex.csv"),
    ("fstat", snpio.write_fstat, "roundtrip.dat"),
    ("popgene", snpio.write_popgene, "roundtrip.pop"),
    ("structure", snpio.write_structure, "roundtrip.str.txt"),
])
def test_export_import_roundtrip(fmt, writer, filename):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, filename)
        written = writer(path, SAMPLES, LOCI, CALLS, POPULATIONS)
        assert written == len(SAMPLES)
        assert os.path.getsize(path) > 0

        detected = snpio.detect_format(path)
        assert detected == fmt, f"{fmt}: auto-detected as {detected}"

        data = snpio.read_genotype_file(path)
        assert data.source_format == fmt
        _compare(fmt, data)


def test_hierfstat_layout():
    """hierfstat output must be a header + one row per individual table."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "h.txt")
        snpio.write_hierfstat(path, SAMPLES, LOCI, CALLS, POPULATIONS)
        with open(path, encoding="utf-8") as handle:
            lines = [ln.rstrip("\n") for ln in handle if ln.strip()]
        assert lines[0].split("\t")[0] == "pop"
        assert lines[0].split("\t")[1:] == LOCI
        assert len(lines) == len(SAMPLES) + 1
        # A/T heterozygote -> codes 1 and 2 -> "12"
        assert lines[1].split("\t")[1] == "11"      # TL01 SNP01 = AA
        assert lines[2].split("\t")[1] == "13"      # TL02 SNP01 = AC -> A=1, C=3
        assert lines[2].split("\t")[2] == "12"      # TL02 SNP02 = AT -> A=1, T=2
        assert len(lines[1].split("\t")) == len(LOCI) + 1


def test_genalex_layout():
    """GenAlEx output must carry the 3-line header block."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "g.csv")
        snpio.write_genalex(path, SAMPLES, LOCI, CALLS, POPULATIONS)
        with open(path, encoding="utf-8") as handle:
            lines = [ln.rstrip("\n") for ln in handle]
        head = lines[0].split(",")
        assert head[:3] == [str(len(LOCI)), str(len(SAMPLES)), "3"]
        assert head[3:] == ["4", "4", "4"]
        assert lines[2].split(",")[:2] == ["Sample", "Pop"]
        assert lines[2].split(",")[2] == LOCI[0]
        assert lines[2].split(",")[3] == ""
        assert len(lines) == 3 + len(SAMPLES)
        assert len(lines[3].split(",")) == 2 + 2 * len(LOCI)


def test_missing_data_is_preserved():
    """'-' must survive export/import in every format that supports it."""
    for fmt, writer, name in [
        ("genepop", snpio.write_genepop, "m.txt"),
        ("arlequin", snpio.write_arlequin, "m.arp"),
        ("hierfstat", snpio.write_hierfstat, "m.hf.txt"),
        ("genalex", snpio.write_genalex, "m.csv"),
    ]:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, name)
            writer(path, SAMPLES, LOCI, CALLS, POPULATIONS)
            data = snpio.read_genotype_file(path, fmt)
            order = _expected_order(fmt)
            for expected_sample, actual_sample in zip(order, data.samples):
                for locus, recovered in _locus_pairs(fmt, data):
                    if CALLS[locus][expected_sample] == "--":
                        assert data.calls[recovered][actual_sample] == "--", \
                            f"{fmt}: missing call lost at {locus}/{expected_sample}"


def test_single_population_falls_back_to_pop1():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "single.arp")
        snpio.write_arlequin(path, SAMPLES, LOCI, CALLS, None)
        data = snpio.read_genotype_file(path, "arlequin")
        assert list(data.populations) == ["Pop1"]
        assert data.n_samples == len(SAMPLES)


def test_forced_format_overrides_detection():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "forced.txt")
        snpio.write_hierfstat(path, SAMPLES, LOCI, CALLS, POPULATIONS)
        # Force the wrong parser on purpose -> must raise a clear error.
        with pytest.raises(snpio.SnpFormatError):
            snpio.read_genotype_file(path, "genepop")


def test_allele_coding_convention():
    assert snpio.DEFAULT_CODE_TO_BASE == {1: "A", 2: "T", 3: "C", 4: "G"}
    assert snpio._allele_code("A") == 1
    assert snpio._allele_code("G") == 4
    assert snpio._allele_code(None) == 0


def test_imported_data_feeds_population_analysis():
    """Re-imported hierfstat data must be analysable (2+ pops, parsed calls)."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "analysis.txt")
        snpio.write_hierfstat(path, SAMPLES, LOCI, CALLS, POPULATIONS)
        data = snpio.read_genotype_file(path, "hierfstat")
        matrix = data.to_matrix_dict()

        parsed = 0
        for locus in data.loci:
            by_pop = {}
            for pop, members in data.populations.items():
                for sample in members:
                    call = matrix[locus][sample]
                    if call != "--":
                        by_pop.setdefault(pop, []).append(call)
                        parsed += 1
            assert len(by_pop) == 3
        assert parsed > 0
