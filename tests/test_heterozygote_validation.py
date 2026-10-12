# -*- coding: utf-8 -*-
"""Tests for heterozygote_validation.py (Major 3(b) tooling)."""

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import heterozygote_validation as hv  # noqa: E402


def test_parse_genotype_accepts_common_notations():
    assert hv.parse_genotype("A/T") == ("A", "T")
    assert hv.parse_genotype("T/A") == ("T", "A")
    assert hv.parse_genotype("AT") == ("A", "T")
    assert hv.parse_genotype("W") == ("A", "T")
    assert hv.parse_genotype("y") == ("C", "T")
    assert hv.parse_genotype("G") == ("G", "G")
    for token in ("", "--", "N", "NA", ".", None):
        assert hv.parse_genotype(token) == (None, None)
    assert hv.parse_genotype("ZX") == (None, None)      # not a DNA call


def test_zygosity_and_labels():
    assert hv.zygosity(("A", "T")) == hv.HET
    assert hv.zygosity(("A", "A")) == hv.HOM
    assert hv.zygosity((None, None)) == hv.MISS
    assert hv.genotype_label(("T", "A")) == "A/T"


def test_wilson_interval_is_conservative_for_small_counts():
    low, high = hv.wilson_interval(1, 1)
    assert low is not None and 0.20 < low < 1.0      # not a fake 100%
    low, high = hv.wilson_interval(0, 880)
    assert low == 0.0 and high < 0.01


def _truth_and_software():
    """8 calls: 3 truth-HET (1 missed), 4 truth-HOM (1 false het), 1 no-call."""
    truth = {
        "S1": {"L1": ("A", "T"), "L2": ("A", "A"), "L3": ("C", "C")},
        "S2": {"L1": ("A", "A"), "L2": ("C", "G"), "L3": ("G", "G")},
        "S3": {"L1": ("T", "T"), "L2": ("A", "T"), "L3": (None, None)},
    }
    software = {
        "S1": {"L1": ("A", "A"), "L2": ("A", "A"), "L3": ("C", "C")},   # missed het
        "S2": {"L1": ("A", "T"), "L2": ("C", "G"), "L3": ("G", "G")},   # false het
        "S3": {"L1": ("T", "T"), "L2": ("A", "A"), "L3": (None, None)},  # missed het
    }
    return truth, software


def test_confusion_matrix_and_class_metrics():
    result = hv.compare_genotypes(*_truth_and_software())

    assert result.confusion[hv.HET][hv.HET] == 1
    assert result.confusion[hv.HET][hv.HOM] == 2
    assert result.confusion[hv.HOM][hv.HET] == 1
    assert result.confusion[hv.HOM][hv.HOM] == 4
    assert result.confusion[hv.MISS][hv.MISS] == 1

    assert result.n_comparable == 8
    assert result.n_concordant == 5
    assert abs(result.concordance - 5 / 8) < 1e-9

    assert result.het.truth_total == 3
    assert result.het.predicted_total == 2
    assert result.het.tp == 1
    assert result.het.fn == 2
    assert result.het.fp == 1
    assert abs(result.het.sensitivity - 1 / 3) < 1e-9
    assert abs(result.het.precision - 1 / 2) < 1e-9

    assert result.hom.truth_total == 5
    assert result.hom.predicted_total == 6
    assert result.hom.tp == 4
    assert result.hom.fp == 2
    assert result.hom.fn == 1


def test_discordant_calls_are_classified():
    result = hv.compare_genotypes(*_truth_and_software())
    categories = {(row["Sample"], row["Locus"]): row["Category"]
                  for row in result.discordant}
    assert categories[("S1", "L1")] == \
        "missed heterozygote (software HOM, truth HET)"
    assert categories[("S2", "L1")] == \
        "false heterozygote (software HET, truth HOM)"
    assert categories[("S3", "L2")] == \
        "missed heterozygote (software HOM, truth HET)"
    assert len(result.discordant) == 3


def test_perfect_agreement_gives_unit_metrics():
    truth, _ = _truth_and_software()
    result = hv.compare_genotypes(truth, truth)
    assert result.concordance == 1.0
    assert result.het.sensitivity == 1.0 and result.het.precision == 1.0
    assert result.discordant == []


def test_truth_table_readers_wide_and_long(tmp_path):
    import pandas as pd

    wide = pd.DataFrame({
        "Sample": ["S1", "S2"],
        "L1": ["A/T", "A/A"],
        "L2": ["W", "G"],
    })
    wide_path = tmp_path / "truth_wide.csv"
    wide.to_csv(wide_path, index=False)
    truth_wide = hv.read_truth_table(str(wide_path))
    assert truth_wide["S1"]["L1"] == ("A", "T")
    assert truth_wide["S1"]["L2"] == ("A", "T")     # IUPAC W
    assert truth_wide["S2"]["L2"] == ("G", "G")

    long = pd.DataFrame({
        "Sample": ["S1", "S1", "S2"],
        "Locus": ["L1", "L2", "L1"],
        "Genotype": ["A/T", "A/A", "T/T"],
    })
    long_path = tmp_path / "truth_long.csv"
    long.to_csv(long_path, index=False)
    truth_long = hv.read_truth_table(str(long_path))
    assert truth_long["S1"]["L2"] == ("A", "A")
    assert truth_long["S2"]["L1"] == ("T", "T")


def test_report_is_written(tmp_path):
    import openpyxl

    truth, software = _truth_and_software()
    result = hv.compare_genotypes(truth, software)
    out = tmp_path / "validation_report.xlsx"
    hv.write_report(str(out), result, "truth.csv", "genotype_matrix.xlsx")
    workbook = openpyxl.load_workbook(str(out))
    assert workbook.sheetnames == ["Summary", "Confusion_Matrix", "Metrics",
                                   "Discordant_Calls"]
    assert workbook["Discordant_Calls"].max_row == 1 + len(result.discordant)


def test_confusion_figure_is_written(tmp_path):
    pytest.importorskip("matplotlib")
    truth, software = _truth_and_software()
    result = hv.compare_genotypes(truth, software)
    plot = tmp_path / "confusion.png"
    hv.plot_confusion(result, str(plot))
    assert os.path.getsize(str(plot)) > 0


def test_cli_end_to_end(tmp_path):
    import pandas as pd

    truth, software = _truth_and_software()
    truth_path = tmp_path / "truth.xlsx"
    pd.DataFrame([
        {"Sample": s, **{locus: hv.genotype_label(calls[locus])
                         for locus in ("L1", "L2", "L3")}}
        for s, calls in truth.items()
    ]).to_excel(str(truth_path), index=False)

    software_path = tmp_path / "genotype_matrix.xlsx"
    pd.DataFrame([
        {"Sample": s, **{locus: hv.genotype_label(calls[locus])
                         for locus in ("L1", "L2", "L3")}}
        for s, calls in software.items()
    ]).to_excel(str(software_path), index=False)

    out = tmp_path / "report.xlsx"
    exit_code = hv.main(["--software", str(software_path),
                         "--truth", str(truth_path), "--out", str(out)])
    assert exit_code == 0
    assert out.exists()
