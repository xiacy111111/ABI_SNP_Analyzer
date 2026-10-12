# -*- coding: utf-8 -*-
"""Tests for the VCF 4.2 export added in v1.1.0.

VCF is a write-only format (the GUI does not re-import it), so these tests
validate the emitted text directly: header structure, coordinate resolution
(genome-mapping metadata > reference-forward SNP name > amplicon fallback),
genotype encoding, missing calls, multi-allelic ALTs and the minus-strand
complementation of alleles.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import snp_format_io as snpio  # noqa: E402


def _write(tmp_path, samples, loci, calls, **kwargs):
    path = os.path.join(str(tmp_path), "out.vcf")
    snpio.write_vcf(path, samples, loci, calls, **kwargs)
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read().splitlines()


def _records(lines):
    return [line.split("\t") for line in lines if line and not
            line.startswith("#")]


def _headers(lines):
    return [line for line in lines if line.startswith("##")]


def test_header_and_genome_named_locus(tmp_path):
    calls = {"NC_031982.2:26437556C>T": {
        "s1": "CC", "s2": "CT", "s3": "--", "s4": "TT"}}
    lines = _write(tmp_path, ["s1", "s2", "s3", "s4"],
                   ["NC_031982.2:26437556C>T"], calls,
                   source="ABI-SNP_Analyzer_v1.1.0")
    assert lines[0] == "##fileformat=VCFv4.2"
    assert "##source=ABI-SNP_Analyzer_v1.1.0" in lines
    assert any(line.startswith("##contig=<ID=NC_031982.2>") for line in lines)
    header = [line for line in lines if line.startswith("#CHROM")][0]
    assert header.split("\t")[:9] == ["#CHROM", "POS", "ID", "REF", "ALT",
                                     "QUAL", "FILTER", "INFO", "FORMAT"]
    assert header.split("\t")[9:] == ["s1", "s2", "s3", "s4"]
    record = _records(lines)[0]
    assert record[:5] == ["NC_031982.2", "26437556",
                          "NC_031982.2:26437556C>T", "C", "T"]
    assert record[7] == "NS=3;AF=0.500000"
    assert record[9:] == ["0/0", "0/1", "./.", "1/1"]


def test_snp_index_fallback_contig(tmp_path):
    calls = {"SNP_62": {"s1": "CC", "s2": "CT", "s3": "TT"}}
    lines = _write(tmp_path, ["s1", "s2", "s3"], ["SNP_62"], calls)
    assert any(line.startswith("##ABI-SNP_Analyzer_note=") for line in lines)
    record = _records(lines)[0]
    assert record[0] == snpio.VCF_FALLBACK_CONTIG
    assert record[1] == "62"
    assert record[3] == "C"
    assert record[4] == "T"
    assert record[9:] == ["0/0", "0/1", "1/1"]


def test_metadata_mapping_minus_strand_complement(tmp_path):
    """A minus-strand mapping must complement alleles to reference-forward."""
    mapping = {31: {
        "mapping_status": "unique",
        "reference_sequence": "Chr1",
        "genome_position": 100,
        "reference_base": "A",
        "allele_complement_required": True,
    }}
    # sequencing strand: T = reference A on the forward strand,
    #                    G = C on the forward strand
    calls = {"SNP_31": {"s1": "TT", "s2": "TG", "s3": "GG"}}
    metadata = snpio.vcf_locus_metadata_from_mapping(
        ["SNP_31"], calls, mapping)
    assert metadata["SNP_31"] == {
        "chrom": "Chr1", "pos": 100, "ref": "A", "alts": ["C"],
        "complement": True}
    lines = _write(tmp_path, ["s1", "s2", "s3"], ["SNP_31"], calls,
                   locus_metadata=metadata)
    record = _records(lines)[0]
    assert record[:5] == ["Chr1", "100", "SNP_31", "A", "C"]
    assert record[9:] == ["0/0", "0/1", "1/1"]


def test_metadata_low_confidence_mapping_is_ignored(tmp_path):
    mapping = {7: {"mapping_status": "low_confidence",
                   "reference_sequence": "Chr1", "genome_position": 5,
                   "reference_base": "A",
                   "allele_complement_required": False}}
    calls = {"SNP_7": {"s1": "AA", "s2": "AG"}}
    metadata = snpio.vcf_locus_metadata_from_mapping(["SNP_7"], calls, mapping)
    assert metadata == {}
    lines = _write(tmp_path, ["s1", "s2"], ["SNP_7"], calls,
                   locus_metadata=metadata)
    record = _records(lines)[0]
    assert record[0] == snpio.VCF_FALLBACK_CONTIG


def test_multiallelic_locus(tmp_path):
    calls = {"SNP_10": {"s1": "AA", "s2": "AG", "s3": "GG", "s4": "GT"}}
    lines = _write(tmp_path, ["s1", "s2", "s3", "s4"], ["SNP_10"], calls)
    record = _records(lines)[0]
    assert record[3] == "A"
    assert record[4] == "G,T"
    # AF: G = 4/8, T = 1/8
    assert record[7] == "NS=4;AF=0.500000,0.125000"
    assert record[9:] == ["0/0", "0/1", "1/1", "1/2"]


def test_reference_line(tmp_path):
    calls = {"NC_0:9A>G": {"s1": "AA", "s2": "AG"}}
    lines = _write(tmp_path, ["s1", "s2"], ["NC_0:9A>G"], calls,
                   reference_name="GCF_001858045.2_O_niloticus")
    assert "##reference=GCF_001858045.2_O_niloticus" in lines
    header = [line for line in lines if line.startswith("#CHROM")][0]
    assert header.split("\t")[9:] == ["s1", "s2"]


def test_missing_locus_data_skipped(tmp_path):
    lines = _write(tmp_path, ["s1"], ["SNP_1"], {"SNP_1": {"s1": "--"}})
    assert _records(lines) == []
