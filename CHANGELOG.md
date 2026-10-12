# Changelog

All notable changes to ABI-SNP Analyzer will be documented in this file.

## [1.1.0] - 2026-10-11

### Added
- **VCF 4.2 export** (reviewer request / Major comment on data export).
  "Export SNP Formats" now also writes `.vcf`: one sample column per
  individual, unphased `0/0` / `0/1` / `1/1` genotypes (`./.` = missing) and
  `NS` plus one `AF` value per ALT allele in INFO.  CHROM/POS/REF/ALT come
  from the unique BLASTN mapping (minus-strand loci are complemented onto the
  reference-forward strand); reference-forward SNP names such as
  `NC_031982.2:26437556C>T` are parsed directly, and loci without a unique
  mapping fall back to the virtual `amplicon` contig with POS = consensus
  coordinate (declared in the VCF header).  VCF is export-only and is not
  auto-detected by "Load Existing Data".
- **Version 1.1.0.**  Executable name, Windows version resource,
  `CITATION.cff` and documentation updated from 1.0.0.
- **Genotype-call validation tool** (`heterozygote_validation.py`, Major 3b).
  Compares the software output against a manually curated truth table and
  reports heterozygote sensitivity/precision/F1 with Wilson 95% intervals, the
  full confusion matrix, and every discordant call classified as missed
  heterozygote / false heterozygote / allele mismatch / no-call. Accepts any
  exported genotype format as the software side and wide- or long-format truth
  tables (.xlsx/.csv, plain alleles or IUPAC codes).
- **Linked-SNP dual-mode analysis** (Major 4v). SNPs from the same amplicon
  are usually in strong LD; treating them as independent loci biases
  multi-locus estimates. `linkage_analysis.py` adds:
  - pairwise r² (unphased dosage estimator, Hill & Weir 1994 / PLINK `--r2`),
    D' (EM haplotype frequencies, Excoffier & Slatkin 1995) and a permutation
    p-value for every SNP pair, **estimated within each population** (pooling
    diverged populations inflates LD through the Wahlburg effect);
  - single-linkage clustering of SNPs into linkage blocks at a user-set r²
    threshold (default 0.80) with one representative SNP per block (lowest
    missing rate, then highest MAF);
  - a new analysis-mode dialog: *All SNPs*, *One representative SNP per
    linkage block* (default) or *Run both modes and compare*;
  - `linkage_report.xlsx` (Pairwise_LD / LD_Summary / Linkage_Blocks /
    Method) and a "Linkage / Modes" tab in the results window with the
    side-by-side mode comparison.
- **hierfstat export** (`.txt`): header row `pop <loci...>`, one row per
  individual, `A=1, T=2, C=3, G=4`, `00` = missing. Read in R with
  `read.table(file, header = TRUE)` and consumed by `hierfstat::basic.stats()`.
- **GenAlEx export** (`.csv`): GenAlEx 6.5 diploid codominant layout
  (3 header rows, 2 allele columns per locus, `0` = missing); opens directly
  in Excel and is read by `poppr::read.genalex()`.
- **Generic genotype import** — "Load Existing Data" now accepts every format
  written by "Export SNP Formats" (GENEPOP, Arlequin `.arp`, hierfstat,
  GenAlEx, FSTAT, POPGENE32, STRUCTURE) in addition to `genotype_matrix.xlsx`.
  The format is auto-detected, with a manual selection dialog as fallback,
  and population structure is restored so imported data feeds the
  population-genetics module directly.
- `snp_format_io.py`: standalone, GUI-agnostic import/export module with
  round-trip tests (`tests/test_snp_format_io.py`,
  `tests/test_gui_import_roundtrip.py`).

### Changed
- All export writers now share one implementation (`snp_format_io`), so an
  exported file is guaranteed to be re-importable.
- GENEPOP: locus names are no longer stripped of variant annotations
  (`SNP01_A>C` stays intact); population names and the allele coding are
  carried in the title line for loss-less re-import.
- Arlequin: locus names are carried in the `Title` tag (.arp has no locus
  field) and are restored on import.
- Files carrying >= 2 populations enable "Statistical Analysis" immediately
  after import, without visiting the Population Manager first.

### Fixed
- FSTAT missing genotypes were written as a single `0`, an invalid two-allele
  code; they are now written as `00`.
- GenAlEx/CSV headers containing commas are now parsed with a real CSV reader
  (a quoted title cell previously corrupted allele-code recovery).

## [1.0.0] — 2026-07-21

### Added
- End-to-end GUI workflow: ABIF parsing → MSA → genotyping → population genetics
- Multi-sample batch processing of `.ab1` Sanger sequencing files
- Base-quality filtering for heterozygote calling (min base quality = Q10)
- Trace plot (PNG) generation for quality control
- IUPAC consensus sequence assembly
- MUSCLE v5 multiple-sequence alignment integration
- BLASTN reference-genome mapping (optional)
- Diploid genotype matrix generation with chromosome coordinates
- Population genetic statistics: Ho, He, PIC, Hardy-Weinberg exact test,
  F-statistics (WC84), AMOVA (Excoffier 1992)
- Export formats: Arlequin (.arp), GENEPOP, FSTAT, STRUCTURE, POPGENE32 (.pop)
- MIT License
