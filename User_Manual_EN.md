# ABI-SNP Analyzer v1.1.0 — User Manual

---

## Table of Contents

1. [Overview](#1-overview)
2. [System Requirements](#2-system-requirements)
3. [Installation & Launch](#3-installation--launch)
4. [GUI Layout](#4-gui-layout)
5. [Parameter Reference](#5-parameter-reference)
6. [Workflow](#6-workflow)
7. [Output Files](#7-output-files)
8. [Population Genetic Analysis](#8-population-genetic-analysis)
9. [Data Export Formats](#9-data-export-formats)
10. [Walkthrough Example](#10-walkthrough-example)
11. [FAQ](#11-faq)

---

## 1. Overview

ABI-SNP Analyzer is a free, open-source desktop application for automated SNP
genotyping and population genetic analysis directly from ABI-format Sanger
sequencing chromatogram files (.ab1). It provides an end-to-end graphical
workflow that integrates ABIF file parsing, heterozygote detection, multiple
sequence alignment, SNP calling, diploid genotype matrix generation, optional
reference genome anchoring via BLASTN, and built-in population genetic analysis
— all within a single graphical interface.

**Key functional modules:**

| Module | Function |
|---|---|
| ABIF Parser | Reads .ab1 binary format, extracts base sequences, quality scores, and four-channel fluorescence signals |
| Heterozygote Detection | Multi-channel signal ratio analysis to identify heterozygous sites and generate IUPAC ambiguity codes |
| Sequence Alignment & SNP Calling | MUSCLE alignment, similarity filtering, SNP discovery, and statistics |
| Reference Genome Mapping | BLASTN alignment of the consensus sequence to a reference genome for chromosome coordinates |
| Genotype Matrix | Sample-by-locus diploid genotype matrix in Excel format |
| Population Genetics | Population assignment, genetic diversity statistics, HWE exact tests, F-statistics, AMOVA |
| Data Export | Supports Genepop, Arlequin, hierfstat, GenAlEx, STRUCTURE, FSTAT and POPGENE32 — all round-trip importable — plus VCF 4.2 export (new in v1.1.0, export-only) |

---

## 2. System Requirements

| Item | Minimum | Recommended |
|---|---|---|
| OS | Windows 7 (64-bit) | Windows 10/11 (64-bit) |
| RAM | 4 GB | 8 GB or more |
| Disk | 500 MB | 1 GB+ (excluding reference genomes) |
| External tools | MUSCLE (v3.8.31 or v5) | BLAST+ 2.17.0+ (optional, for genome coordinate mapping) |

**External dependencies (download separately):**
- **MUSCLE** (v3.8.31 or v5, auto-detected): Multiple sequence aligner — https://www.drive5.com/muscle/
- **BLAST+ 2.17.0+**: Reference genome alignment (optional) — https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

> Tip: Place `muscle.exe` (and optionally `blastn.exe`) in the same directory
> as the application, or specify their paths via the GUI.

---

## 3. Installation & Launch

### 3.1 Standalone Executable (Recommended)

1. Download `ABI_SNP_Analyzer_v1.1.0.exe` from
   [GitHub Releases](https://github.com/xiacy111111/abi-snp-analyzer/releases).
2. Place `muscle.exe` (and optionally `blastn.exe`) in the same folder.
3. Double-click `ABI_SNP_Analyzer_v1.1.0.exe` to launch.

### 3.2 Run from Source

```bash
git clone https://github.com/xiacy111111/abi-snp-analyzer.git
cd abi-snp-analyzer && git checkout v1.1.0
python -m venv .venv && .venv\Scripts\activate
pip install -r requirements.txt
python abi_to_snp_analyzer_v1.1.0.py
```

---

## 4. GUI Layout

The interface is organized into four functional zones:

```
┌─────────────────────────────────────────────────────────┐
│  1. File Selection                                      │
│  ├─ ABI Files: Select .ab1 files or folder              │
│  ├─ MUSCLE Path: Specify muscle.exe location            │
│  ├─ Reference Genome: Reference FASTA file (optional)   │
│  ├─ BLASTN Path: Specify blastn.exe location (optional) │
│  └─ Output Directory: Select result output folder       │
├─────────────────────────────────────────────────────────┤
│  2. Filter Threshold Settings                           │
│  ├─ Sequence Similarity / MAF Threshold                 │
│  ├─ Min Observed Het / Signal Ratio                     │
│  └─ Consensus Quality / Min Base Quality                │
├─────────────────────────────────────────────────────────┤
│  3. Run Controls                                        │
│  ├─ [Start Analysis] [Stop] [Load Existing Data]       │
│  └─ [Population Manager] [Statistical Analysis]        │
│     [Export SNP Formats]                               │
├─────────────────────────────────────────────────────────┤
│  4. Run Log                                             │
│     Real-time display of progress, parameters, and      │
│     error messages                                      │
└─────────────────────────────────────────────────────────┘
```

---

## 5. Parameter Reference

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Sequence Similarity** | 0.90 | 0–1 | Minimum average pairwise identity. Sequences below this threshold are excluded. Recommended 0.85–0.95 for same-amplicon samples |
| **MAF Threshold** | 0.10 | 0–1 | Minimum minor allele frequency. Only SNPs with MAF >= this value are retained |
| **Min Observed Het** | 0.10 | 0–1 | Minimum observed heterozygosity per locus. Set to 0 to disable filtering |
| **Signal Ratio** | 0.30 | 0–1 | Heterozygote detection threshold (secondary/primary peak). Lower = more sensitive; higher = more specific |
| **Consensus Quality** | 0.20 | 0–1 | Consensus calling threshold. Positions with consensus proportion below this are set to N |
| **Min Base Quality** | 5 | 0–100 | Phred quality score threshold. Bases with quality below this are excluded from genotyping |

---

## 6. Workflow

### 6.1 Basic Pipeline: SNP Genotyping

**Step 1 — Prepare input**
Place all .ab1 files in a single folder. Ensure unique file names and that all
files come from the same target amplicon region.

**Step 2 — Configure paths**
- Click **Select Files/Folder** and choose the folder (or select individual .ab1 files)
- Click **Browse** next to **MUSCLE Path** and select `muscle.exe`
- (Optional) Click **Browse** next to **Reference Genome** and select a FASTA file
- (Optional, required if reference genome is selected) Click **Browse** next to **BLASTN Path** and select `blastn.exe`
- Click **Select** next to **Output Directory** and choose the output folder

**Step 3 — Set parameters**
Adjust the six parameters according to your experiment (see [Section 5](#5-parameter-reference)).
Default values are suitable for most routine analyses.

**Step 4 — Run analysis**
Click **Start Analysis**. The software automatically performs these sub-steps:

1. **ABIF parsing & heterozygote detection**: Read each .ab1 file, extract
   chromatogram signals, identify heterozygous sites, and generate
   IUPAC-encoded FASTA sequences
2. **Sequence similarity filtering**: Remove sequences that are too divergent
   from the group
3. **MUSCLE multiple sequence alignment**: Align all retained sequences
4. **Consensus sequence & trace plots**: Generate the consensus sequence;
   produce chromatogram trace plots and heterozygote reports for samples
   with heterozygous calls
5. **SNP discovery & filtering**: Identify polymorphic columns and filter by
   MAF and heterozygosity thresholds
6. **(Optional) BLASTN genome mapping**: Align consensus to reference genome
7. **Genotype matrix generation**: Output the sample-by-locus diploid genotype
   matrix (Excel)

**Step 5 — Review results**
When analysis completes, an "Analysis completed" dialog appears. Inspect the
generated files in the output directory (see [Section 7](#7-output-files)).

### 6.2 Resuming: Load Existing Data

If genotyping was previously completed, you can reload results for population
genetics analysis:

1. Click **Load Existing Data**
2. Select a data file. Besides the native `genotype_matrix.xlsx`, **any file
   written by "Export SNP Formats" can be read back**: GENEPOP (.txt),
   Arlequin (.arp), hierfstat (.txt), GenAlEx (.csv), FSTAT (.dat),
   POPGENE32 (.pop), STRUCTURE (.txt)
3. The format is detected automatically; if the content is ambiguous a dialog
   lets you pick the format explicitly
4. A preview shows the detected format, sample/loci counts and populations —
   confirm to load
5. After loading, **Population Manager** and **Export SNP Formats** become
   active. **When the file itself carries ≥2 populations (Arlequin, GENEPOP,
   hierfstat, GenAlEx, FSTAT and POPGENE32 all do), the Statistical Analysis
   button is enabled as well, so you can go straight to the analysis.**

**Round-trip fidelity:**

| Format | Individual IDs | Locus names | Population names |
|---|---|---|---|
| `genotype_matrix.xlsx` | kept | kept | n/a |
| GENEPOP (.txt) | kept | kept | kept (carried in the title tag) |
| Arlequin (.arp) | kept | kept (carried in the `Title` tag) | kept |
| hierfstat (.txt) | rebuilt as `<pop>_Ind####` | kept | rebuilt as Pop1…PopN |
| GenAlEx (.csv) | kept | kept | kept |
| FSTAT (.dat) | rebuilt | kept | rebuilt as Pop1…PopN |
| POPGENE32 (.pop) | rebuilt | kept | rebuilt as Pop1…PopN |
| STRUCTURE (.txt) | kept | rebuilt as Locus_1…N | rebuilt as Pop1…PopN |

> hierfstat and FSTAT encode populations as bare numbers and store no
> individual labels; POPGENE32 and STRUCTURE have no name fields either. This
> is a property of those formats — genotype calls are always preserved.

**Allele coding:** every numeric format uses `1=A, 2=T, 3=C, 4=G, 0=missing`.
The mapping is written into the GENEPOP title line, the Arlequin `Title` and
the GenAlEx dataset title, and is read back on import instead of being assumed.

---

## 7. Output Files

After analysis, the output directory contains:

| File | Format | Description |
|---|---|---|
| `{sample}_sequence_with_iupac.fasta` | FASTA | IUPAC-encoded sequence per sample |
| `{sample}_trace_plot_combined.png` | PNG | Chromatogram trace with heterozygote annotations (samples with het calls only) |
| `{sample}_heterozygote_report.xlsx` | Excel | Per-sample heterozygote details (positions, signal ratios) |
| `combined_all_forward_sequences_with_iupac.fasta` | FASTA | All IUPAC sequences combined |
| `alignment.aln` | text | MUSCLE multiple sequence alignment |
| `consensus_sequence.fasta` | FASTA | Consensus sequence |
| `consensus_vs_reference.aln` | text | BLASTN alignment to reference genome (optional) |
| `genotype_matrix.xlsx` | Excel | Sample-by-SNP diploid genotype matrix (sheet: Genotype_Matrix) |
| `genotype_matrix_snp_analysis.xlsx` | Excel | Per-SNP statistics: MAF, heterozygosity, genome coordinates (sheet: SNP_Analysis_Results) |
| `analysis_metadata.json` | JSON | Reproducibility metadata: parameters, input files, dependency versions |

**Genotype encoding:**
- Homozygous: `AA`, `TT`, `CC`, `GG`
- Heterozygous: `AG`, `CT`, `AC`, etc.
- Missing: `--`

**SNP naming:**
- Without reference genome: `SNP_586` (position on consensus sequence)
- With reference genome: `chr1:1234567G>A` (chromosome:position reference>alternate)

---

## 8. Population Genetic Analysis

### 8.1 Population Manager

1. After genotyping completes, click **Population Manager**
2. In the dialog:
   - The left panel shows all samples
   - Click **Add Population** to create a new population and enter its name
   - Select samples from the left and click **>** to assign them to the population on the right
   - Each sample belongs to exactly one population
   - At least two populations must be defined for statistical analysis

### 8.2 Statistical Analysis

After population assignment, click **Statistical Analysis**. The software
first asks how linked SNPs should be handled (§8.3) and then computes:

| Analysis | Output |
|---|---|
| Data overview | Sample and locus counts per population |
| Genetic diversity | Observed heterozygosity (Ho), expected heterozygosity (He), polymorphic information content (PIC) |
| HWE exact test | Per-locus, per-population Hardy-Weinberg equilibrium exact test p-values [Guo & Thompson 1992; Wigginton et al. 2005] |
| WC84 F-statistics | Fis (inbreeding coefficient), Fst (differentiation), Fit (total) with 1,000 permutation tests [Weir & Cockerham 1984] |
| Null allele frequency | Brookfield estimator 3 [Brookfield 1996] |
| AMOVA | Analysis of molecular variance with Phi-statistic and significance p-value [Excoffier et al. 1992] |
| Genotype frequency chi-square | Inter-population genotype frequency difference test with Benjamini-Hochberg FDR-corrected p-values (Bonferroni also reported) |
| Interactive visualization | Per-locus genotype frequency bar charts (300 DPI PNG/PDF export) |
| Summary report | All statistical results in Excel format with raw and corrected p-values |

### 8.3 Linked SNPs: dual-mode handling

SNPs called from the same amplicon are physically linked and usually in strong
linkage disequilibrium (LD). Treating them as independent loci inflates the
degrees of freedom of multi-locus tests, underestimates the variance and
violates the independent-loci assumption of software such as STRUCTURE. The
analysis therefore starts with a mode dialog:

| Mode | Behaviour | When to use |
|---|---|---|
| **All SNPs** | Keep every SNP (independent-loci assumption) | Loci come from different amplicons |
| **One representative SNP per linkage block** (default) | SNPs with r² ≥ threshold are clustered into blocks; only the most reliable SNP of each block enters the statistics | Several SNPs share an amplicon |
| **Run both modes and compare** | Runs the pipeline twice and shows a side-by-side comparison | Reporting the effect of the linkage correction |

**Algorithm and output**

1. **Pairwise LD**: r² is the squared Pearson correlation of minor-allele
   dosages (0/1/2) — the standard unphased estimator (Hill & Weir 1994;
   PLINK `--r2`, Purcell et al. 2007). D' is computed from two-locus
   haplotype frequencies estimated by EM (Lewontin 1964; Excoffier & Slatkin
   1995). The p-value is a permutation test (default 999) that shuffles the
   genotypes of the second locus.
2. **LD is estimated within each population** and summarised across
   populations (the maximum r² drives the clustering). Pooling diverged
   populations inflates LD through the Wahlburg effect, so it is avoided.
3. **Blocks**: single-linkage clustering on r² ≥ threshold (default 0.80,
   adjustable 0.50–1.00 in the dialog).
4. **Representative**: per block, the SNP with the lowest missing rate, then
   the highest MAF, then the original order is kept.
5. **`linkage_report.xlsx`** with four sheets: `Pairwise_LD` (per-population
   r², D', p), `LD_Summary` (cross-population summary), `Linkage_Blocks`
   (blocks and representatives), `Method` (estimators and references).
6. The results window gains a **"Linkage / Modes"** tab with the block table,
   the pairwise LD table and the mode-comparison table (loci analysed, mean
   Ho/He/PIC/Fst, AMOVA Phi-ST and p).

> D' is strongly upward-biased in small samples — define blocks with r²; D' is
> reported for comparability with Arlequin/GENEPOP output.

---

## 9. Data Export Formats

SNP genotype matrices can be exported to eight formats: seven common
population-genetics formats, all read back with **Load Existing Data**, plus
VCF 4.2 (new in v1.1.0, export-only, for downstream tools such as
bcftools / vcftools / tabix).

| Format | Extension | Encoding | Target software |
|---|---|---|---|
| **Genepop** | .txt | 2-digit integer (A=01, T=02, C=03, G=04); 0000=missing | Genepop, Fstat |
| **Arlequin** | .arp | Arlequin 3.5 project file, STANDARD data + [Structure] block | Arlequin 3.5 |
| **hierfstat** | .txt | Header `pop` + locus names, then `pop_id genotype` rows (A=1,T=2,C=3,G=4; 00=missing) | R hierfstat (`read.table(header=TRUE)`) |
| **GenAlEx** | .csv | 3 header rows + 2 allele columns per locus (0=missing) | GenAlEx 6.5, poppr/adegenet |
| **STRUCTURE** | .txt | 2 columns per SNP; -9=missing | STRUCTURE |
| **FSTAT** | .dat | Integer allele codes with locus headers | FSTAT 2.9.x |
| **POPGENE32** | .pop | 2-letter genotype codes; `..`=missing | POPGENE 1.32 |
| **VCF 4.2** | .vcf | One genotype column per sample (0/0, 0/1, 1/1; ./. = missing); CHROM/POS from the unique BLASTN mapping (minus-strand loci complemented to the reference-forward strand); loci without a unique mapping fall back to the virtual contig `amplicon` | bcftools, vcftools, tabix |

hierfstat example (tab-delimited):

```
pop     SNP01   SNP02   SNP03
1       11      12      00
1       12      00      11
2       11      12      00
```

GenAlEx example (comma-delimited):

```
3,4,2,2,2
"ABI-SNP Analyzer v1.1.0 - SNP genotypes | ALLELE_CODE:1=A,2=T,3=C,4=G",,,PopA,PopB
Sample,Pop,SNP01,,SNP02,,SNP03,
S01,PopA,1,1,1,2,0,0
S02,PopA,1,2,0,0,1,1
S03,PopB,1,1,1,2,0,0
S04,PopB,1,2,0,0,1,1
```

**Export steps:**
1. After genotyping or loading data, click **Export SNP Formats**
2. Select the target format from the dialog
3. Check "Include population assignment info" to include population grouping
4. Click **Export** and choose the save location

---

## 10. Walkthrough Example

This example uses Sanger sequencing data from a salt tolerance-associated QTL
region validation study in Red tilapia (*Oreochromis* spp.).

### 10.1 Data Preparation
- 88 .ab1 Sanger sequencing files (44 salt-tolerant + 44 salt-sensitive individuals)
- MUSCLE v3.8.31 (`muscle.exe`)
- Nile tilapia reference genome FASTA (~1 GB)
- NCBI BLAST+ (`blastn.exe`)

### 10.2 Parameter Settings
| Parameter | Value |
|---|---|
| Sequence Similarity | 0.90 |
| MAF Threshold | 0.10 |
| Min Observed Het | 0.10 |
| Signal Ratio | 0.30 |
| Consensus Quality | 0.20 |
| Min Base Quality | 5 |

### 10.3 Execution
1. Select the folder containing 88 .ab1 files
2. Configure MUSCLE, reference genome, and BLASTN paths
3. Choose an output directory
4. Set parameters as above
5. Click **Start Analysis** — processing completes in approximately 8 minutes

### 10.4 Results
- Consensus sequence length: 861 bp
- SNPs detected: 10
- Genotype data points: 880 (88 samples × 10 loci)
- Genotyping accuracy: 99.8% (878/880 concordant with manual curation)

### 10.5 Population Genetic Analysis
1. Click **Population Manager**, create two populations:
   - "Salt-Tolerant" (44 individuals)
   - "Salt-Sensitive" (44 individuals)
2. Click **Statistical Analysis** to view population genetic statistics
3. Click **Export SNP Formats** to export data for downstream software

---

## 11. FAQ

**Q: The software reports "MUSCLE file not found." What should I do?**
A: Ensure the `muscle.exe` path is correct. Download MUSCLE (v3.8.31 or v5) from
https://www.drive5.com/muscle/.

**Q: How long does reference genome mapping take?**
A: It depends on genome size. With the Nile tilapia genome (~1 GB),
BLASTN alignment typically takes a few seconds to tens of seconds.

**Q: How long does batch processing of 400 samples take?**
A: Approximately 37 minutes with BLASTN genome mapping enabled
(tested on i7-10700F, 16 GB RAM).

**Q: How accurate is heterozygote detection?**
A: At the default signal ratio threshold of 0.30, concordance with manual
curation is approximately 99.8%. Adjust the Signal Ratio parameter to change
detection sensitivity.

**Q: Does the software support polyploid species?**
A: The current version only supports diploid genotyping and population
genetics analysis.

**Q: Where are the output files?**
A: All output files are located in the Output Directory selected before
starting the analysis. Sort by modification time to find the latest files.

**Q: How do I stop a running analysis?**
A: Click the **Stop** button to safely abort the current analysis.

**Q: Can I share analysis settings across computers?**
A: Each run generates an `analysis_metadata.json` file that records all
parameters and dependency versions for reproducibility.

---

## Contact & Support

- **GitHub**: https://github.com/xiacy111111/abi-snp-analyzer
- **Issue Tracker**: https://github.com/xiacy111111/abi-snp-analyzer/issues
- **Email**: pbning@xidian.edu.cn

---

*ABI-SNP Analyzer v1.1.0 | MIT License | Copyright (c) 2026 Chengyu Xia, Wenhui Chu, Pengbo Ning*
