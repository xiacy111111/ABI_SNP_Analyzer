# ABI-SNP Analyzer v1.0.0

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%2064--bit-lightgrey.svg)](#system-requirements)

**ABI-SNP Analyzer** is a free, open-source desktop application that automates
SNP genotyping and population genetic analysis directly from ABI-format Sanger
sequencing chromatograms (.ab1 files). It provides an end-to-end graphical
workflow: ABIF file parsing, heterozygote detection via multi-channel signal
ratio analysis, MUSCLE-based multiple sequence alignment, SNP calling, diploid
genotype matrix generation, optional BLASTN-based reference genome anchoring,
and built-in population genetic analyses.

This repository accompanies a manuscript submitted to *SoftwareX*.
The source code is released under the MIT License.

---

## Table of Contents

- [Overview](#overview)
- [Motivation](#motivation)
- [Features](#features)
- [Software Architecture](#software-architecture)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Parameters](#parameters)
- [Input and Output](#input-and-output)
- [Population Genetic Analyses](#population-genetic-analyses)
- [Export Formats](#export-formats)
- [Illustrative Example](#illustrative-example)
- [Performance](#performance)
- [Limitations](#limitations)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Citation](#citation)
- [Code Metadata](#code-metadata)
- [License](#license)
- [Third-Party Software](#third-party-software)
- [References](#references)

---

## Overview

SNPs are the most abundant class of genetic markers in eukaryotic genomes and
serve as the cornerstone of modern animal and plant genetics, genomics, and
molecular breeding. Applications such as marker-assisted selection (MAS),
genomic selection (GS), and genome-wide association studies (GWAS) all rely on
accurate and scalable SNP genotyping [1, 2].

Despite the widespread use of Sanger sequencing for SNP genotyping, the
analytical workflow from raw .ab1 files to a finalized sample-by-locus genotype
matrix remains fragmented and technically demanding. The typical workflow
involves parsing binary ABIF files, visually inspecting chromatogram traces for
heterozygous sites, performing multiple sequence alignment, manually
constructing a genotype matrix, and computing population genetic statistics —
all using disconnected tools.

ABI-SNP Analyzer integrates these operations into a single graphical application:

1. Parse ABI/ABIF chromatogram files with a custom binary parser.
2. Detect heterozygous sites using an adjustable peak-signal ratio.
3. Encode mixed calls as IUPAC ambiguity codes.
4. Filter anomalous sequences and run MUSCLE multiple-sequence alignment.
5. Discover and filter SNP loci (MAF, heterozygosity thresholds).
6. Generate diploid genotype matrices and per-locus statistics (Excel).
7. Optionally map SNPs to reference genome coordinates via BLASTN.
8. Assign samples to populations and compute genetic statistics.
9. Export data in Genepop, STRUCTURE, FSTAT, Arlequin, and POPGENE32 formats.

## Motivation

Existing tools address parts of this pipeline but fall short of providing a
fully integrated, user-friendly solution. Commercial platforms such as Geneious
Prime are expensive and lack specialized batch processing for breeding-oriented
genotype matrix generation. Mutation Surveyor targets clinical diagnostics
rather than population genetics. Open-source tools such as Phred/Phrap [5] and
SAGPEK [6] require command-line expertise. Furthermore, no existing Sanger
sequencing analysis tool provides a built-in, GUI-driven population genetics
module that operates directly on the generated genotype matrix.

ABI-SNP Analyzer fills this gap by offering a complete, open-source pipeline
from raw chromatogram files to population genetic statistics in a single
graphical interface, targeting researchers in animal and plant genetics,
conservation biology, and molecular breeding.

## Features

### ABIF File Parsing
- Custom Python ABIF parser — no external commercial libraries required.
- Reads binary ABIF format (big-endian multi-byte integers, FWO_ channel-order
  mapping).
- Extracts base sequences, quality scores, positions, and four-channel
  fluorescence signals (A/C/G/T).

### Heterozygote Detection
- Adaptive base-interval peak detection with configurable signal-ratio threshold
  (default 0.30).
- IUPAC ambiguity encoding for valid two-base combinations (R, Y, S, W, K, M).
- Per-sample chromatogram trace plots (PNG) and heterozygote reports (Excel).

### Sequence Alignment and SNP Genotyping
- Pairwise sequence similarity filtering (default 0.90).
- MUSCLE [13] multiple-sequence alignment.
- SNP discovery with MAF and observed heterozygosity filtering.
- Diploid genotype matrix output (Excel) with per-locus QC statistics.

### Reference Genome Mapping
- BLASTN+ [14] in `-subject` mode — no pre-built database required.
- Positive/negative strand coordinate handling.
- Coordinate-aware SNP identifiers (e.g., `chr1:1234567G>A`).

### Population Genetics
- Graphical population assignment via tree-view dialog.
- Summary statistics: observed heterozygosity (Ho), expected heterozygosity
  (He), polymorphic information content (PIC) [15, 16].
- Hardy-Weinberg equilibrium exact tests [17, 18].
- WC84 F-statistics (Fis, Fst, Fit) with permutation testing [19].
- AMOVA [21] with Phi-statistic.
- Genotype frequency chi-square tests with Bonferroni correction.
- Interactive genotype frequency bar charts (Matplotlib).
- Multi-format export: Genepop [8], STRUCTURE, FSTAT [10], Arlequin [9],
  POPGENE32 [11].

## Software Architecture

ABI-SNP Analyzer follows a three-layer architecture:

1. **Data Layer** — ABIF file parsing, fluorescence signal extraction, and
   sequence quality assessment.
2. **Processing Layer** — core algorithms for heterozygote detection, sequence
   alignment, SNP identification, genome coordinate mapping, and population
   genetics analysis.
3. **Presentation Layer** — Tkinter-based GUI with parameter controls, real-time
   log output, progress reporting, and result visualization.

| Module | Responsibility |
|---|---|
| `src/abi_to_snp_analyzer_v1.0.0.py` | GUI, ABIF parser, pipeline orchestration |
| `population_genetics.py` | Statistical calculations and export formatting |
| `chromatogram.py` | Chromatogram processing and visualization helpers |

**Technology stack:** Python 3, Tkinter, NumPy, SciPy, Pandas, Matplotlib,
Biopython, openpyxl, Pillow. External tools: MUSCLE v5, BLAST+ 2.12+.

## System Requirements

| Requirement | Minimum | Recommended |
|---|---|---|
| Operating system | Windows 7–11, 64-bit | Windows 10/11, 64-bit |
| RAM | 4 GB | 8 GB+ |
| Disk space | 500 MB | 1 GB+ (excluding reference genomes) |
| Python | Not required (standalone .exe) | 3.8+ for running from source |

**External dependencies (must be installed separately):**
- [MUSCLE v5](https://www.drive5.com/muscle/) — multiple sequence alignment
- [BLAST+ 2.12+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) — optional, for reference genome mapping

## Installation

### Standalone Executable

1. Download the latest `ABI_SNP_Analyzer_v1.0.0.exe` from
   [Releases](https://github.com/xiacy111111/ABI_SNP_Analyzer/releases).
2. Place `muscle.exe` (and optionally `blastn.exe`) in the same directory
   as the executable, or ensure they are on the system PATH.
3. Double-click `ABI_SNP_Analyzer_v1.0.0.exe` to launch.

### From Source

```bash
git clone https://github.com/xiacy111111/ABI_SNP_Analyzer.git
cd abi-snp-analyzer && git checkout v1.0.0
python -m venv .venv && .venv\Scripts\activate
pip install -r requirements.txt
python src/abi_to_snp_analyzer_v1.0.0.py
```

## Quick Start

1. Launch the application.
2. Click **Select Files/Folder** and choose a folder of `.ab1` files.
3. Select `muscle.exe`.
4. Choose an output directory.
5. (Optional) Select a reference genome FASTA and `blastn.exe`.
6. Click **Start Analysis**.
7. Review trace plots (PNG) and heterozygote reports (Excel).
8. Open `genotype_matrix.xlsx` for the sample-by-SNP genotype table.
9. Click **Population Manager**, assign samples to at least two populations,
   then click **Statistical Analysis** for population genetic statistics.
10. Use **Export SNP Formats** to export for Genepop / STRUCTURE / FSTAT /
    Arlequin / POPGENE32.

To resume a previous session, click **Load Existing Data** and select
`genotype_matrix.xlsx`.

## Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| Sequence Similarity | 0.90 | 0–1 | Minimum average pairwise identity; excludes divergent sequences |
| MAF Threshold | 0.05 | 0–1 | Minimum minor allele frequency for SNP inclusion |
| Min Observed Heterozygosity | 0.10 | 0–1 | Minimum observed heterozygosity per locus |
| Signal Ratio | 0.30 | 0–1 | Heterozygote detection threshold; lower = more sensitive |
| Consensus Quality | 0.20 | 0–1 | Bases below this consensus proportion are called as N |
| Min Base Quality | 5 | 0–100 | Phred quality score threshold for base inclusion in genotyping |

## Input and Output

### Input

**Required:** `.ab1` files (same target amplicon region), unique file names,
`muscle.exe`.

**Optional:** Reference genome FASTA, `blastn.exe`. All data is processed
locally.

### Output Files

| File | Format | Description |
|---|---|---|
| `{sample}_trace_plot_combined.png` | PNG | Chromatogram trace with heterozygote annotations |
| `{sample}_heterozygote_report.xlsx` | Excel | Per-sample heterozygote positions and signal ratios |
| `{sample}_sequence_with_iupac.fasta` | FASTA | IUPAC-encoded sequence per sample |
| `combined_all_forward_sequences_with_iupac.fasta` | FASTA | All IUPAC sequences merged |
| `alignment.aln` | text | MUSCLE multiple sequence alignment |
| `consensus_sequence.fasta` | FASTA | Consensus sequence |
| `genotype_matrix.xlsx` | Excel | Sample × SNP diploid genotype matrix |
| `genotype_matrix_snp_analysis.xlsx` | Excel | Per-locus SNP statistics (MAF, coordinates, etc.) |
| `consensus_vs_reference.aln` | text | BLASTN alignment to reference genome (optional) |

Genotype encoding: `AA`, `AG`, `--` (missing). SNP identifiers:
`SNP_586` (consensus-position) or `chr1:1234567G>A` (genome coordinate).

## Population Genetic Analyses

| Analysis | Details |
|---|---|
| Summary statistics | Ho, He, PIC [15, 16] per population and per locus |
| HWE exact test | Per-locus Hardy-Weinberg equilibrium exact test [17, 18] |
| F-statistics | Fis, Fst, Fit with 999-permutation test (WC84) [19] |
| Null allele frequency | Brookfield estimator [20] |
| AMOVA | Pairwise-difference variance components, Phi-statistic, p-value [21] |
| Genotype frequency test | Chi-square test with Bonferroni correction [22] |
| Visualization | Interactive genotype frequency bar charts (300 DPI) |

## Export Formats

| Format | Extension | Details |
|---|---|---|
| Genepop | .txt | 3-digit integer encoding (A=1, T=2, C=3, G=4); 000 = missing |
| STRUCTURE | .txt | Two columns per SNP; -9 = missing |
| FSTAT | .dat | Integer allele encoding with locus headers |
| Arlequin | .arp | Excoffier 1992 AMOVA format; `[[Loci]]` + `[Structure]` blocks |
| POPGENE32 | .pop | Two-letter code; `..` = missing |

## Illustrative Example

We validated ABI-SNP Analyzer using Sanger sequencing data from a salt
tolerance-associated QTL region validation study in Red tilapia
(*Oreochromis* spp.) [23]. The dataset comprised 88 individuals
(44 salt-tolerant and 44 salt-sensitive).

**Analysis parameters:** sequence similarity = 0.9, signal ratio = 0.3,
MAF = 0.10, min observed heterozygosity = 0.10, consensus quality = 0.2,
min base quality = 5. The *O. niloticus* reference genome (~1 Gb) was used
for BLASTN mapping.

**Results:** The software completed SNP genotyping in 8 min 01 s, generating
an 861-bp consensus sequence with 10 detected SNP loci (880 genotype data
points). Genotype accuracy was **99.8%** (878/880 concordant calls) compared
with manually verified genotypes.

**Population analysis:** Salt-tolerant individuals exhibited slightly higher
genetic diversity (mean Ho = 0.4273, mean PIC = 0.3072) than salt-sensitive
individuals (mean He = 0.3207, mean PIC = 0.2460). Overall Fst = 0.0443.
AMOVA attributed 5.93% of variation to between-population differences
(Phi = 0.0593, p = 0.0030), consistent with independent Arlequin validation
(Fst = 0.04919, p = 0.00293). Two SNPs showed significant differentiation
at p < 0.05.

## Performance

Benchmarks on Windows 11 (Intel i7-10700F, 16 GB RAM) with the
*O. niloticus* genome (~1 Gb) as the BLASTN reference:

| Samples | Genotype detection | Population genetics | Total |
|---|---:|---:|---:|
| 10 | 2 min 28 s | 2 s | 2 min 30 s |
| 50 | 4 min 33 s | 16 s | 4 min 49 s |
| 100 | 9 min 0 s | 21 s | 9 min 21 s |
| 200 | 17 min 57 s | 1 min 10 s | 19 min 7 s |
| 400 | 35 min 59 s | 1 min 13 s | 37 min 12 s |

Per-sample genotype detection time decreased from 14.8 s/sample (10 samples)
to approximately 5.4 s/sample (50–400 samples) as initialization overhead was
amortized. Population genetics analysis time grows modestly with sample size.

## Limitations

- Requires MUSCLE as an external dependency for multiple sequence alignment.
- Population genetic analyses assume biallelic diploid markers.
- GUI currently runs on Windows 64-bit only.
- Processing time scales with sample count and reference genome size.
- Heterozygote detection via peak-signal ratio may be affected by dye blobs,
  pull-up peaks, and baseline noise.
- Automated genotyping results should be reviewed against chromatogram traces.

## Repository Structure

```
ABI_SNP_Analyzer/
  src/
    abi_to_snp_analyzer_v1.0.0.py # Main application (GUI + pipeline)
    population_genetics.py         # Population genetics statistics module
    chromatogram.py                # Chromatogram processing helpers
  requirements.txt                 # Python dependencies
  README.md                        # This file
  LICENSE                          # MIT License
  THIRD_PARTY_NOTICES.md           # Third-party software notices
  CHANGELOG.md                     # Release history
  CITATION.cff                     # Citation metadata (CFF format)
  CONTRIBUTING.md                  # Contribution guidelines
  version_info.txt                 # Windows EXE version metadata
  bin/                             # Bundled external tools (MUSCLE, BLASTN+)
  docx/                            # User manuals (EN/CN)
  examples/                        # Example .ab1 input files
  genome/                          # Example reference FASTA
  tests/                           # Unit tests (pytest)
    test_*.py
  .github/workflows/               # CI configuration
    tests.yml
  example/                         # Example .ab1 input files
```

> **Note:** External tools (MUSCLE, BLAST+) and reference genome FASTA files
> are not bundled — see [System Requirements](#system-requirements).

## Contributing

Bug reports, feature suggestions, documentation improvements, and pull requests
are welcome. Before submitting an issue, please:

1. Confirm the problem occurs with the latest release.
2. Provide OS version, software version, and parameter settings.
3. Include a minimal reproducible example if possible.

Use the [GitHub issue tracker](https://github.com/xiacy111111/ABI_SNP_Analyzer/issues)
for public reports. Email [25061300031@stu.xidian.edu.cn](mailto:25061300031@stu.xidian.edu.cn)
for private inquiries.

## Citation

If you use ABI-SNP Analyzer in your research, please cite:

> Xia, C. et al. (2026). ABI-SNP Analyzer: An open-source desktop
> application for automated SNP genotyping and population genetic analysis
> from Sanger sequencing chromatograms. *SoftwareX*.
> [DOI to be inserted upon acceptance]

BibTeX:
```bibtex
@article{xia2026abisnp,
  author = {Xia, Chengyu, Wenhui Chu, Pengbo Ning},
  title = {ABI-SNP Analyzer: An open-source desktop application for
           automated SNP genotyping and population genetic analysis
           from Sanger sequencing chromatograms},
  journal = {SoftwareX},
  year = {2026},
  note = {DOI to be inserted upon acceptance}
}
```

## Code Metadata

| | Description | Metadata |
|---|---|---|
| C1 | Current code version | v1.0.0 |
| C2 | Permanent link to repository | https://github.com/xiacy111111/ABI_SNP_Analyzer |
| C3 | Legal code license | MIT License |
| C4 | Code versioning system | git |
| C5 | Languages and tools | Python 3, Tkinter, NumPy, Pandas, Matplotlib, Biopython, MUSCLE, BLASTN+ |
| C6 | Requirements and dependencies | Windows 7–11 (64-bit); Python 3.8+; MUSCLE v5; BLAST+ 2.12+ |
| C7 | Documentation | https://github.com/xiacy111111/ABI_SNP_Analyzer |
| C8 | Support email | 25061300031@stu.xidian.edu.cn |

## License

ABI-SNP Analyzer is released under the **MIT License**.
Copyright (c) 2026 Chengyu Xia, Wenhui Chu, Pengbo Ning.

See [`LICENSE`](LICENSE) for the full license text.

## Third-Party Software

ABI-SNP Analyzer invokes separately distributed command-line tools at runtime:

- **MUSCLE v5** (https://www.drive5.com/muscle/) — multiple sequence alignment.
- **NCBI BLAST+ 2.17.0** (https://blast.ncbi.nlm.nih.gov/) — reference genome
  mapping.

These tools are not covered by the MIT License. Users must download them
independently and comply with their respective licenses. See
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for details.

## References

[1] Vignal A, et al. A review on SNP and other types of molecular markers. *Genet Sel Evol*. 2002;34(3):275–305.

[2] Davey JW, et al. Genome-wide genetic marker discovery and genotyping using next-generation sequencing. *Nat Rev Genet*. 2011;12(7):499–510.

[3] Sanger F, Nicklen S, Coulson AR. DNA sequencing with chain-terminating inhibitors. *Proc Natl Acad Sci USA*. 1977;74(12):5463–5467.

[4] Kircher M, Kelso J. High-throughput DNA sequencing — concepts and limitations. *Bioessays*. 2010;32(6):524–536.

[5] Ewing B, et al. Base-calling of automated sequencer traces using Phred. *Genome Res*. 1998;8(3):175–185.

[6] Wang J, et al. SAGPEK: fast and flexible approach to identify genotypes of Sanger sequencing data. *BMC Bioinformatics*. 2025;26:246.

[7] Suarez Menendez M, et al. PHFinder: assisted detection of point heteroplasmy in Sanger sequencing chromatograms. *PeerJ*. 2023;11:e16028.

[8] Rousset F. genepop'007. *Mol Ecol Resour*. 2008;8(1):103–106.

[9] Excoffier L, Lischer HEL. Arlequin suite ver 3.5. *Mol Ecol Resour*. 2010;10(3):564–567.

[10] Goudet J. FSTAT (version 1.2). *J Hered*. 1995;86(6):485–486.

[11] Yeh FC, et al. POPGENE 32, version 1.32. University of Alberta; 2000.

[12] Cock PJA, et al. Biopython. *Bioinformatics*. 2009;25(11):1422–1423.

[13] Edgar RC. MUSCLE. *Nucleic Acids Res*. 2004;32(5):1792–1797.

[14] Camacho C, et al. BLAST+. *BMC Bioinformatics*. 2009;10:421.

[15] Nei M. Estimation of average heterozygosity. *Genetics*. 1978;89(3):583–590.

[16] Botstein D, et al. Construction of a genetic linkage map in man using RFLPs. *Am J Hum Genet*. 1980;32(3):314–331.

[17] Guo SW, Thompson EA. Performing the exact test of HWE. *Biometrics*. 1992;48(2):361–372.

[18] Wigginton JE, et al. A note on exact tests of HWE. *Am J Hum Genet*. 2005;76(5):887–893.

[19] Weir BS, Cockerham CC. Estimating F-statistics. *Evolution*. 1984;38(6):1358–1370.

[20] Brookfield JFY. Null allele frequency. *Mol Ecol*. 1996;5(3):453–455.

[21] Excoffier L, et al. AMOVA. *Genetics*. 1992;131(2):479–491.

[22] Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*. 1995;57(1):289–300.

[23] Liu TD, et al. Identification of a novel QTL on LG16 associated with acute salt tolerance in red tilapia. *Mar Biotechnol*. 2025;27(1):39.
