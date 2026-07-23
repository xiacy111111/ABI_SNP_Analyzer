# Third-Party Notices

ABI-SNP Analyzer invokes the following separately distributed command-line
tools at runtime. These tools are **not** part of the ABI-SNP Analyzer source
code and are **not** covered by the MIT License. Users must download them
independently and comply with their respective licenses.

---

## 1. MUSCLE (Multiple Sequence Comparison by Log-Expectation)

| Field | Detail |
|---|---|
| **Purpose** | Multiple sequence alignment of IUPAC-coded FASTA sequences |
| **Version used in development** | v5.3 (5.3) |
| **Author** | Robert C. Edgar |
| **Website** | https://www.drive5.com/muscle/ |
| **Download URL** | https://www.drive5.com/muscle/downloads.htm |
| **License** | MUSCLE is free for non-commercial use. See https://www.drive5.com/muscle/license.htm |
| **Modifications** | None — used as-is, binary unmodified |
| **SHA256** | Not computed (user-provided binary) |
| **Citation** | Edgar, R.C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Research*, 32(5), 1792–1797. |

### How invoked
```
muscle -align input.fasta -output alignment.aln
```
For MUSCLE v3: `muscle -in input.fasta -out alignment.aln` (auto-detected at runtime).

---

## 2. NCBI BLAST+ (Basic Local Alignment Search Tool)

| Field | Detail |
|---|---|
| **Purpose** | Reference genome alignment of consensus sequences for SNP coordinate mapping |
| **Version used in development** | 2.17.0+ (ncbi-blast-2.17.0+-x64-win64) |
| **Publisher** | National Center for Biotechnology Information (NCBI) |
| **Website** | https://blast.ncbi.nlm.nih.gov/ |
| **Download URL** | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| **License** | Public Domain (U.S. Government work). See https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE |
| **Modifications** | None — used as-is, no source modifications |
| **SHA256** | Not computed (distributed as-is from NCBI) |
| **Citation** | Camacho, C. et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. |

### How invoked
```
blastn -query consensus.fasta -subject reference.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -out output.tsv -evalue 1e-10 -dust no -strand both
```
Version auto-detected at runtime via `blastn -version`.

---

## 3. Python Dependencies

The Python source code depends on the following third-party libraries, all
installable via `pip install -r requirements.txt`:

| Package | Typical Version | License | Purpose |
|---|---|---|---|
| numpy | ≥1.24 | BSD-3-Clause | Numerical arrays |
| scipy | ≥1.11 | BSD-3-Clause | Statistical tests (chi2, Fisher) |
| pandas | ≥2.0 | BSD-3-Clause | DataFrames, Excel I/O |
| matplotlib | ≥3.7 | PSF-based | Plotting |
| openpyxl | ≥3.1 | MIT | Excel read/write |
| biopython | ≥1.81 | Biopython License | FASTA/alignment I/O |

The exact versions pinned during development are recorded in `requirements.txt`.

---

## 4. Bundled Third-Party Libraries (in dist/)

The standalone Windows executable (`dist/ABI_SNP_Analyzer_v1.0.0.exe`) is built
with PyInstaller and bundles runtime components from:

- **Tcl/Tk** (BSD-style license) — GUI toolkit
- **Python standard library** (PSF License) — Python runtime

These are included solely for the purpose of running the frozen application and
are not modified. Their original license terms are retained in the distribution.

---

*Last updated: 2026-07-20*
