# External Tools

ABI-SNP Analyzer relies on two external command-line tools that are **not bundled**
in this repository. Please download them separately:

## MUSCLE v5 (Multiple Sequence Alignment)

- **Website:** <https://www.drive5.com/muscle/>
- **File:** `muscle.exe`
- **Windows download:** <https://www.drive5.com/muscle/downloads.htm>
- **Usage:** Place `muscle.exe` in the same directory as the ABI-SNP Analyzer
  executable, or add it to your system PATH.

## BLAST+ (Reference Genome Mapping, optional)

- **Website:** <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
- **File:** `blastn.exe` (part of the BLAST+ suite)
- **Windows download:** `ncbi-blast-*-x64-win64.tar.gz`
- **Usage:** Place the extracted `bin/` folder (or just `blastn.exe`) in the
  same directory as the ABI-SNP Analyzer executable, or add it to your system
  PATH.

## Reference Genome (optional)

A reference genome in FASTA format is required for SNP coordinate mapping.
For Nile tilapia (*Oreochromis niloticus*) analysis, download:

- **GenBank:** GCF_001858045.2
- **URL:** <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/858/045/GCF_001858045.2_O_niloticus_UMD_NMBU/>
