# Changelog

All notable changes to ABI-SNP Analyzer will be documented in this file.

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
