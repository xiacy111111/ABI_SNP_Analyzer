# Contributing to ABI-SNP Analyzer

Thank you for your interest in contributing!

## Report a Bug

1. Check existing [Issues](https://github.com/xiacy111111/ABI_SNP_Analyzer/issues) to avoid duplicates.
2. Open a new issue with:
   - A clear title and description
   - Steps to reproduce
   - Expected vs actual behavior
   - Screenshot of the error (if applicable)
   - Your Windows version and Python version (if running from source)

## Suggest a Feature

Open an issue with the "enhancement" label. Describe the problem this feature
solves and how it fits into the existing workflow.

## Development Setup

```bash
git clone https://github.com/xiacy111111/ABI_SNP_Analyzer.git
cd abi-snp-analyzer
python -m venv .venv && .venv\Scripts\activate
pip install -r requirements.txt
```

## Running Tests

```bash
pytest tests/ -v
```

## Code Style

- Follow PEP 8 conventions
- Use explicit, descriptive variable names
- Include docstrings for all public functions and classes
- Keep the Tkinter GUI code in `src/abi_to_snp_analyzer_v1.0.0.py`
- Place statistical logic in `population_genetics.py`
- Place chromatogram helpers in `chromatogram.py`

## Pull Request Process

1. Fork the repository and create a feature branch.
2. Run `pytest tests/ -v` to ensure all tests pass.
3. Update `CHANGELOG.md` under the "Unreleased" section.
4. Open a PR against the `main` branch with a clear description.
