# -*- coding: utf-8 -*-
"""SNP genotype format readers/writers for ABI-SNP Analyzer.

This module isolates every genotype-file import/export routine from the GUI so
that the round-trip ``export -> import -> analysis`` can be tested head-less.

Supported input formats (auto-detected, optionally forced by the caller)
-----------------------------------------------------------------------
=======================  ================================================
format id                description
=======================  ================================================
``genotype_matrix``      native ``genotype_matrix.xlsx`` (sheet
                         ``Genotype_Matrix``)
``genepop``              GENEPOP 4.x (``Pop`` blocks)
``arlequin``             Arlequin 3.5 project file (``.arp``)
``hierfstat``            hierfstat / R table (header row ``pop <loci...>``)
``genalex``              GenAlEx 6.5 (3 header rows, 2 columns per locus)
``fstat``                FSTAT 2.9.3 (``np nl nu digit_code`` header)
``popgene``              POPGENE32 diploid (``.pop``)
``structure``            STRUCTURE two-column-per-locus (``-9`` = missing)
=======================  ================================================

Allele coding convention
------------------------
Every numeric format written here uses one fixed, documented coding::

    1 = A, 2 = T, 3 = C, 4 = G, 0 = missing

Formats that allow free text (GENEPOP title line, Arlequin ``Title``, GenAlEx
dataset title) additionally carry an ``ALLELE_CODE:1=A,2=T,3=C,4=G`` tag so the
mapping is recovered verbatim on re-import instead of being assumed.

Export-only formats
-------------------
``write_vcf`` (v1.1.0) additionally emits VCF 4.2 with one sample column per
individual.  VCF files are produced for downstream tools (bcftools, vcftools,
tabix) and are *not* re-imported: there is no VCF reader and
``detect_format`` does not recognise ``.vcf``.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "MISSING_CALL",
    "GenotypeData",
    "SnpFormatError",
    "AmbiguousFormatError",
    "detect_format",
    "read_genotype_file",
    "write_hierfstat",
    "write_genalex",
    "write_genepop",
    "write_arlequin",
    "write_fstat",
    "write_popgene",
    "write_structure",
    "write_vcf",
    "vcf_locus_metadata_from_mapping",
    "FORMAT_LABELS",
]

MISSING_CALL = "--"

DEFAULT_CODE_TO_BASE = {1: "A", 2: "T", 3: "C", 4: "G"}
DEFAULT_BASE_TO_CODE = {base: code for code, base in DEFAULT_CODE_TO_BASE.items()}

#: Human readable names used by the GUI.
FORMAT_LABELS = {
    "genotype_matrix": "Genotype Matrix (Excel, native)",
    "genepop": "GENEPOP",
    "arlequin": "Arlequin (.arp)",
    "hierfstat": "hierfstat (R table)",
    "genalex": "GenAlEx (CSV)",
    "fstat": "FSTAT (.dat)",
    "popgene": "POPGENE32 (.pop)",
    "structure": "STRUCTURE",
}

#: Order in which auto-detection tries the parsers.  Cheap / unambiguous
#: signatures come first; ambiguous numeric headers are validated last.
_DETECTION_ORDER = (
    "genotype_matrix",
    "arlequin",
    "genepop",
    "hierfstat",
    "popgene",
    "genalex",
    "fstat",
    "structure",
)

_ALLELE_CODE_TAG_RE = re.compile(
    r"ALLELE_CODE\s*[:=]\s*([0-9A-Za-z=:,;\s]+?)(?:\||$)")
_ALLELE_CODE_ITEM_RE = re.compile(r"^\s*(\d+)\s*[:=]\s*([ACGTacgt])\s*$")
_NAMES_TAG_RE = re.compile(r"\b(LOCI|POPS)\s*[:=]\s*([^|]*)")

#: Maximum length of the machine-readable tag block written into a free-text
#: header.  Keeps the Arlequin ``Title``/GENEPOP title line reasonable even for
#: very large SNP panels.
_TAG_BLOCK_LIMIT = 2000


class SnpFormatError(ValueError):
    """Raised when a file cannot be parsed as the requested SNP format."""


class AmbiguousFormatError(SnpFormatError):
    """Raised when auto-detection cannot settle on a single format."""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _parse_allele_code_map(text: Optional[str]) -> Optional[Dict[int, str]]:
    """Recover an ``{allele_code: base}`` map embedded in a free-text header."""
    if not text:
        return None
    match = _ALLELE_CODE_TAG_RE.search(text)
    if not match:
        return None
    mapping: Dict[int, str] = {}
    for chunk in re.split(r"[,;]", match.group(1)):
        item = _ALLELE_CODE_ITEM_RE.match(chunk)
        if item:
            mapping[int(item.group(1))] = item.group(2).upper()
    return mapping or None


def _allele_code_tag(code_to_base: Optional[Dict[int, str]] = None) -> str:
    """Render the machine-readable allele-code tag for a header line."""
    code_to_base = code_to_base or DEFAULT_CODE_TO_BASE
    return "ALLELE_CODE:" + ",".join(
        f"{code}={code_to_base[code]}" for code in sorted(code_to_base))


def _names_tag(key: str, names: Sequence[str]) -> str:
    """Render a ``LOCI=...`` / ``POPS=...`` tag (``;`` separated)."""
    return f"{key}=" + ";".join(str(name) for name in names)


def _parse_names_tag(text: Optional[str], key: str) -> Optional[List[str]]:
    """Recover a ``;``-separated name list embedded in a free-text header."""
    if not text:
        return None
    for match in _NAMES_TAG_RE.finditer(text):
        if match.group(1) == key:
            names = [part.strip() for part in match.group(2).split(";")]
            names = [n for n in names if n]
            return names or None
    return None


def _read_text(path: str) -> str:
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            with open(path, "r", encoding=encoding) as handle:
                return handle.read()
        except UnicodeDecodeError:
            continue
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def _split(line: str, sep: Optional[str] = None) -> List[str]:
    """Split a data line, honouring an explicit or sniffed delimiter."""
    if sep:
        return [cell.strip() for cell in line.split(sep)]
    if "\t" in line:
        return [cell.strip() for cell in line.split("\t")]
    if "," in line and "=" not in line and '"' not in line:
        return [cell.strip() for cell in line.split(",")]
    return line.split()


def _dedup(names: Sequence[str]) -> List[str]:
    seen: Dict[str, int] = {}
    out: List[str] = []
    for name in names:
        base = str(name) if str(name).strip() else "Ind"
        count = seen.get(base, 0)
        seen[base] = count + 1
        out.append(base if count == 0 else f"{base}#{count + 1}")
    return out


def _make_call(allele1: Optional[str], allele2: Optional[str]) -> str:
    """Build the internal diploid call string (``"AT"``, ``"AA"``...)."""
    if not allele1 and not allele2:
        return MISSING_CALL
    if allele1 and allele2:
        return f"{allele1}{allele2}"
    single = allele1 or allele2
    return f"{single}{single}"


def _decode_allele(token: str, code_to_base: Dict[int, str],
                   missing_tokens: Sequence[str] = ()) -> Optional[str]:
    """Translate a single allele token into a base letter (or ``None``)."""
    token = str(token).strip()
    if not token or token in missing_tokens:
        return None
    if token.upper() in ("A", "C", "G", "T"):
        return token.upper()
    try:
        code = int(float(token))
    except (TypeError, ValueError):
        return None
    if code == 0:
        return None
    return code_to_base.get(code)


def _decode_genotype(token: str, width: int, code_to_base: Dict[int, str],
                     missing_tokens: Sequence[str] = ()) -> str:
    """Decode a concatenated diploid genotype such as ``"0102"`` -> ``"AT"``."""
    raw = str(token).strip()
    if not raw or raw in missing_tokens or set(raw) <= set(missing_tokens or ""):
        return MISSING_CALL
    if width <= 1:
        if len(raw) % 2:
            raise SnpFormatError(f"unexpected genotype token {raw!r}")
        width = len(raw) // 2
    if len(raw) != 2 * width:
        raise SnpFormatError(
            f"genotype token {raw!r} is not {2 * width} characters wide")
    allele1 = _decode_allele(raw[:width], code_to_base, missing_tokens)
    allele2 = _decode_allele(raw[width:], code_to_base, missing_tokens)
    return _make_call(allele1, allele2)


def _infer_width(tokens: Sequence[str], expected_loci: Optional[int] = None) -> int:
    """Infer how many characters encode one allele (1, 2 or 3)."""
    usable = [t for t in tokens if t and set(t) != {"0"}]
    if not usable:
        return 1
    probe = usable[0]
    widths = {len(probe) // 2}
    if expected_loci and len(tokens) == expected_loci:
        candidate = len(probe) // 2
        if candidate >= 1:
            widths.add(candidate)
    width = max(widths)
    return min(max(width, 1), 3)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

class GenotypeData(object):
    """In-memory genotype table recovered from (or bound for) a file."""

    def __init__(self, samples: List[str], loci: List[str],
                 calls: Dict[str, Dict[str, str]],
                 populations: Optional[Dict[str, List[str]]] = None,
                 warnings: Optional[List[str]] = None,
                 source_format: str = ""):
        self.samples = list(samples)
        self.loci = list(loci)
        self.calls = calls                       # locus -> sample -> call
        self.populations = populations or {}     # population -> [samples]
        self.warnings = list(warnings or [])
        self.source_format = source_format

    # -- convenience ---------------------------------------------------
    @property
    def n_samples(self) -> int:
        return len(self.samples)

    @property
    def n_loci(self) -> int:
        return len(self.loci)

    def to_matrix_dict(self) -> Dict[str, Dict[str, str]]:
        """Return ``{locus: {sample: call}}`` suitable for the GUI."""
        matrix: Dict[str, Dict[str, str]] = {}
        for locus in self.loci:
            matrix[locus] = {
                sample: self.calls.get(locus, {}).get(sample, MISSING_CALL)
                for sample in self.samples
            }
        return matrix

    def summary(self) -> str:
        return (f"{self.n_samples} samples x {self.n_loci} loci"
                + (f", {len(self.populations)} populations"
                   if self.populations else ", no population structure"))


def _finalise(samples: List[str], loci: List[str],
              calls: Dict[str, Dict[str, str]],
              populations: Dict[str, List[str]],
              warnings: List[str], source_format: str) -> GenotypeData:
    if not samples:
        raise SnpFormatError("no samples found in file")
    if not loci:
        raise SnpFormatError("no loci found in file")
    known = set(samples)
    for population, members in list(populations.items()):
        kept = [m for m in members if m in known]
        if kept:
            populations[population] = kept
        else:
            del populations[population]
    return GenotypeData(samples, loci, calls, populations, warnings, source_format)


# ===========================================================================
# Readers
# ===========================================================================

def read_genotype_matrix_xlsx(path: str) -> GenotypeData:
    """Read the native ``genotype_matrix.xlsx`` produced by the analyser."""
    import pandas as pd  # imported lazily; keeps the module GUI-agnostic

    warnings: List[str] = []
    try:
        frame = pd.read_excel(path, sheet_name="Genotype_Matrix")
    except ValueError:
        frame = pd.read_excel(path, sheet_name=0)
        warnings.append("Sheet 'Genotype_Matrix' not found — used the first sheet.")

    if frame.shape[1] < 2:
        raise SnpFormatError("genotype matrix needs at least one locus column")

    samples = _dedup([str(v) for v in frame.iloc[:, 0].tolist()])
    loci = [str(c) for c in frame.columns[1:]]
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    for locus in loci:
        column = frame[locus]
        for sample, value in zip(samples, column.tolist()):
            text = "" if value is None else str(value).strip()
            if text.lower() in ("", "nan", "none", "na", "n/a", "--", "nn"):
                calls[locus][sample] = MISSING_CALL
            else:
                calls[locus][sample] = text.upper()
    return _finalise(samples, loci, calls, {}, warnings, "genotype_matrix")


def read_genepop(path: str) -> GenotypeData:
    """Read a GENEPOP 4.x file (locus list, then ``Pop`` blocks)."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines()]
    if not lines:
        raise SnpFormatError("empty GENEPOP file")

    code_to_base = _parse_allele_code_map(lines[0]) or DEFAULT_CODE_TO_BASE
    warnings: List[str] = []
    pop_hint = _parse_names_tag(lines[0], "POPS")

    idx = 1                                   # line 0 is the free-text title
    loci: List[str] = []
    while idx < len(lines):
        stripped = lines[idx].strip()
        if stripped.lower() == "pop":
            break
        if stripped:
            loci.append(stripped)
        idx += 1
    if not loci:
        raise SnpFormatError("no locus names found before the first 'Pop' line")
    if idx >= len(lines):
        raise SnpFormatError("no 'Pop' section found")

    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}
    pop_counter = 0
    width: Optional[int] = None

    for line in lines[idx:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.lower().startswith("pop"):
            pop_counter += 1
            if pop_hint and pop_counter <= len(pop_hint):
                pop_name = pop_hint[pop_counter - 1]
            else:
                pop_name = f"Pop{pop_counter}"
            populations.setdefault(pop_name, [])
            continue
        if not populations:
            raise SnpFormatError("genotype data found before any 'Pop' line")
        pop_name = list(populations.keys())[-1]
        if "," in stripped:
            sample_id, payload = stripped.split(",", 1)
        else:
            sample_id, payload = "", stripped
        sample_id = sample_id.strip() or f"Ind_{len(samples) + 1:04d}"
        tokens = payload.split()
        if len(tokens) != len(loci):
            raise SnpFormatError(
                f"sample {sample_id!r} has {len(tokens)} genotypes, "
                f"expected {len(loci)}")
        if width is None:
            width = _infer_width(tokens, len(loci))
        for locus, token in zip(loci, tokens):
            try:
                calls[locus][sample_id] = _decode_genotype(
                    token, width, code_to_base, missing_tokens=("?", "0"))
            except SnpFormatError as exc:
                raise SnpFormatError(f"{exc} (sample {sample_id!r}, locus {locus!r})")
        samples.append(sample_id)
        populations[pop_name].append(sample_id)

    return _finalise(samples, loci, calls, populations, warnings, "genepop")


def read_arlequin(path: str) -> GenotypeData:
    """Read an Arlequin 3.5 ``.arp`` STANDARD genotypic data file."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines()]
    warnings: List[str] = []

    profile: Dict[str, str] = {}
    in_profile = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            in_profile = stripped.lower().startswith("[profile]")
            if stripped.lower().startswith("[data]"):
                break
            continue
        if in_profile and "=" in stripped:
            key, value = stripped.split("=", 1)
            profile[key.strip().lower()] = value.strip().strip('"')

    title = profile.get("title", "")
    code_to_base = _parse_allele_code_map(title) or DEFAULT_CODE_TO_BASE
    loci_hint = _parse_names_tag(title, "LOCI")
    missing_char = profile.get("missingdata", "?").strip('"')[:1] or "?"
    separator = profile.get("locusseparator", "NONE").upper()
    if separator not in ("NONE", "WHITESPACE", "TAB"):
        separator = "NONE"

    # ---- walk the [[Samples]] blocks -------------------------------------
    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {}
    populations: Dict[str, List[str]] = {}
    loci: List[str] = []

    pop_name = ""
    block: List[str] = []
    in_block = False
    declared_size: Optional[int] = None

    def flush_block(name: str, rows: List[str], size: Optional[int]):
        if not rows:
            return
        drop_last = bool(rows) and rows[-1].strip().startswith("}")
        body = rows[:-1] if drop_last else rows
        body = [r for r in body if r.strip()]
        if not body:
            return
        one_line_mode = size is not None and len(body) == size
        members: List[str] = []
        index = 0
        while index < len(body):
            if one_line_mode:
                tokens = body[index].split()
                sample_id = tokens[0] if tokens else f"Ind_{len(samples) + 1:04d}"
                payload = tokens[2:] if len(tokens) >= 3 else tokens[1:]
                index += 1
            else:
                first = body[index].split()
                sample_id = first[0] if first else f"Ind_{len(samples) + 1:04d}"
                payload = first[2:] if len(first) >= 3 else first[1:]
                index += 1
                if index < len(body):
                    payload2 = body[index].split()
                    index += 1
                else:
                    payload2 = []
            if separator == "NONE":
                strand1 = "".join(payload)
                strand2 = "".join(payload2)
                n_loci = max(len(strand1), len(strand2))
                tokens1 = list(strand1) if strand1 else []
                tokens2 = list(strand2) if strand2 else []
            else:
                tokens1 = list(payload)
                tokens2 = list(payload2)
                n_loci = len(tokens1)
            if one_line_mode:
                merged = tokens1 + tokens2
                if len(merged) % 2:
                    raise SnpFormatError(
                        f"odd number of alleles for sample {sample_id!r}")
                tokens1 = merged[0::2]
                tokens2 = merged[1::2]
                n_loci = len(tokens1)
            if not loci:
                if loci_hint and len(loci_hint) == n_loci:
                    loci.extend(loci_hint)
                else:
                    for position in range(n_loci):
                        loci.append(f"Locus_{position + 1}")
            elif n_loci != len(loci):
                raise SnpFormatError(
                    f"sample {sample_id!r} has {n_loci} loci, expected {len(loci)}")
            for position, (tok1, tok2) in enumerate(zip(tokens1, tokens2)):
                locus = loci[position]
                calls.setdefault(locus, {})
                allele1 = _decode_allele(tok1, code_to_base, (missing_char, "0"))
                allele2 = _decode_allele(tok2, code_to_base, (missing_char, "0"))
                calls[locus][sample_id] = _make_call(allele1, allele2)
            samples.append(sample_id)
            members.append(sample_id)
        if members:
            populations.setdefault(name, []).extend(members)

    for line in lines:
        stripped = line.strip()
        if stripped.lower().startswith("samplename="):
            if in_block:
                flush_block(pop_name, block, declared_size)
            pop_name = stripped.split("=", 1)[1].strip().strip('"') or \
                f"Pop{len(populations) + 1}"
            block = []
            declared_size = None
            continue
        if stripped.lower().startswith("samplesize="):
            try:
                declared_size = int(stripped.split("=", 1)[1].strip())
            except ValueError:
                declared_size = None
            continue
        if "sampledata" in stripped.lower() and "{" in stripped:
            in_block = True
            block = []
            continue
        if in_block:
            if stripped.startswith("}"):
                flush_block(pop_name, block, declared_size)
                in_block = False
                block = []
                pop_name = ""
                continue
            block.append(line)
    if in_block:
        flush_block(pop_name, block, declared_size)

    if loci and all(locus.startswith("Locus_") for locus in loci):
        warnings.append(
            "Arlequin (.arp) does not store locus names — loci were renamed "
            "Locus_1..Locus_N. Export GENEPOP/hierfstat/GenAlEx to keep names.")
    return _finalise(samples, loci, calls, populations, warnings, "arlequin")


def read_hierfstat(path: str) -> GenotypeData:
    """Read a hierfstat / R table: header ``pop <locus...>`` + one row/individual."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise SnpFormatError("hierfstat file needs a header row and data rows")

    sep = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    header = _split(lines[0], sep)
    if not header or header[0].strip().lower() not in ("pop", "population", "popid"):
        raise SnpFormatError("hierfstat header must start with a 'pop' column")
    loci = [str(name).strip() for name in header[1:]]
    if not loci:
        raise SnpFormatError("no locus columns in hierfstat header")

    warnings: List[str] = []
    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}
    width: Optional[int] = None

    for line in lines[1:]:
        cells = _split(line, sep)
        if len(cells) < len(loci) + 1:
            raise SnpFormatError(
                f"row has {len(cells)} fields, expected {len(loci) + 1}")
        pop_label = cells[0].strip() or "1"
        tokens = cells[1:1 + len(loci)]
        if width is None:
            width = _infer_width(tokens, len(loci))
        sample_id = f"{pop_label}_Ind{len(samples) + 1:04d}"
        for locus, token in zip(loci, tokens):
            calls[locus][sample_id] = _decode_genotype(
                token, width, DEFAULT_CODE_TO_BASE, missing_tokens=("?", "NA"))
        samples.append(sample_id)
        populations.setdefault(f"Pop{pop_label}", []).append(sample_id)

    warnings.append(
        "hierfstat files carry no individual IDs — sample names were "
        "regenerated as <pop>_Ind####.")
    return _finalise(samples, loci, calls, populations, warnings, "hierfstat")


def read_genalex(path: str) -> GenotypeData:
    """Read a GenAlEx 6.5 file (3 header rows, two columns per diploid locus)."""
    text = _read_text(path)
    delimiter = _sniff_delimiter(text)
    table = _csv_rows(text, delimiter)
    if len(table) < 4:
        raise SnpFormatError("GenAlEx file needs 3 header rows plus data rows")

    head = table[0]
    if len(head) < 3 or not all(_is_int(cell) for cell in head[:3]):
        raise SnpFormatError("GenAlEx header row must start with 3 integers")
    n_loci, n_samples_declared, n_pops = (int(head[0]), int(head[1]), int(head[2]))
    if n_loci < 1 or n_samples_declared < 1:
        raise SnpFormatError("GenAlEx header declares no data")

    title = table[1][0] if len(table) > 1 and table[1] else ""
    code_to_base = _parse_allele_code_map(title) or DEFAULT_CODE_TO_BASE

    locus_row = table[2] if len(table) > 2 else []
    loci: List[str] = []
    position = 2
    while len(loci) < n_loci and position < len(locus_row):
        name = locus_row[position].strip()
        loci.append(name if name else f"Locus_{len(loci) + 1}")
        position += 2
    if len(loci) != n_loci:
        raise SnpFormatError("GenAlEx locus-name row does not match the locus count")

    warnings: List[str] = []
    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}

    data_rows = [row for row in table[3:] if any(cell.strip() for cell in row)]
    data_rows = data_rows[:n_samples_declared]
    if len(data_rows) != n_samples_declared:
        raise SnpFormatError("GenAlEx data rows do not match the declared sample count")

    for cells in data_rows:
        if len(cells) < 2 + 2 * n_loci:
            raise SnpFormatError(
                f"GenAlEx row has {len(cells)} fields, expected {2 + 2 * n_loci}")
        sample_id = cells[0].strip() or f"Ind_{len(samples) + 1:04d}"
        pop_name = cells[1].strip() or f"Pop{len(populations) + 1}"
        for locus_index, locus in enumerate(loci):
            allele1 = _decode_allele(cells[2 + 2 * locus_index], code_to_base, ("?", "0"))
            allele2 = _decode_allele(cells[3 + 2 * locus_index], code_to_base, ("?", "0"))
            calls[locus][sample_id] = _make_call(allele1, allele2)
        samples.append(sample_id)
        populations.setdefault(pop_name, []).append(sample_id)

    if len(populations) != n_pops:
        warnings.append(
            f"GenAlEx header declares {n_pops} populations, "
            f"{len(populations)} were recovered.")
    return _finalise(samples, loci, calls, populations, warnings, "genalex")


def read_fstat(path: str) -> GenotypeData:
    """Read an FSTAT 2.9.x ``.dat`` file."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise SnpFormatError("empty FSTAT file")
    head = lines[0].split()
    if len(head) < 4 or not all(_is_int(cell) for cell in head[:4]):
        raise SnpFormatError("FSTAT header must be 'np nl nu digit_code'")
    n_pops, n_loci, _max_allele, digit = (int(head[0]), int(head[1]),
                                          int(head[2]), int(head[3]))
    if n_loci < 1 or digit < 1:
        raise SnpFormatError("invalid FSTAT header values")
    loci = [ln.strip() for ln in lines[1:1 + n_loci]]
    if len(loci) != n_loci:
        raise SnpFormatError("FSTAT locus-name block is truncated")

    warnings: List[str] = []
    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}

    for line in lines[1 + n_loci:]:
        cells = line.split()
        if len(cells) < n_loci + 1:
            raise SnpFormatError(
                f"FSTAT row has {len(cells)} fields, expected {n_loci + 1}")
        pop_label = cells[0]
        sample_id = f"Pop{pop_label}_Ind{len(samples) + 1:04d}"
        for locus, token in zip(loci, cells[1:1 + n_loci]):
            calls[locus][sample_id] = _decode_genotype(
                token, digit, DEFAULT_CODE_TO_BASE, missing_tokens=("?",))
        samples.append(sample_id)
        populations.setdefault(f"Pop{pop_label}", []).append(sample_id)

    warnings.append(
        "FSTAT files carry no individual IDs — sample names were regenerated.")
    return _finalise(samples, loci, calls, populations, warnings, "fstat")


def read_popgene(path: str) -> GenotypeData:
    """Read a POPGENE32 diploid ``.pop`` file."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines()]
    warnings: List[str] = []

    n_loci = None
    for line in lines:
        match = re.search(r"Number\s+of\s+loci\s*=\s*(\d+)", line, re.I)
        if match:
            n_loci = int(match.group(1))
            break
    if not n_loci:
        raise SnpFormatError("POPGENE header 'Number of loci =' not found")

    # locus names sit on the line right after 'Locus name :'
    loci: List[str] = []
    for index, line in enumerate(lines):
        if re.search(r"Locus\s*name\s*:", line, re.I):
            if index + 1 < len(lines):
                loci = [token for token in lines[index + 1].split()
                        if token.lower() != "and"]
            break
    if len(loci) != n_loci:
        raise SnpFormatError("POPGENE locus names do not match the locus count")

    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}

    started = False
    pop_index = 0
    current: List[str] = []
    pending_blank = False
    for line in lines:
        stripped = line.strip()
        if not started:
            if stripped.lower().startswith("locus") or re.match(
                    r"^\s*[ACGTacgt.]{2}(\s+[ACGTacgt.]{2})*\s*$", stripped):
                started = bool(re.match(
                    r"^\s*[ACGTacgt.]{2}(\s+[ACGTacgt.]{2})*\s*$", stripped))
                if not started:
                    continue
            else:
                continue
        if not stripped:
            pending_blank = True
            continue
        if pending_blank and current:
            pop_index += 1
            populations[f"Pop{pop_index}"] = current
            current = []
            pending_blank = False
        tokens = stripped.split()
        if len(tokens) != n_loci:
            raise SnpFormatError(
                f"POPGENE genotype row has {len(tokens)} calls, expected {n_loci}")
        sample_id = f"Pop{pop_index + 1}_Ind{len(current) + 1:04d}"
        for locus, token in zip(loci, tokens):
            text_value = token.upper()
            if set(text_value) <= {"."}:
                calls[locus][sample_id] = MISSING_CALL
            else:
                calls[locus][sample_id] = text_value
        samples.append(sample_id)
        current.append(sample_id)
    if current:
        populations[f"Pop{pop_index + 1}"] = current

    warnings.append(
        "POPGENE files carry no individual IDs or population names — both "
        "were regenerated.")
    return _finalise(samples, loci, calls, populations, warnings, "popgene")


def read_structure(path: str) -> GenotypeData:
    """Read the STRUCTURE export (id, pop, then two allele columns per locus)."""
    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise SnpFormatError("empty STRUCTURE file")
    sep = "\t" if "\t" in lines[0] else None
    warnings: List[str] = []

    first = _split(lines[0], sep)
    if len(first) < 4 or (len(first) - 2) % 2:
        raise SnpFormatError("STRUCTURE rows need id, pop and 2 columns per locus")
    n_loci = (len(first) - 2) // 2
    loci = [f"Locus_{i + 1}" for i in range(n_loci)]

    samples: List[str] = []
    calls: Dict[str, Dict[str, str]] = {locus: {} for locus in loci}
    populations: Dict[str, List[str]] = {}

    for line in lines:
        cells = _split(line, sep)
        if len(cells) < 2 + 2 * n_loci:
            raise SnpFormatError("STRUCTURE row length mismatch")
        sample_id = cells[0].strip() or f"Ind_{len(samples) + 1:04d}"
        pop_label = cells[1].strip() or "1"
        for index, locus in enumerate(loci):
            allele1 = _decode_allele(cells[2 + 2 * index], DEFAULT_CODE_TO_BASE,
                                     ("-9", "?"))
            allele2 = _decode_allele(cells[3 + 2 * index], DEFAULT_CODE_TO_BASE,
                                     ("-9", "?"))
            calls[locus][sample_id] = _make_call(allele1, allele2)
        samples.append(sample_id)
        populations.setdefault(f"Pop{pop_label}", []).append(sample_id)

    warnings.append(
        "STRUCTURE files carry no locus names — loci were renamed Locus_1..N.")
    return _finalise(samples, loci, calls, populations, warnings, "structure")


_READERS = {
    "genotype_matrix": read_genotype_matrix_xlsx,
    "genepop": read_genepop,
    "arlequin": read_arlequin,
    "hierfstat": read_hierfstat,
    "genalex": read_genalex,
    "fstat": read_fstat,
    "popgene": read_popgene,
    "structure": read_structure,
}


def _sniff_delimiter(text: str) -> str:
    """Return the field delimiter used by a delimited text file."""
    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")
    if "\t" in first_line:
        return "\t"
    if ";" in first_line and "," not in first_line:
        return ";"
    return ","


def _csv_rows(text: str, delimiter: str) -> List[List[str]]:
    """Parse a delimited text block into rows, honouring quoted fields."""
    import csv as _csv
    import io

    reader = _csv.reader(io.StringIO(text), delimiter=delimiter,
                         quotechar='"', skipinitialspace=True)
    return [[cell.strip() for cell in row] for row in reader if row]


def _is_int(text: str) -> bool:
    try:
        int(str(text).strip())
        return True
    except (TypeError, ValueError):
        return False


# ===========================================================================
# Format detection
# ===========================================================================

def detect_format(path: str) -> str:
    """Guess the SNP genotype format of ``path``.

    Raises ``AmbiguousFormatError`` when no parser accepts the file.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension in (".xlsx", ".xls"):
        return "genotype_matrix"

    text = _read_text(path)
    lines = [ln.rstrip("\r\n") for ln in text.splitlines()]
    non_empty = [ln for ln in lines if ln.strip()]
    if not non_empty:
        raise AmbiguousFormatError("file is empty")

    candidates: List[str] = []
    if "[Profile]" in text or "[Data]" in text:
        candidates.append("arlequin")
    if any(ln.strip().lower().startswith("pop") and
           len(ln.strip()) <= 12 for ln in non_empty):
        candidates.append("genepop")
    header = _split(non_empty[0])
    if header and header[0].strip().lower() in ("pop", "population", "popid"):
        candidates.append("hierfstat")
    if re.search(r"Number\s+of\s+(populations|loci)\s*=", text, re.I):
        candidates.append("popgene")
    if header and len(header) >= 3 and all(_is_int(cell) for cell in header[:3]):
        candidates.extend(["genalex", "fstat"])
    candidates.append("structure")

    ordered = [fmt for fmt in _DETECTION_ORDER if fmt in candidates]
    ordered += [fmt for fmt in _DETECTION_ORDER if fmt not in ordered]

    failures: List[str] = []
    for fmt in ordered:
        try:
            data = _READERS[fmt](path)
        except Exception as exc:                     # noqa: BLE001 - probe only
            failures.append(f"{fmt}: {exc}")
            continue
        if data.n_samples and data.n_loci:
            return fmt
    raise AmbiguousFormatError(
        "could not identify the SNP genotype format. Tried: "
        + "; ".join(failures[:4]))


def read_genotype_file(path: str, fmt: Optional[str] = None) -> GenotypeData:
    """Read ``path``; ``fmt`` forces a parser (see :data:`FORMAT_LABELS`)."""
    if fmt:
        key = fmt.lower()
        if key not in _READERS:
            raise SnpFormatError(f"unsupported SNP genotype format: {fmt!r}")
        return _READERS[key](path)
    return _READERS[detect_format(path)](path)


# ===========================================================================
# Writers
# ===========================================================================

def _ordered_populations(samples: Sequence[str], loci: Sequence[str],
                         calls: Dict[str, Dict[str, str]]) -> Dict[str, List[str]]:
    return {"Pop1": list(samples)}


def _resolve_populations(samples: Sequence[str],
                         populations: Optional[Dict[str, List[str]]]
                         ) -> Tuple[Dict[str, List[str]], List[str]]:
    """Return ``({pop: [samples]}, ordered_samples)`` with every sample placed."""
    groups: Dict[str, List[str]] = {}
    if populations:
        for name, members in populations.items():
            members = [m for m in members if m in set(samples)]
            if members:
                groups[str(name)] = sorted(members)
    assigned = {m for members in groups.values() for m in members}
    leftover = [s for s in samples if s not in assigned]
    if leftover:
        if groups:                      # keep partial assignments visible
            groups.setdefault("Unassigned", []).extend(leftover)
        else:
            groups = {"Pop1": list(leftover)}
    if not groups:
        groups = {"Pop1": list(samples)}
    ordered: List[str] = []
    for name in sorted(groups):
        ordered.extend(groups[name])
    return groups, ordered


def _allele_code(base: Optional[str]) -> int:
    if not base:
        return 0
    return DEFAULT_BASE_TO_CODE.get(str(base).upper(), 0)


def _call_alleles(call: str) -> Tuple[Optional[str], Optional[str]]:
    text = (call or "").strip().upper()
    if text in ("", MISSING_CALL, "NN", "N/A", "NA", "NAN", "0"):
        return (None, None)
    alleles = [c for c in text if c in "ACGT"]
    if not alleles:
        return (None, None)
    if len(alleles) == 1:
        return (alleles[0], alleles[0])
    return (alleles[0], alleles[1])


def write_hierfstat(path: str, samples: Sequence[str], loci: Sequence[str],
                    calls: Dict[str, Dict[str, str]],
                    populations: Optional[Dict[str, List[str]]] = None,
                    title: Optional[str] = None,
                    digits_per_allele: int = 1) -> int:
    """Write a hierfstat / R genotype table.

    Layout::

        pop <TAB> locus1 <TAB> locus2 ...
        1    <TAB> 11     <TAB> 12    ...

    Alleles are concatenated with ``digits_per_allele`` characters each
    (A=1, T=2, C=3, G=4); a fully missing call is written as ``0`` repeated
    ``2 * digits_per_allele`` times.  The table is read in R with
    ``read.table(file, header = TRUE)`` and consumed directly by
    ``hierfstat::basic.stats()`` / ``hierfstat::wc()``.
    """
    groups, ordered = _resolve_populations(samples, populations)
    pop_index = {name: i + 1 for i, name in enumerate(sorted(groups))}
    sample_pop = {}
    for name, members in groups.items():
        for member in members:
            sample_pop[member] = pop_index[name]

    lines = ["\t".join(["pop"] + [str(locus) for locus in loci])]
    for sample in ordered:
        row = [str(sample_pop.get(sample, 1))]
        for locus in loci:
            allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
            code1, code2 = _allele_code(allele1), _allele_code(allele2)
            row.append(f"{code1:0{digits_per_allele}d}{code2:0{digits_per_allele}d}")
        lines.append("\t".join(row))

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(ordered)


def write_genalex(path: str, samples: Sequence[str], loci: Sequence[str],
                  calls: Dict[str, Dict[str, str]],
                  populations: Optional[Dict[str, List[str]]] = None,
                  title: Optional[str] = None) -> int:
    """Write a GenAlEx 6.5 diploid codominant file (CSV).

    Layout (Peakall & Smouse 2012)::

        row 1 : nLoci, nInd, nPop, size_pop1, size_pop2, ...
        row 2 : dataset title, <blank>, <blank>, pop1, pop2, ...
        row 3 : "Sample", "Pop", locus1, <blank>, locus2, <blank>, ...
        row 4+: sample id, population, allele1, allele2, allele1, allele2, ...

    Alleles use the exported coding A=1, T=2, C=3, G=4; missing = 0.
    """
    groups, ordered = _resolve_populations(samples, populations)
    pop_names = sorted(groups)
    n_loci, n_ind, n_pops = len(loci), len(ordered), len(pop_names)

    dataset_title = title or "ABI-SNP Analyzer SNP genotypes"
    dataset_title = f"{dataset_title} | {_allele_code_tag()}"

    rows: List[List[str]] = []
    rows.append([str(n_loci), str(n_ind), str(n_pops)]
                + [str(len(groups[name])) for name in pop_names])
    rows.append([dataset_title, "", ""] + list(pop_names))
    locus_header = ["Sample", "Pop"]
    for locus in loci:
        locus_header.extend([str(locus), ""])
    rows.append(locus_header)

    sample_pop = {}
    for name, members in groups.items():
        for member in members:
            sample_pop[member] = name

    for sample in ordered:
        row = [str(sample), str(sample_pop.get(sample, ""))]
        for locus in loci:
            allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
            row.append(str(_allele_code(allele1)))
            row.append(str(_allele_code(allele2)))
        rows.append(row)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        for row in rows:
            handle.write(",".join(_csv_cell(cell) for cell in row) + "\n")
    return n_ind


def _csv_cell(value: str) -> str:
    text = "" if value is None else str(value)
    if any(ch in text for ch in (",", '"', "\n")):
        return '"' + text.replace('"', '""') + '"'
    return text


def write_genepop(path: str, samples: Sequence[str], loci: Sequence[str],
                  calls: Dict[str, Dict[str, str]],
                  populations: Optional[Dict[str, List[str]]] = None,
                  title: Optional[str] = None) -> int:
    """Write a GENEPOP 4.x file (2-digit allele codes, ``0000`` = missing)."""
    groups, ordered = _resolve_populations(samples, populations)
    title_text = title or "SNP Genotype Data (GENEPOP)"
    tags = [_allele_code_tag()]
    if len(_names_tag("POPS", sorted(groups))) <= _TAG_BLOCK_LIMIT:
        tags.append(_names_tag("POPS", sorted(groups)))
    header = f"{title_text} | " + " | ".join(tags)

    # GENEPOP reads one locus name per line, so a name only needs to be a
    # single whitespace-free token — sanitising is kept minimal (blank -> "_")
    # so that locus names survive the export -> import round trip intact.
    clean_loci = [_clean_genepop_locus(locus) for locus in loci]
    lines = [header] + clean_loci
    for name in sorted(groups):
        lines.append("Pop")
        for sample in groups[name]:
            tokens = []
            for locus in loci:
                allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
                code1, code2 = _allele_code(allele1), _allele_code(allele2)
                tokens.append(f"{code1:02d}{code2:02d}")
            safe_id = str(sample).replace(",", "_").replace(" ", "_")
            lines.append(f"{safe_id}, " + " ".join(tokens))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(ordered)


def _clean_genepop_locus(name: str) -> str:
    """Make a locus name a single whitespace-free token for the GENEPOP parser.

    Only blanks are replaced: GENEPOP reads one locus name per line, so every
    other printable character is safe and is preserved verbatim (this keeps
    names such as ``SNP01_A>C`` recoverable on re-import).
    """
    cleaned = re.sub(r"\s+", "_", str(name).strip())
    if not cleaned:
        cleaned = "Locus"
    if not cleaned[0].isalpha():
        cleaned = "Locus_" + cleaned
    return cleaned


def write_arlequin(path: str, samples: Sequence[str], loci: Sequence[str],
                   calls: Dict[str, Dict[str, str]],
                   populations: Optional[Dict[str, List[str]]] = None,
                   title: Optional[str] = None) -> int:
    """Write an Arlequin 3.5 STANDARD ``.arp`` file (two lines per individual)."""
    groups, ordered = _resolve_populations(samples, populations)
    tab = "\t"
    title_text = title or "ABI-SNP Analyzer SNP genotypes"
    tags = [_allele_code_tag()]
    if len(_names_tag("LOCI", loci)) <= _TAG_BLOCK_LIMIT:
        tags.append(_names_tag("LOCI", loci))   # .arp has no locus-name field
    dataset_title = f"{title_text} | " + " | ".join(tags)

    lines = [
        "[Profile]",
        f'{tab}Title="{dataset_title}"',
        f"{tab}NbSamples={len(groups)}",
        f"{tab}GenotypicData=1",
        f"{tab}GameticPhase=0",
        f"{tab}RecessiveData=0",
        f"{tab}DataType=STANDARD",
        f"{tab}LocusSeparator=NONE",
        f'{tab}MissingData="?"',
        "",
        "[Data]",
        f"{tab}[[Samples]]",
    ]
    for name in sorted(groups):
        lines.append(f'{tab}{tab}SampleName="{name}"')
        lines.append(f"{tab}{tab}SampleSize={len(groups[name])}")
        lines.append(f"{tab}{tab}SampleData={{")
        for sample in groups[name]:
            strand1, strand2 = [], []
            for locus in loci:
                allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
                if allele1 is None or allele2 is None:
                    strand1.append("?")
                    strand2.append("?")
                else:
                    strand1.append(str(_allele_code(allele1)))
                    strand2.append(str(_allele_code(allele2)))
            safe_id = str(sample).replace(".", "_")
            lines.append(f"{tab}{tab}{tab}{safe_id}{tab}1{tab}{''.join(strand1)}")
            lines.append(f"{tab}{tab}{tab}{tab}{''.join(strand2)}")
        lines.append(f"{tab}{tab}}}")
    lines += [
        "",
        "[Structure]",
        f'{tab}StructureName="All populations"',
        f"{tab}NbGroups=1",
        f"{tab}Group={{",
    ]
    for name in sorted(groups):
        lines.append(f'{tab}{tab}"{name}"')
    lines.append(f"{tab}}}")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(ordered)


def write_fstat(path: str, samples: Sequence[str], loci: Sequence[str],
                calls: Dict[str, Dict[str, str]],
                populations: Optional[Dict[str, List[str]]] = None,
                title: Optional[str] = None,
                digits_per_allele: int = 1) -> int:
    """Write an FSTAT 2.9.x ``.dat`` file (``0`` repeated = missing)."""
    groups, ordered = _resolve_populations(samples, populations)
    pop_index = {name: i + 1 for i, name in enumerate(sorted(groups))}
    sample_pop = {}
    for name, members in groups.items():
        for member in members:
            sample_pop[member] = pop_index[name]

    lines = [f"{len(groups)} {len(loci)} 4 {digits_per_allele}"]
    lines += [str(locus) for locus in loci]
    missing = "0" * (2 * digits_per_allele)
    for sample in ordered:
        row = [str(sample_pop.get(sample, 1))]
        for locus in loci:
            allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
            code1, code2 = _allele_code(allele1), _allele_code(allele2)
            row.append(missing if not code1 and not code2
                       else f"{code1:0{digits_per_allele}d}{code2:0{digits_per_allele}d}")
        lines.append(" ".join(row))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(ordered)


def write_popgene(path: str, samples: Sequence[str], loci: Sequence[str],
                  calls: Dict[str, Dict[str, str]],
                  populations: Optional[Dict[str, List[str]]] = None,
                  title: Optional[str] = None) -> int:
    """Write a POPGENE32 diploid ``.pop`` file (``..`` = missing)."""
    groups, ordered = _resolve_populations(samples, populations)
    lines = [
        f"/* {title or 'SNP Genotype Data'} */",
        f"Number of populations = {len(groups)}",
        f"Number of loci = {len(loci)}",
        "Locus name :",
        " " + " ".join(str(locus) for locus in loci),
        "",
    ]
    n_samples = 0
    for index, name in enumerate(sorted(groups)):
        if index:
            lines.append("")
        for sample in groups[name]:
            tokens = []
            for locus in loci:
                allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
                if allele1 is None or allele2 is None:
                    tokens.append("..")
                else:
                    tokens.append("".join(sorted([allele1.upper(), allele2.upper()])))
            lines.append(" " + " ".join(tokens))
            n_samples += 1
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return n_samples


def write_structure(path: str, samples: Sequence[str], loci: Sequence[str],
                    calls: Dict[str, Dict[str, str]],
                    populations: Optional[Dict[str, List[str]]] = None,
                    title: Optional[str] = None) -> int:
    """Write the STRUCTURE two-column-per-locus file (``-9`` = missing)."""
    groups, ordered = _resolve_populations(samples, populations)
    pop_index = {name: i + 1 for i, name in enumerate(sorted(groups))}
    sample_pop = {}
    for name, members in groups.items():
        for member in members:
            sample_pop[member] = pop_index[name]
    with open(path, "w", encoding="utf-8") as handle:
        for sample in ordered:
            row = [str(sample), str(sample_pop.get(sample, 1))]
            for locus in loci:
                allele1, allele2 = _call_alleles(calls.get(locus, {}).get(sample, ""))
                code1, code2 = _allele_code(allele1), _allele_code(allele2)
                row.append(str(code1) if code1 else "-9")
                row.append(str(code2) if code2 else "-9")
            handle.write("\t".join(row) + "\n")
    return len(ordered)


# ===========================================================================
# VCF 4.2 export (write-only; VCF files are not re-imported by the GUI)
# ===========================================================================

#: Virtual contig used for loci that have no unique reference-genome mapping.
VCF_FALLBACK_CONTIG = "amplicon"

#: Genotype string used for a call with no (or a half) allele.
VCF_MISSING_GT = "./."

#: ``<chrom>:<pos><REF>><ALT[,ALT]>`` -- the reference-forward SNP name the
#: analyser generates from a unique BLASTN mapping (``get_genome_snp_name``).
_GENOME_SNP_NAME_RE = re.compile(
    r"^(?P<chrom>.+):(?P<pos>\d+)(?P<ref>[ACGT])>"
    r"(?P<alts>[ACGT](?:,[ACGT])*)$")
_LOCUS_INDEX_RE = re.compile(r"(\d+)\s*$")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _complement_allele(allele: Optional[str]) -> Optional[str]:
    if not allele:
        return allele
    return str(allele).upper().translate(_COMPLEMENT)


def vcf_locus_metadata_from_mapping(loci: Sequence[str],
                                    calls: Dict[str, Dict[str, str]],
                                    mapping: Optional[dict]) -> Dict[str, dict]:
    """Build VCF locus metadata from SNP names and the BLASTN genome mapping.

    Precedence per locus:

    1. reference-forward SNP name ``<chrom>:<pos><REF>><ALT[,ALT]>`` (already
       carries every coordinate);
    2. ``SNP_<consensus position>`` plus a *unique* entry of the mapping
       produced by ``align_consensus_to_reference``.  When the mapping was
       built on the minus strand (``allele_complement_required``) the observed
       alleles are complemented so REF/ALT are reported on the
       reference-forward strand.

    Returns ``{locus: {"chrom", "pos", "ref", "alts", "complement"}}``;
    loci without usable information are simply absent, and :func:`write_vcf`
    then falls back to amplicon coordinates.
    """
    metadata: Dict[str, dict] = {}
    for locus in loci:
        name = str(locus).strip()
        name_match = _GENOME_SNP_NAME_RE.match(name)
        if name_match:
            metadata[locus] = {
                "chrom": name_match.group("chrom"),
                "pos": int(name_match.group("pos")),
                "ref": name_match.group("ref"),
                "alts": name_match.group("alts").split(","),
                "complement": False,
            }
            continue
        index_match = _LOCUS_INDEX_RE.search(name)
        info = None
        if index_match:
            info = (mapping or {}).get(int(index_match.group(1)))
        if not isinstance(info, dict):
            continue
        if str(info.get("mapping_status", "")) != "unique":
            continue
        reference = str(info.get("reference_base", "")).strip().upper()
        position = info.get("genome_position")
        if reference not in "ACGT" or not position:
            continue
        chrom = str(info.get("reference_sequence") or "").split()[0]
        if not chrom:
            continue
        complement_required = bool(info.get("allele_complement_required"))
        observed = set()
        for call in calls.get(locus, {}).values():
            for allele in _call_alleles(call):
                if allele is None:
                    continue
                observed.add(_complement_allele(allele)
                             if complement_required else allele)
        alts = [allele for allele in observed if allele != reference]
        if not alts:
            continue
        metadata[locus] = {
            "chrom": chrom,
            "pos": int(position),
            "ref": reference,
            "alts": sorted(alts, key="ACGT".index),
            "complement": complement_required,
        }
    return metadata


def write_vcf(path: str, samples: Sequence[str], loci: Sequence[str],
              calls: Dict[str, Dict[str, str]],
              populations: Optional[Dict[str, List[str]]] = None,
              locus_metadata: Optional[Dict[str, dict]] = None,
              reference_name: Optional[str] = None,
              source: Optional[str] = None) -> int:
    """Write SNP genotypes as a VCF 4.2 file (one column per sample).

    Coordinates come from ``locus_metadata`` (see
    :func:`vcf_locus_metadata_from_mapping`); when a locus has no reference
    mapping, its reference-forward SNP name ``<chrom>:<pos><REF>><ALT>`` is
    parsed, and as a last resort the record is emitted on the virtual contig
    ``amplicon`` with POS = consensus position (declared in the header).

    Genotypes are encoded ``0/0`` / ``0/1`` / ``1/1`` (unphased, sorted) with
    ``./.`` for a missing call; multi-allelic loci are written as comma
    separated ALT alleles with one AF value per ALT.  Each record carries
    ``NS`` (samples with data) and ``AF`` (alternate allele frequencies) in
    INFO.  Population structure is not part of the VCF specification and is
    ignored here.
    """
    import datetime

    ordered = _dedup(samples)
    metadata = dict(locus_metadata or {})
    fallback_loci: List[str] = []
    contigs: List[str] = []
    contig_lengths: Dict[str, int] = {}
    body: List[str] = []

    for locus in loci:
        observed = sorted(
            {allele
             for call in calls.get(locus, {}).values()
             for allele in _call_alleles(call)
             if allele is not None},
            key="ACGT".index)

        entry = metadata.get(locus)
        if entry and entry.get("chrom") and entry.get("pos"):
            chrom = str(entry["chrom"])
            position = int(entry["pos"])
            reference = str(entry.get("ref", "")).upper()
            alternates = [str(a).upper() for a in (entry.get("alts") or [])
                          if str(a).upper() != str(entry.get("ref", "")).upper()]
            complement = bool(entry.get("complement"))
            name = str(locus)
            if entry.get("length"):
                contig_lengths.setdefault(chrom, int(entry["length"]))
        else:
            name_match = _GENOME_SNP_NAME_RE.match(str(locus).strip())
            if name_match:
                chrom = name_match.group("chrom")
                position = int(name_match.group("pos"))
                reference = name_match.group("ref")
                alternates = name_match.group("alts").split(",")
                complement = False
                name = str(locus)
            else:
                if not observed:
                    continue
                chrom = VCF_FALLBACK_CONTIG
                index_match = _LOCUS_INDEX_RE.search(str(locus).strip())
                position = (int(index_match.group(1)) if index_match
                            else len(body) + 1)
                reference = observed[0]
                alternates = observed[1:]
                complement = False
                name = str(locus)
                fallback_loci.append(locus)
        if chrom not in contigs:
            contigs.append(chrom)

        alleles = [reference] + [a for a in alternates if a and a != reference]
        for allele in observed:
            if complement:
                allele = _complement_allele(allele)
            if allele and allele not in alleles:
                alleles.append(allele)
        allele_index = {allele: idx for idx, allele in enumerate(alleles)}

        genotypes: List[str] = []
        n_called = 0
        total_alleles = 0
        alt_counts = [0] * (len(alleles) - 1)
        for sample in ordered:
            first, second = _call_alleles(calls.get(locus, {}).get(sample, ""))
            if complement:
                first = _complement_allele(first)
                second = _complement_allele(second)
            if (first is None or first not in allele_index or
                    second is None or second not in allele_index):
                genotypes.append(VCF_MISSING_GT)
                continue
            n_called += 1
            total_alleles += 2
            idx1 = allele_index[first]
            idx2 = allele_index[second]
            for idx in (idx1, idx2):
                if idx:
                    alt_counts[idx - 1] += 1
            genotypes.append("/".join(str(v) for v in sorted((idx1, idx2))))

        alt_field = ",".join(alleles[1:]) if len(alleles) > 1 else "."
        if total_alleles:
            frequencies = ",".join(f"{count / total_alleles:.6f}"
                                   for count in alt_counts)
        else:
            frequencies = "0"
        info_field = f"NS={n_called};AF={frequencies}"
        body.append("\t".join([
            chrom, str(position), name, reference, alt_field, ".", "PASS",
            info_field, "GT"] + genotypes))

    lines = [
        "##fileformat=VCFv4.2",
        "##fileDate=" + datetime.date.today().strftime("%Y%m%d"),
        "##source=" + (source or "ABI-SNP_Analyzer"),
    ]
    if reference_name:
        lines.append("##reference=" + str(reference_name).strip())
    lines += [
        '##INFO=<ID=NS,Number=1,Type=Integer,Description='
        '"Number of samples with data">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description='
        '"Alternate allele frequency in the sample set">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    for chrom in contigs:
        if chrom in contig_lengths:
            lines.append(f"##contig=<ID={chrom},length={contig_lengths[chrom]}>")
        else:
            lines.append(f"##contig=<ID={chrom}>")
    if fallback_loci:
        names = ", ".join(str(locus) for locus in fallback_loci[:20])
        more = "" if len(fallback_loci) <= 20 else " ... ({} total)".format(
            len(fallback_loci))
        lines.append(
            '##ABI-SNP_Analyzer_note=Loci without a unique reference-genome '
            'mapping are emitted on the virtual contig "'
            + VCF_FALLBACK_CONTIG + '" with POS = consensus (amplicon) '
            'coordinate: ' + names + more)
    lines.append("\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
         "FORMAT"] + [str(sample) for sample in ordered]))
    lines += body

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return len(ordered)
