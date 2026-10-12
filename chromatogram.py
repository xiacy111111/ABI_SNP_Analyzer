"""
Chromatogram processing module for ABI-SNP Analyzer.

Classes:
    Trace           — Raw ABIF trace data container
    BaseCalls       — Processed base call results
    AbifParser      — ABIF binary format reader and peak detector

IUPAC utilities:
    _IUPAC_TO_BASES        — IUPAC code → base set mapping
    _iupac_code             — Base set → IUPAC code
    _is_reliable_diploid_iupac — Check if IUPAC code is a reliable diploid
    _compute_allele_freqs_from_iupac — Decompose IUPAC codes to allele freqs
    convert_iupac_to_diploid — IUPAC code → diploid genotype string
"""

import struct
from typing import Dict, List, Optional, Tuple
from collections import Counter

import numpy as np


# ---------- Data Containers ----------


class Trace:
    def __init__(self):
        self.leadingGaps = 0
        self.trailingGaps = 0
        self.basecalls1 = ""
        self.basecalls2 = ""
        self.qual = []
        self.basecallpos = []
        self.traceACGT = {
            'A': np.array([]),
            'C': np.array([]),
            'G': np.array([]),
            'T': np.array([]),
        }


class BaseCalls:
    def __init__(self):
        self.consensus = ""
        self.primary = ""
        self.secondary = ""
        self.secDecompose = ""
        self.bcPos = []
        self.estQual = []
        self.heterozygote_info = []


# ---------- ABIF Parser ----------


class AbifParser:
    def __init__(self):
        pass

    @staticmethod
    def _read_bin_str(buffer: bytes, pos: int, length: int) -> str:
        return buffer[pos:pos+length].decode('ascii', errors='ignore').strip('\x00')

    @staticmethod
    def _read_bin_i16(buffer: bytes, pos: int) -> int:
        return struct.unpack('>h', buffer[pos:pos+2])[0]

    @staticmethod
    def _read_bin_i32(buffer: bytes, pos: int) -> int:
        return struct.unpack('>i', buffer[pos:pos+4])[0]

    @staticmethod
    def _read_bin_ui8(buffer: bytes, pos: int) -> int:
        return buffer[pos]

    @staticmethod
    def _replace_non_dna(sequence: str) -> str:
        sequence = sequence.upper()
        return "".join(base if base in "ACGT" else "N" for base in sequence)

    def read_abi(self, filename: str) -> Optional[Trace]:
        """Read an ABIF file with bounds checks and graceful failure."""
        try:
            with open(filename, "rb") as file_handle:
                buffer = file_handle.read()

            if len(buffer) < 30:
                print(f"Invalid ABIF file (header too short): {filename}")
                return None
            if self._read_bin_str(buffer, 0, 4) != "ABIF":
                print(f"File is not in ABIF format: {filename}")
                return None

            entry_size = self._read_bin_i16(buffer, 16)
            entry_count = self._read_bin_i32(buffer, 18)
            directory_offset = self._read_bin_i32(buffer, 26)
            if entry_size < 28 or entry_count < 0 or directory_offset < 0:
                print(f"Invalid ABIF directory metadata: {filename}")
                return None
            if entry_count > 1_000_000:
                print(f"Unreasonable ABIF directory-entry count: {entry_count}")
                return None
            directory_end = directory_offset + entry_count * entry_size
            if directory_offset > len(buffer) or directory_end > len(buffer):
                print(f"ABIF directory extends beyond end of file: {filename}")
                return None

            trace = Trace()
            channel_order = ""
            raw_traces = {9: [], 10: [], 11: [], 12: []}

            for entry_index in range(entry_count):
                entry_offset = directory_offset + entry_index * entry_size
                if entry_offset < 0 or entry_offset + 28 > len(buffer):
                    print(f"Truncated ABIF directory entry {entry_index}: {filename}")
                    return None

                name = self._read_bin_str(buffer, entry_offset, 4)
                number = self._read_bin_i32(buffer, entry_offset + 4)
                element_type = self._read_bin_i16(buffer, entry_offset + 8)
                data_size = self._read_bin_i32(buffer, entry_offset + 16)
                data_offset = self._read_bin_i32(buffer, entry_offset + 20)
                if data_size < 0:
                    print(f"Negative ABIF data size in {name}.{number}")
                    return None

                key = f"{name}.{number}"
                data_start = data_offset if data_size > 4 else entry_offset + 20
                data_end = data_start + data_size
                if data_start < 0 or data_end < data_start or data_end > len(buffer):
                    print(f"ABIF tag {key} points outside the file")
                    return None

                if name == "PCON":
                    element_type = 1

                if element_type == 2:
                    data_string = self._read_bin_str(buffer, data_start, data_size)
                    if key == "PBAS.2":
                        trace.basecalls1 = self._replace_non_dna(data_string)
                    elif key == "P2BA.1":
                        trace.basecalls2 = self._replace_non_dna(data_string)
                    elif key == "FWO_.1":
                        channel_order = data_string.upper()
                elif element_type == 4:
                    if data_size % 2:
                        print(f"Odd byte count for 16-bit ABIF tag {key}")
                        return None
                    values = [
                        self._read_bin_i16(buffer, data_start + index * 2)
                        for index in range(data_size // 2)
                    ]
                    if key == "PLOC.2":
                        trace.basecallpos = values
                    elif key.startswith("DATA.") and number in raw_traces:
                        raw_traces[number] = values
                elif element_type == 1 and key == "PCON.2":
                    trace.qual = [
                        self._read_bin_ui8(buffer, data_start + index)
                        for index in range(data_size)
                    ]

            if len(channel_order) != 4 or set(channel_order) != set("ACGT"):
                if channel_order:
                    print(
                        f"Invalid FWO_ channel order {channel_order!r}; using GATC"
                    )
                channel_order = "GATC"

            for offset, base in enumerate(channel_order):
                trace.traceACGT[base] = np.asarray(raw_traces[9 + offset])

            if not trace.basecalls1 or not trace.basecallpos:
                print(f"ABIF file lacks PBAS.2 or PLOC.2 base calls: {filename}")
                return None
            if any(trace.traceACGT[base].size == 0 for base in "ACGT"):
                print(f"ABIF file lacks one or more DATA.9-12 trace channels: {filename}")
                return None

            trace_length = min(trace.traceACGT[base].size for base in "ACGT")
            for base in "ACGT":
                trace.traceACGT[base] = trace.traceACGT[base][:trace_length]

            base_count = min(len(trace.basecalls1), len(trace.basecallpos))
            valid_positions = 0
            for position in trace.basecallpos[:base_count]:
                if 0 <= position < trace_length:
                    valid_positions += 1
                else:
                    break
            base_count = min(base_count, valid_positions)
            if base_count == 0:
                print(f"ABIF base positions are outside trace bounds: {filename}")
                return None

            trace.basecalls1 = trace.basecalls1[:base_count]
            trace.basecallpos = trace.basecallpos[:base_count]
            if len(trace.qual) < base_count:
                print(
                    f"ABIF quality scores missing/short; padding {filename} with zeros"
                )
                trace.qual = list(trace.qual) + [0] * (base_count - len(trace.qual))
            else:
                trace.qual = list(trace.qual[:base_count])
            return trace
        except (OSError, struct.error, IndexError, TypeError, ValueError) as exc:
            print(f"Failed to parse ABIF file {filename}: {exc}")
            return None

    def _find_peak(self, signal: np.ndarray, start: float, end: float) -> Tuple[int, int]:
        start_i = max(1, int(np.floor(start)))
        end_i = min(len(signal) - 1, int(np.floor(end)))

        if start_i >= end_i:
            return int(np.floor(start)), 0

        best_idx = start_i
        best_val = 0

        for i in range(start_i, end_i):
            if ((signal[i-1] <= signal[i] and signal[i] > signal[i+1]) or
                (signal[i-1] < signal[i] and signal[i] >= signal[i+1])):
                if signal[i] > best_val:
                    best_val = signal[i]
                    best_idx = i

        # Note: we intentionally do NOT fall back to the interval max when
        # no local maximum is found.  Monotonically changing signals within
        # an interval are almost always spillover from an adjacent base
        # (e.g. the tail of the previous base's peak or the rising edge of
        # the next base's peak).  Counting those edge values as genuine
        # peaks produces false heterozygous calls.
        return best_idx, best_val

    @staticmethod
    def _base_intervals(base_positions):
        """Return midpoint-bounded signal intervals for every called base."""
        if not base_positions:
            return [], []
        starts = []
        ends = []
        for index, position in enumerate(base_positions):
            if index == 0:
                left_spacing = (
                    base_positions[1] - position
                    if len(base_positions) > 1
                    else max(position, 1)
                )
            else:
                left_spacing = position - base_positions[index - 1]
            if index == len(base_positions) - 1:
                right_spacing = left_spacing
            else:
                right_spacing = base_positions[index + 1] - position
            starts.append(position - 0.5 * left_spacing)
            ends.append(position + 0.5 * right_spacing)
        return starts, ends

    def _peak_call_state(
        self,
        trace,
        base_index,
        interval_start,
        interval_end,
        signal_ratio,
        minimum_base_quality=5,
    ):
        """Classify one trace position using a shared, intensity-sorted rule."""
        quality = trace.qual[base_index] if base_index < len(trace.qual) else 0
        peaks = {}
        for base in "ACGT":
            peak_index, peak_value = self._find_peak(
                trace.traceACGT[base], interval_start, interval_end
            )
            peaks[base] = {
                "position": peak_index,
                "intensity": max(float(peak_value), 0.0),
            }

        ordered = sorted(
            peaks.items(), key=lambda item: (-item[1]["intensity"], item[0])
        )
        maximum = ordered[0][1]["intensity"] if ordered else 0.0
        for _, peak in ordered:
            peak["ratio"] = peak["intensity"] / maximum if maximum > 0 else 0.0

        if quality < minimum_base_quality:
            status = "low_quality"
            valid = []
        elif maximum <= 0:
            status = "no_signal"
            valid = []
        else:
            valid = [item for item in ordered if item[1]["ratio"] >= signal_ratio]
            # If the dominant peak is extremely weak (< 150 absolute intensity),
            # treat the position as homozygous regardless of secondary-peak ratios.
            # Very weak signals are prone to false heterozygous calls because the
            # true base's peak may be truncated by the interval boundary (leaving
            # only noise-level secondary peaks inside the interval).
            if maximum > 0 and maximum < 150:
                status = "homozygous"
                valid = valid[:1]
            elif len(valid) == 1:
                status = "homozygous"
            elif len(valid) == 2:
                status = "heterozygous"
            else:
                status = "ambiguous_3plus_peaks"

        return {
            "status": status,
            "quality": quality,
            "ordered_peaks": ordered,
            "valid_peaks": valid,
            "maximum_intensity": maximum,
        }

    @staticmethod
    def _filter_false_peaks(all_states, false_peak_threshold=0.2, neighbor_count=5):
        """Filter out false peaks where max signal is too low vs neighbors.

        For each position, compare its maximum peak intensity against the
        average of the *neighbor_count* positions immediately before and
        after.  If the ratio falls below *false_peak_threshold*, reclassify
        the position as ``false_peak`` (no genuine base call).

        This catches noise-only intervals where the trace shows no real
        signal relative to the surrounding true peaks.
        """
        n = len(all_states)
        for i in range(n):
            max_intensity = all_states[i].get("maximum_intensity", 0)
            if max_intensity <= 0:
                continue

            # Collect neighbor max intensities (5 before + 5 after)
            neighbors = []
            for j in range(max(0, i - neighbor_count),
                           min(n, i + neighbor_count + 1)):
                if j == i:
                    continue
                neighbor_max = all_states[j].get("maximum_intensity", 0)
                if neighbor_max > 0:
                    neighbors.append(neighbor_max)

            if len(neighbors) < 3:
                continue  # Not enough neighbours for a reliable comparison

            avg_neighbor = sum(neighbors) / len(neighbors)
            if avg_neighbor > 0 and max_intensity < avg_neighbor * false_peak_threshold:
                all_states[i]["status"] = "false_peak"
                all_states[i]["valid_peaks"] = []

        return all_states

    def detect_heterozygotes(
        self,
        tr: Trace,
        sigratio: float = 0.3,
        trim_bp: int = 20,
        minimum_base_quality: int = 5,
    ) -> List[Dict]:
        """Detect two-peak heterozygotes and flag three-plus-peak ambiguity."""
        if not 0.0 <= sigratio <= 1.0:
            raise ValueError("sigratio must be between 0 and 1")
        if trim_bp < 0 or minimum_base_quality < 0:
            raise ValueError("trim_bp and minimum_base_quality must be non-negative")

        starts, ends = self._base_intervals(tr.basecallpos)
        total_positions = len(starts)

        # --- Pass 1: compute peak-call state for EVERY position ---
        all_states = []
        for index in range(total_positions):
            state = self._peak_call_state(
                tr,
                index,
                starts[index],
                ends[index],
                sigratio,
                minimum_base_quality,
            )
            all_states.append(state)

        # --- Pass 2: filter out false peaks (noise-only positions) ---
        all_states = self._filter_false_peaks(all_states,
                                              false_peak_threshold=0.2,
                                              neighbor_count=5)

        # --- Pass 3: collect detections (skip trim_bp edges) ---
        detections = []
        for index in range(total_positions):
            if index < trim_bp or index >= total_positions - trim_bp:
                continue
            state = all_states[index]
            if state["status"] not in {
                "heterozygous",
                "ambiguous_3plus_peaks",
            }:
                continue

            valid = state["valid_peaks"]
            bases = [base for base, _ in valid]
            detections.append(
                {
                    "position": index,
                    "peak_position": tr.basecallpos[index],
                    "bases": bases,
                    "intensities": [peak["intensity"] for _, peak in valid],
                    "ratios": [peak["ratio"] for _, peak in valid],
                    "iupac_code": _iupac_code(bases),
                    "is_heterozygote": state["status"] == "heterozygous",
                    "qc_status": (
                        "ok"
                        if state["status"] == "heterozygous"
                        else "ambiguous_3plus_peaks"
                    ),
                    "base_quality": state["quality"],
                    "primary_base": bases[0] if bases else "",
                    "secondary_base": bases[1] if len(bases) > 1 else "",
                }
            )
        return detections

    def generate_iupac_sequence(self, tr: Trace, heterozygotes: List[Dict], trim_bp: int = 20) -> str:
        """
        Generate IUPAC-coded sequence from basecalls and heterozygote detections.

        [FIX v1.1] For positions with 3+ peaks (qc_status='ambiguous_3plus_peaks'),
        the base is left as the original basecall (not converted to IUPAC) because
        the three-peak code (V/H/D/B) is not reliably interpretable as a diploid genotype.

        [FIX v1.2] trim_bp bases at each end are set to 'N' to exclude edge
        positions from downstream SNP analysis.
        """
        sequence = list(tr.basecalls1)
        total_positions = len(sequence)

        # Mask trim_bp edges to N
        for i in range(trim_bp):
            if i < total_positions:
                sequence[i] = 'N'
        for i in range(total_positions - trim_bp, total_positions):
            if i >= 0:
                sequence[i] = 'N'

        for het in heterozygotes:
            pos = het['position']
            if pos < len(sequence):
                qc = het.get('qc_status', 'ok')
                if qc == 'ambiguous_3plus_peaks':
                    # Three-plus-peak calls are not reliable diploid genotypes.
                    sequence[pos] = 'N'
                    continue
                iupac_code = het['iupac_code']
                sequence[pos] = iupac_code

        return ''.join(sequence)

    def basecall(
        self,
        tr: Trace,
        sigratio: float = 0.3,
        minimum_base_quality: int = 5,
        trim_bp: int = 20,
    ) -> BaseCalls:
        """Generate primary/secondary calls using the shared peak classifier.

        trim_bp: number of bases at each end to exclude from SNP analysis
                 (set to 'N' in the output consensus).
        """
        if not 0.0 <= sigratio <= 1.0:
            raise ValueError("sigratio must be between 0 and 1")
        starts, ends = self._base_intervals(tr.basecallpos)
        total_positions = len(starts)

        # --- Pass 1: compute peak-call state for EVERY position ---
        all_states = []
        for index in range(total_positions):
            state = self._peak_call_state(
                tr,
                index,
                starts[index],
                ends[index],
                sigratio,
                minimum_base_quality,
            )
            all_states.append(state)

        # --- Pass 2: filter out false peaks (noise-only positions) ---
        all_states = self._filter_false_peaks(all_states,
                                              false_peak_threshold=0.2,
                                              neighbor_count=5)

        # --- Pass 3: generate calls ---
        calls = BaseCalls()
        for index in range(total_positions):
            state = all_states[index]
            valid = state["valid_peaks"]
            selected_position = (
                valid[0][1]["position"]
                if valid
                else int((starts[index] + ends[index]) / 2)
            )

            # Exclude trim_bp edges and false peaks from analysis
            if (index < trim_bp or index >= total_positions - trim_bp
                    or state["status"] == "false_peak"):
                calls.primary += "N"
                calls.secondary += "N"
                calls.consensus += "N"
                calls.heterozygote_info.append((selected_position, [], []))
            elif state["status"] == "homozygous":
                base = valid[0][0]
                calls.primary += base
                calls.secondary += base
                calls.consensus += base
                calls.heterozygote_info.append((selected_position, [base], [1.0]))
            elif state["status"] == "heterozygous":
                bases = [base for base, _ in valid]
                ratios = [peak["ratio"] for _, peak in valid]
                calls.primary += bases[0]
                calls.secondary += _iupac_code(bases)
                calls.consensus += "N"
                calls.heterozygote_info.append((selected_position, bases, ratios))
            else:
                calls.primary += "N"
                calls.secondary += "N"
                calls.consensus += "N"
                calls.heterozygote_info.append((selected_position, [], []))
            calls.bcPos.append(selected_position)

        calls.estQual = list(tr.qual[: len(calls.primary)])
        if len(calls.estQual) < len(calls.primary):
            calls.estQual.extend([0] * (len(calls.primary) - len(calls.estQual)))
        return calls

    def _estimate_qualities_simple(self, bc: BaseCalls, tr: Trace):
        bc.estQual = [30] * len(bc.primary)


# ---------- IUPAC Utilities ----------


def _iupac_code(bases: List[str]) -> str:
    iupac_map = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['C', 'G']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['A', 'C', 'G']): 'V', frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'G', 'T']): 'D', frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }
    key = frozenset(bases)
    return iupac_map.get(key, 'N')


# IUPAC code to set of bases mapping (for decomposing heterozygous codes)
_IUPAC_TO_BASES = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['C', 'G'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'],
    'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'],
    'D': ['A', 'G', 'T'], 'B': ['C', 'G', 'T'],
    'N': [], '-': []
}


def _is_reliable_diploid_iupac(code: str) -> bool:
    """[FIX v1.1] Check if IUPAC code represents a reliable diploid genotype.
    Only standard two-allele IUPAC codes (R/Y/S/W/K/M) and single bases are
    considered reliable for diploid genotype calling.
    Three-allele codes (V/H/D/B) and N are not reliable diploid genotypes.
    """
    return code in {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M'}


def _compute_allele_freqs_from_iupac(bases_list):
    """Compute per-allele frequencies by decomposing IUPAC codes."""
    counter = Counter()
    total = 0
    for b in bases_list:
        if _is_reliable_diploid_iupac(b):
            alleles = _IUPAC_TO_BASES.get(b, [])
            for a in alleles:
                counter[a] += 1
                total += 1
    if total == 0:
        return {}
    return {k: v / total for k, v in counter.items()}


def convert_iupac_to_diploid(iupac_code, allele_freqs=None):
    """Convert IUPAC code to diploid genotype string."""
    diploid_map = {
        'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT',
        'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
        'K': 'GT', 'M': 'AC'
    }
    if iupac_code in diploid_map:
        return diploid_map[iupac_code]
    return '--'
