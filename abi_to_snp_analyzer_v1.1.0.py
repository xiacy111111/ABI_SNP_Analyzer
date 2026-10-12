# ========== Required fix for EXE packaging ==========
import sys
import os
import hashlib
import importlib.metadata
import json
import platform
import time
if getattr(sys, "frozen", False):
    # PyInstaller's Matplotlib hook locates mpl-data; do not override MATPLOTLIBDATA.
    os.environ["MPLBACKEND"] = "TkAgg"

# [FIX v1.1] Use TkAgg backend consistently for GUI+savefig compatibility.
# Agg backend prevents FigureCanvasTkAgg from working in the frequency plot tab.
# All plot saving uses savefig() which works with any backend.
import matplotlib
matplotlib.use("TkAgg")
# ====================================
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
from typing import List, Tuple, Dict, Optional
import os
import glob
from datetime import datetime
import re
import gc
import subprocess
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter
import threading
import queue
import csv
from collections import defaultdict
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

# Population genetics module
try:
    from population_genetics import PopulationGenetics, LocusSummary, AMOVAResult, GenotypeFreqTestResult
    POPGEN_AVAILABLE = True
except ImportError:
    POPGEN_AVAILABLE = False

# SNP genotype import/export (Arlequin / GENEPOP / hierfstat / GenAlEx / ...)
try:
    import snp_format_io as snpio
    SNP_IO_AVAILABLE = True
except ImportError:
    snpio = None
    SNP_IO_AVAILABLE = False

# Linkage disequilibrium / one-SNP-per-amplicon-block analysis
try:
    import linkage_analysis as linkage
    LINKAGE_AVAILABLE = True
except ImportError:
    linkage = None
    LINKAGE_AVAILABLE = False

# Get the base runtime path
def get_runtime_paths(
    frozen=None,
    executable=None,
    module_file=None,
    meipass=None,
):
    """Return (application_directory, bundle_directory).

    External executables belong beside the installed application, while files
    collected by PyInstaller are read from the temporary bundle directory.
    """
    if frozen is None:
        frozen = bool(getattr(sys, "frozen", False))
    executable = executable or sys.executable
    module_file = module_file or __file__
    meipass = meipass or getattr(sys, "_MEIPASS", None)

    if frozen:
        application_directory = os.path.dirname(os.path.abspath(executable))
        bundle_directory = os.path.abspath(meipass or application_directory)
    else:
        application_directory = os.path.dirname(os.path.abspath(module_file))
        bundle_directory = application_directory
    return application_directory, bundle_directory


def find_external_tool(file_name, application_directory=None):
    """Find a separately distributed executable near the application."""
    application_directory = application_directory or APP_DIR
    candidates = [
        os.path.join(application_directory, file_name),
        os.path.join(application_directory, "bin", file_name),
        os.path.join(application_directory, "blast", "bin", file_name),
    ]
    if file_name.lower() == "blastn.exe":
        candidates.extend(
            sorted(
                glob.glob(
                    os.path.join(
                        application_directory,
                        "ncbi-blast*",
                        "bin",
                        file_name,
                    )
                )
            )
        )
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return ""


APP_DIR, BUNDLE_DIR = get_runtime_paths()
# Backward-compatible alias for internal bundled resources only.
base_path = BUNDLE_DIR

# [FIX v1.1] Software version
SOFTWARE_VERSION = "1.1.0"


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


    @staticmethod
    @staticmethod
    def _replace_non_dna(sequence: str) -> str:
        sequence = sequence.upper()
        return "".join(base if base in "ACGT" else "N" for base in sequence)


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
            if len(valid) == 1:
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
                    "iupac_code": self._iupac_code(bases),
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
                calls.secondary += self._iupac_code(bases)
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

    @staticmethod
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


def list_ab1_files(directory):
    ab1_files = []
    for file in os.listdir(directory):
        if file.lower().endswith('.ab1'):
            ab1_files.append(os.path.join(directory, file))
    return ab1_files


def plot_raw_trace_data(trace_data: Trace, filename: str, heterozygotes: List[Dict] = None, output_dir: str = None, position_map: dict = None):
    """Plot raw trace data with optional consensus-position mapping.

    position_map: {1-based_sample_pos: 1-based_consensus_pos}.  When provided,
    x-axis tick labels and heterozygous position markers use consensus positions
    instead of sample positions.

    [FIX v1.0.0] Uses Figure + FigureCanvasAgg (headless Agg backend) instead
    of pyplot to avoid Windows GDI handle exhaustion with hundreds of samples.
    """
    if not trace_data.traceACGT['A'].size:
        print("No raw trace data available for plotting")
        return

    num_bases = len(trace_data.basecalls1)
    num_rows = (num_bases + 99) // 100

    fig = Figure(figsize=(15, 5*num_rows))
    _canvas = FigureCanvasAgg(fig)  # headless — zero GDI handles

    axes = [fig.add_subplot(num_rows, 1, r + 1) for r in range(num_rows)]

    for row in range(num_rows):
        start_pos = row * 100
        end_pos = min((row + 1) * 100, num_bases)

        ax = axes[row]

        start_signal = trace_data.basecallpos[start_pos] if start_pos < len(trace_data.basecallpos) else 0
        end_signal = trace_data.basecallpos[end_pos - 1] if end_pos <= len(trace_data.basecallpos) else len(trace_data.traceACGT['A'])

        signal_start_idx = max(0, int(start_signal) - 50)
        signal_end_idx = min(len(trace_data.traceACGT['A']), int(end_signal) + 50)

        x_range = range(signal_start_idx, signal_end_idx)
        if signal_end_idx > signal_start_idx:
            ax.plot(x_range, trace_data.traceACGT['A'][signal_start_idx:signal_end_idx], label='A', color='green', alpha=0.7)
            ax.plot(x_range, trace_data.traceACGT['C'][signal_start_idx:signal_end_idx], label='C', color='blue', alpha=0.7)
            ax.plot(x_range, trace_data.traceACGT['G'][signal_start_idx:signal_end_idx], label='G', color='black', alpha=0.7)
            ax.plot(x_range, trace_data.traceACGT['T'][signal_start_idx:signal_end_idx], label='T', color='red', alpha=0.7)

        for i in range(start_pos, end_pos):
            if i < len(trace_data.basecallpos) and i < len(trace_data.basecalls1):
                pos = trace_data.basecallpos[i]
                base = trace_data.basecalls1[i]
                color_map = {'A': 'green', 'C': 'blue', 'G': 'black', 'T': 'red'}
                color = color_map.get(base, 'gray')

                is_heterozygote = False
                if heterozygotes:
                    for het in heterozygotes:
                        if het['position'] == i:
                            is_heterozygote = True
                            break

                ax.axvline(x=pos, color=color, linestyle='--', alpha=0.5)

                # [FIX v1.1] Safe access: handle N, gap, and non-ACGT bases
                y_pos = 0
                if base in trace_data.traceACGT:
                    arr = trace_data.traceACGT[base]
                    if int(pos) < len(arr):
                        y_pos = arr[int(pos)]

                text = base
                bbox_color = "white"
                if is_heterozygote:
                    text = f"{base}*"
                    bbox_color = "yellow"

                ax.text(pos, y_pos, text, fontsize=8, ha='center', va='bottom',
                       bbox=dict(boxstyle="round,pad=0.1", facecolor=bbox_color, alpha=0.8))

        base_positions = trace_data.basecallpos[start_pos:end_pos]
        tick_positions = []
        tick_labels = []
        for i, pos in enumerate(base_positions):
            if i % 10 == 0:
                tick_positions.append(pos)
                seq_pos_1based = start_pos + i + 1
                label = seq_pos_1based
                if position_map is not None:
                    label = position_map.get(seq_pos_1based, seq_pos_1based)
                tick_labels.append(label)

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', rotation=45)

        if position_map is not None:
            _s, _e = start_pos + 1, end_pos
            _cs = position_map.get(_s, _s)
            _ce = position_map.get(_e, _e)
            ax.set_title(f'Raw Trace Data - Row {row+1} (Consensus Bases {_cs}-{_ce})')
        else:
            ax.set_title(f'Raw Trace Data - Row {row+1} (Bases {start_pos+1}-{end_pos})')
        ax.set_xlabel('Base Number (Consensus Position)' if position_map is not None else 'Base Number')
        ax.set_ylabel('Signal Intensity')
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'Raw Trace Data (A/C/G/T Signals) - {os.path.basename(filename)}\n(*) Indicates heterozygous positions', fontsize=14)

    fig.tight_layout()

    if output_dir:
        basename = os.path.splitext(os.path.basename(filename))[0]
        img_filename = os.path.join(output_dir, f"{basename}_trace_plot_combined.png")
        _canvas.print_figure(img_filename, dpi=300, bbox_inches='tight')

    # Agg backend: just clearing Python objects — no GDI handles to leak.
    del _canvas
    plt.close("all")


def save_fasta_file_with_iupac(trace_data: Trace, iupac_sequence: str, filename: str, output_dir: str):
    basename = os.path.splitext(os.path.basename(filename))[0]
    fasta_filename = os.path.join(output_dir, f"{basename}_sequence_with_iupac.fasta")

    with open(fasta_filename, 'w') as f:
        f.write(f">Forward_sequence_with_IUPAC_{basename}\n")
        f.write(f"{iupac_sequence}\n")

    print(f"  -> IUPAC FASTA file saved: {fasta_filename}")


def save_heterozygote_report(heterozygotes: List[Dict], filename: str, output_dir: str, position_map: dict = None):
    """Save heterozygote report, optionally remapping positions to consensus.

    position_map: {0-based_sample_pos: consensus_1based_pos}.  When provided,
    Position is rewritten to the consensus coordinate; Sample_Position records
    the original 0-based sample position.
    """
    # [FIX v1.1] Always generate report (even empty) for traceability
    basename = os.path.splitext(os.path.basename(filename))[0]
    report_filename = os.path.join(output_dir, f"{basename}_heterozygote_report.xlsx")

    columns = [
        'Position', 'Sample_Position', 'Peak_Position', 'Detected_Bases', 'Signal_Intensities',
        'Signal_Ratios', 'Base_Quality', 'IUPAC_Code', 'Is_Heterozygote', 'QC_Status'
    ] if position_map is not None else [
        'Position', 'Peak_Position', 'Detected_Bases', 'Signal_Intensities',
        'Signal_Ratios', 'Base_Quality', 'IUPAC_Code', 'Is_Heterozygote', 'QC_Status'
    ]

    if not heterozygotes:
        # [FIX v1.1] Generate empty report instead of skipping
        empty_df = pd.DataFrame(columns=columns)
        empty_df.to_excel(report_filename, index=False)
        print(f"  -> No heterozygotes detected, empty report saved: {report_filename}")
        return

    report_data = {col: [] for col in columns}

    for het in heterozygotes:
        sample_pos_0based = het['position']
        sample_pos_1based = sample_pos_0based + 1
        if position_map is not None:
            consensus_pos = position_map.get(sample_pos_1based, sample_pos_1based)
            report_data['Position'].append(consensus_pos)
            report_data['Sample_Position'].append(sample_pos_0based)
        else:
            report_data['Position'].append(sample_pos_0based)
        report_data['Peak_Position'].append(het['peak_position'])
        report_data['Detected_Bases'].append('/'.join(het['bases']))
        report_data['Signal_Intensities'].append(str(het['intensities']))
        report_data['Signal_Ratios'].append(str([f"{r:.2f}" for r in het['ratios']]))
        report_data['Base_Quality'].append(het.get('base_quality', ''))
        report_data['IUPAC_Code'].append(het['iupac_code'])
        report_data['Is_Heterozygote'].append(het['is_heterozygote'])
        report_data['QC_Status'].append(het.get('qc_status', 'ok'))  # [FIX v1.1]

    df = pd.DataFrame(report_data)
    df.to_excel(report_filename, index=False)

    print(f"  -> Heterozygote report saved: {report_filename}")
    print(f"  -> Detected {len(heterozygotes)} heterozygous positions")


def combine_iupac_fastas(output_dir: str):
    iupac_pattern = os.path.join(output_dir, "*_sequence_with_iupac.fasta")
    iupac_files = glob.glob(iupac_pattern, recursive=True)

    if not iupac_files:
        print("No IUPAC FASTA files found to merge")
        return None

    combined_filename = os.path.join(output_dir, "combined_all_forward_sequences_with_iupac.fasta")

    with open(combined_filename, 'w') as outfile:
        for file_path in sorted(iupac_files):
            with open(file_path, 'r') as infile:
                content = infile.read()
                outfile.write(content)
                outfile.write("\n")

    print(f"  -> Merged {len(iupac_files)} IUPAC FASTA files into: {combined_filename}")
    return combined_filename


# [FIX v1.1] Detect MUSCLE version and use appropriate command syntax
def _detect_muscle_version(muscle_exe: str) -> Tuple[int, str]:
    """Detect MUSCLE version. Returns (major_version, version_string)."""
    try:
        result = subprocess.run([muscle_exe, "-version"], capture_output=True, text=True, timeout=10)
        output = result.stdout + result.stderr
        # Try to find version number like "5.1" or "3.8"
        match = re.search(r'(\d+)\.(\d+)', output)
        if match:
            major = int(match.group(1))
            version_str = match.group(0)
        else:
            major = 3  # Default to v3 syntax
            version_str = "unknown"
        return major, version_str
    except Exception:
        return 3, "unknown (defaulting to v3 syntax)"


def run_muscle_alignment(input_fasta, output_aln, muscle_exe):
    """Run MUSCLE multiple sequence alignment.

    [FIX v1.1] Detects MUSCLE version and uses appropriate command syntax:
    - MUSCLE v5+: -align input -output alignment
    - MUSCLE v3: -in input -out alignment
    """
    try:
        major, ver_str = _detect_muscle_version(muscle_exe)
        print(f"MUSCLE version detected: {ver_str} (major={major})")

        if major >= 5:
            cmd = [muscle_exe, "-align", input_fasta, "-output", output_aln]
        else:
            cmd = [muscle_exe, "-in", input_fasta, "-out", output_aln]

        print(f"MUSCLE command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"MUSCLE failed (rc={result.returncode}): {result.stderr}")
            # Fallback: try v3 syntax if v5 command failed
            if major >= 5:
                print("Trying v3 syntax as fallback...")
                cmd_v3 = [muscle_exe, "-in", input_fasta, "-out", output_aln]
                result = subprocess.run(cmd_v3, capture_output=True, text=True, timeout=300)
                if result.returncode != 0:
                    print(f"MUSCLE v3 fallback also failed: {result.stderr}")
                    return False

        print("MUSCLE alignment complete")
        return True
    except subprocess.TimeoutExpired:
        print("MUSCLE alignment timed out (300s)")
        return False
    except FileNotFoundError:
        print(f"MUSCLE program not found: {muscle_exe}")
        return False


# [FIX v1.1] IUPAC compatibility mapping for sequence similarity calculation
_IUPAC_COMPAT = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'C', 'G'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'},
    'V': {'A', 'C', 'G'}, 'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'}, 'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'},
}

def _bases_compatible(b1: str, b2: str) -> bool:
    """Check if two bases (possibly IUPAC) are compatible.
    E.g., R (A/G) is compatible with A, G, and R."""
    set1 = _IUPAC_COMPAT.get(b1, set())
    set2 = _IUPAC_COMPAT.get(b2, set())
    if not set1 or not set2:
        return False
    return bool(set1 & set2)  # Intersection non-empty → compatible


def calculate_sequence_similarity(
    alignment_file, similarity_threshold, minimum_overlap_coverage=0.8
):
    """Calculate IUPAC-aware identity and usable overlap coverage.

    Usable overlap coverage is the number of columns where both sequences have
    non-gap, non-N symbols divided by the larger non-gap sequence span. A pair
    with high identity but only a short usable overlap therefore cannot pass.
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        sequence_count = len(alignment)
        alignment_length = alignment.get_alignment_length()
        print(
            f"Alignment loaded: {sequence_count} sequences, length "
            f"{alignment_length}"
        )

        if sequence_count < 2:
            excluded = []
            for record in alignment:
                excluded.append(
                    {
                        "id": record.id,
                        "seq": str(record.seq).upper(),
                        "similarity": 0.0,
                        "overlap_coverage": 0.0,
                        "reason": "at least two sequences are required",
                    }
                )
            print("At least two sequences are required for similarity filtering")
            return [], excluded

        sequence_ids = [record.id for record in alignment]
        identity_matrix = {sequence_id: {} for sequence_id in sequence_ids}
        coverage_matrix = {sequence_id: {} for sequence_id in sequence_ids}

        for i in range(sequence_count):
            sequence1 = str(alignment[i].seq).upper()
            id1 = alignment[i].id
            span1 = sum(base != "-" for base in sequence1)
            identity_matrix[id1][id1] = 1.0
            coverage_matrix[id1][id1] = 1.0

            for j in range(i + 1, sequence_count):
                sequence2 = str(alignment[j].seq).upper()
                id2 = alignment[j].id
                span2 = sum(base != "-" for base in sequence2)
                comparable = 0
                compatible = 0

                for base1, base2 in zip(sequence1, sequence2):
                    if base1 in {"-", "N"} or base2 in {"-", "N"}:
                        continue
                    comparable += 1
                    if _bases_compatible(base1, base2):
                        compatible += 1

                identity = compatible / comparable if comparable else 0.0
                larger_span = max(span1, span2)
                coverage = comparable / larger_span if larger_span else 0.0
                identity_matrix[id1][id2] = identity
                identity_matrix[id2][id1] = identity
                coverage_matrix[id1][id2] = coverage
                coverage_matrix[id2][id1] = coverage

        retained = []
        excluded = []
        for record in alignment:
            sequence_id = record.id
            other_ids = [item for item in sequence_ids if item != sequence_id]
            identities = [identity_matrix[sequence_id][item] for item in other_ids]
            coverages = [coverage_matrix[sequence_id][item] for item in other_ids]
            average_identity = sum(identities) / len(identities)
            average_coverage = sum(coverages) / len(coverages)

            reasons = []
            if average_identity < similarity_threshold:
                reasons.append(
                    f"avg_identity={average_identity:.4f} < "
                    f"threshold={similarity_threshold:.4f}"
                )
            if average_coverage < minimum_overlap_coverage:
                reasons.append(
                    f"avg_usable_overlap={average_coverage:.4f} < "
                    f"threshold={minimum_overlap_coverage:.4f}"
                )

            result = {
                "id": sequence_id,
                "seq": str(record.seq).upper(),
                "similarity": average_identity,
                "overlap_coverage": average_coverage,
            }
            if reasons:
                result["reason"] = "; ".join(reasons)
                excluded.append(result)
            else:
                retained.append(result)

            print(
                f"Sequence {sequence_id}: average compatible identity="
                f"{average_identity:.4f}, average usable overlap="
                f"{average_coverage:.4f}, status="
                f"{'excluded' if reasons else 'retained'}"
            )

        print(
            f"Similarity filter: retained {len(retained)}, excluded "
            f"{len(excluded)}; identity threshold={similarity_threshold}, "
            f"usable-overlap threshold={minimum_overlap_coverage}"
        )
        if excluded:
            print(f"Excluded: {[(item['id'], item['reason']) for item in excluded]}")
        return retained, excluded
    except Exception as exc:
        print(f"Error calculating sequence similarity: {exc}")
        import traceback

        traceback.print_exc()
        return [], []




def create_filtered_fasta(similar_sequences, original_fasta, filtered_fasta):
    original_records = list(SeqIO.parse(original_fasta, "fasta"))

    keep_ids = {s['id'] for s in similar_sequences}

    filtered_records = [rec for rec in original_records if rec.id in keep_ids]

    SeqIO.write(filtered_records, filtered_fasta, "fasta")
    print(f"Filtered FASTA created: {filtered_fasta}, {len(filtered_records)} sequences")
    return filtered_records


# [FIX v1.1] IUPAC-aware consensus sequence generation
def create_consensus_sequence(alignment_file, output_fasta, quality_threshold=0.2):
    """
    Generate consensus sequence from multiple sequence alignment.

    [FIX v1.1] IUPAC-aware consensus: decomposes IUPAC ambiguity codes into
    their constituent alleles (e.g., R → A+G) and counts each allele with
    fractional weight (e.g., R contributes 0.5 to A and 0.5 to G).
    The base with the highest weighted count is selected as consensus.
    N and gap characters are excluded from counting.
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        print(f"Read alignment: {len(alignment)} sequences, length {alignment.get_alignment_length()}")

        consensus = []
        positions = alignment.get_alignment_length()
        total_sequences = len(alignment)
        min_valid_sequences = total_sequences * quality_threshold

        for pos in range(positions):
            # [FIX v1.1] IUPAC-aware allele counting with fractional weights
            base_weights = Counter()
            for record in alignment:
                base = str(record.seq[pos]).upper()
                if base in ['-', 'N']:
                    continue
                if base in _IUPAC_COMPAT:
                    alleles = list(_IUPAC_COMPAT[base])
                    # Distribute weight evenly among constituent alleles
                    weight = 1.0 / len(alleles) if alleles else 0.0
                    for a in alleles:
                        base_weights[a] += weight
                else:
                    base_weights[base] += 1.0

            valid_bases_count = sum(base_weights.values())
            if valid_bases_count >= min_valid_sequences:
                most_common_base = base_weights.most_common(1)[0][0]
                consensus.append(most_common_base)
            else:
                consensus.append('N')

        consensus_seq = ''.join(consensus)
        print(f"Consensus sequence length: {len(consensus_seq)}")
        print(f"First 100 bases: {consensus_seq[:100]}")

        with open(output_fasta, 'w') as f:
            f.write(f">Consensus_Sequence\n{consensus_seq}\n")

        print(f"Consensus sequence saved to: {output_fasta}")
        return consensus_seq
    except Exception as e:
        print(f"Error creating consensus sequence: {e}")
        import traceback
        traceback.print_exc()
        return None


def calculate_minor_allele_freq(genotypes: List[str]) -> float:
    """Calculate MAF from complete A/C/G/T diploid genotype strings."""
    allele_counts = Counter()
    for genotype in genotypes or []:
        if genotype is None:
            continue
        text = str(genotype).strip().upper()
        if len(text) < 2 or text in {"--", "NN", "N/A", "NA", "NAN"}:
            continue
        allele1, allele2 = text[0], text[1]
        if allele1 not in "ACGT" or allele2 not in "ACGT":
            continue
        allele_counts[allele1] += 1
        allele_counts[allele2] += 1

    total_alleles = sum(allele_counts.values())
    if total_alleles == 0 or len(allele_counts) < 2:
        return 0.0
    return min(count / total_alleles for count in allele_counts.values())



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
    """Compute per-allele frequencies by decomposing IUPAC codes.

    For each base in bases_list, decompose the IUPAC code into individual
    alleles (A/C/G/T) and count. Returns dict {allele: frequency}.
    [FIX v1.1] Only decomposes reliable diploid IUPAC codes; ignores
    tri-allelic and N codes.
    """
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
    """Convert IUPAC code to diploid genotype string.

    [FIX v1.1] Only converts reliable diploid codes (A/C/G/T/R/Y/S/W/K/M).
    Tri-allelic codes (V/H/D/B) return '--' (missing) as they cannot be
    reliably resolved to a diploid genotype from chromatogram data alone.
    """
    diploid_map = {
        'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT',
        'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
        'K': 'GT', 'M': 'AC'
    }
    if iupac_code in diploid_map:
        return diploid_map[iupac_code]

    # [FIX v1.1] N, -, V, H, D, B: cannot reliably call diploid genotype
    return '--'

def _decode_alignment_symbol(symbol: str) -> Optional[str]:
    """Decode one alignment symbol into a reliable diploid genotype string."""
    if symbol is None:
        return None
    code = str(symbol).strip().upper()
    genotype = convert_iupac_to_diploid(code)
    return None if genotype == "--" else genotype


def _analyze_snp_column(alignment, retained_ids, position):
    """Return allele/genotype statistics for one alignment column.

    SNP discovery is based on decomposed A/C/G/T alleles, not on the number of
    distinct IUPAC symbols in the alignment.
    """
    sample_genotypes = {}
    symbol_counts = Counter()
    allele_counts = Counter()
    genotype_counts = Counter()
    heterozygote_count = 0

    for record in alignment:
        if record.id not in retained_ids:
            continue
        symbol = str(record.seq[position]).upper()
        symbol_counts[symbol] += 1
        genotype = _decode_alignment_symbol(symbol)
        sample_genotypes[record.id] = genotype or "--"
        if genotype is None:
            continue
        allele1, allele2 = genotype[0], genotype[1]
        allele_counts[allele1] += 1
        allele_counts[allele2] += 1
        canonical = "".join(sorted((allele1, allele2)))
        genotype_counts[canonical] += 1
        if allele1 != allele2:
            heterozygote_count += 1

    total_sequences = len(retained_ids)
    called_genotypes = sum(genotype_counts.values())
    missing_genotypes = total_sequences - called_genotypes
    call_rate = called_genotypes / total_sequences if total_sequences else 0.0
    heterozygote_rate = (
        heterozygote_count / called_genotypes if called_genotypes else 0.0
    )
    total_alleles = sum(allele_counts.values())
    allele_frequencies = {
        allele: count / total_alleles for allele, count in allele_counts.items()
    } if total_alleles else {}
    maf = min(allele_frequencies.values()) if len(allele_frequencies) >= 2 else 0.0

    return {
        "sample_genotypes": sample_genotypes,
        "symbol_counts": symbol_counts,
        "allele_counts": allele_counts,
        "genotype_counts": genotype_counts,
        "observed_alleles": sorted(allele_counts, key="ACGT".index),
        "total_sequences": total_sequences,
        "called_genotypes": called_genotypes,
        "missing_genotypes": missing_genotypes,
        "call_rate": call_rate,
        "heterozygote_count": heterozygote_count,
        "heterozygote_rate": heterozygote_rate,
        "maf": maf,
    }



def create_snp_analysis_report(
    snp_data,
    output_excel,
    min_observed_heterozygosity=0.1,
    genome_snp_mapping=None,
    ref_seq_name=None,
    snp_genotype_matrix=None,
):
    """Write candidate SNP statistics and explicit filtering decisions."""
    del min_observed_heterozygosity, genome_snp_mapping, ref_seq_name
    del snp_genotype_matrix  # retained for backward API compatibility

    columns = [
        "SNP_ID",
        "Position",
        "Genome_Position",
        "Consensus_Base",
        "Reference_Base",
        "Sequencing_Orientation_Alleles",
        "Reference_Orientation_Alleles",
        "Mapping_Status",
        "Mapping_Strand",
        "Mapping_Identity",
        "Mapping_Query_Coverage",
        "Mapping_Evalue",
        "Mapping_Bitscore",
        "Alleles",
        "Total_Sequences",
        "Valid_Genotypes_Count",
        "Missing_Genotypes_Count",
        "Call_Rate",
        "Percentage_Valid",
        "Allele_Counts",
        "Genotype_Counts",
        "Alignment_Symbol_Counts",
        "Minor_Allele_Frequency",
        "Heterozygote_Rate",
        "Pass_Call_Rate",
        "Pass_MAF",
        "Pass_Min_Heterozygosity",
        "Retained",
    ]

    rows = []
    for snp in snp_data or []:
        rows.append(
            {
                "SNP_ID": snp.get("SNP_ID", ""),
                "Position": snp.get("Position", ""),
                "Genome_Position": snp.get("Genome_Position", ""),
                "Consensus_Base": snp.get("Consensus_Base", ""),
                "Reference_Base": snp.get("Reference_Base", ""),
                "Sequencing_Orientation_Alleles": snp.get("Sequencing_Alleles", ""),
                "Reference_Orientation_Alleles": snp.get("Reference_Alleles", ""),
                "Mapping_Status": snp.get("Mapping_Status", "unmapped"),
                "Mapping_Strand": snp.get("Mapping_Strand", ""),
                "Mapping_Identity": snp.get("Mapping_Identity", ""),
                "Mapping_Query_Coverage": snp.get("Mapping_Query_Coverage", ""),
                "Mapping_Evalue": snp.get("Mapping_Evalue", ""),
                "Mapping_Bitscore": snp.get("Mapping_Bitscore", ""),
                "Alleles": snp.get("All_Bases", ""),
                "Total_Sequences": snp.get("Total_Sequences", 0),
                "Valid_Genotypes_Count": snp.get("Valid_Bases_Count", 0),
                "Missing_Genotypes_Count": snp.get("Missing_Genotypes_Count", 0),
                "Call_Rate": snp.get("Call_Rate", 0.0),
                "Percentage_Valid": snp.get("Call_Rate", 0.0) * 100.0,
                "Allele_Counts": str(snp.get("Allele_Counts", {})),
                "Genotype_Counts": str(snp.get("Genotype_Counts", {})),
                "Alignment_Symbol_Counts": str(snp.get("Base_Counts", {})),
                "Minor_Allele_Frequency": snp.get("Minor_Allele_Frequency", 0.0),
                "Heterozygote_Rate": snp.get("Heterozygote_Rate", 0.0),
                "Pass_Call_Rate": snp.get("Pass_Call_Rate", False),
                "Pass_MAF": snp.get("Pass_MAF", False),
                "Pass_Min_Heterozygosity": snp.get(
                    "Pass_Min_Heterozygosity", False
                ),
                "Retained": snp.get("Retained", False),
            }
        )

    dataframe = pd.DataFrame(rows, columns=columns)
    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        dataframe.to_excel(writer, sheet_name="SNP_Analysis_Results", index=False)

    retained_count = int(dataframe["Retained"].sum()) if not dataframe.empty else 0
    print(f"SNP analysis report saved to: {output_excel}")
    print(f"Candidate SNPs: {len(dataframe)}; retained SNPs: {retained_count}")




def _align_two_sequences(seq1_id, seq1, seq2_id, seq2, muscle_exe):
    """Align two sequences using MUSCLE, return alignment object."""
    import tempfile
    tmp_input = None
    tmp_output = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
            tmp_in.write(f">{seq1_id}\n{seq1}\n")
            tmp_in.write(f">{seq2_id}\n{seq2}\n")
            tmp_input = tmp_in.name

        fd, tmp_output = tempfile.mkstemp(suffix='.aln')
        os.close(fd)

        major, _ = _detect_muscle_version(muscle_exe)
        if major >= 5:
            cmd = [muscle_exe, "-align", tmp_input, "-output", tmp_output]
        else:
            cmd = [muscle_exe, "-in", tmp_input, "-out", tmp_output]
        subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=True)

        return AlignIO.read(tmp_output, "fasta")
    except Exception as e:
        if tmp_output and os.path.exists(tmp_output):
            try:
                os.remove(tmp_output)
            except Exception:
                pass
        return None
    finally:
        if tmp_input and os.path.exists(tmp_input):
            try:
                os.remove(tmp_input)
            except Exception:
                pass


def _score_alignment(alignment, ref_id, consensus_id):
    """Evaluate alignment quality, return (matches, identity_ratio)."""
    try:
        ref_idx = cons_idx = None
        for i, rec in enumerate(alignment):
            if rec.id == ref_id:
                ref_idx = i
            elif rec.id == consensus_id:
                cons_idx = i
        if ref_idx is None or cons_idx is None:
            return 0, 0.0

        ref_seq = str(alignment[ref_idx].seq)
        cons_seq = str(alignment[cons_idx].seq)

        matches = sum(1 for r, c in zip(ref_seq, cons_seq)
                      if r == c and r not in '-')
        aligned_cols = sum(1 for r, c in zip(ref_seq, cons_seq)
                           if r not in '-' and c not in '-')
        identity = matches / aligned_cols if aligned_cols > 0 else 0.0

        return matches, identity
    except Exception:
        return 0, 0.0


def _detect_blastn_version(blastn_exe: str) -> Tuple[int, str]:
    """[FIX v1.1] Detect BLASTN version. Returns (major_version, version_string)."""
    try:
        result = subprocess.run([blastn_exe, "-version"], capture_output=True, text=True, timeout=10)
        output = result.stdout + result.stderr
        match = re.search(r'(\d+)\.(\d+)', output)
        if match:
            major = int(match.group(1))
            version_str = match.group(0)
        else:
            major = 2
            version_str = "unknown"
        return major, version_str
    except Exception:
        return 2, "unknown (defaulting to v2 syntax)"


def _align_with_blastn(consensus_file, ref_genome_file, blastn_exe):
    """Run BLASTN to align consensus against reference genome.

    Uses BLASTN -subject mode (no makeblastdb needed), searches all sequences
    in the reference genome file in a single command.
    """
    import tempfile

    # [FIX v1.1] Detect BLASTN version for compatibility
    blast_major, blast_ver = _detect_blastn_version(blastn_exe)
    print(f"BLASTN version detected: {blast_ver} (major={blast_major})")

    tmp_output = None
    try:
        fd, tmp_output = tempfile.mkstemp(suffix='.tsv')
        os.close(fd)

        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
        cmd = [blastn_exe, "-query", consensus_file, "-subject", ref_genome_file,
               "-outfmt", outfmt, "-out", tmp_output, "-evalue", "1e-10",
               "-dust", "no", "-strand", "both"]  # [FIX v1.1] Use "both" strands

        print(f"Running BLASTN: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            print(f"BLASTN failed (returncode={result.returncode}): {result.stderr}")
            return []

        if not os.path.exists(tmp_output) or os.path.getsize(tmp_output) == 0:
            print("BLASTN produced no output (no hits found)")
            return []

        hsps = []
        with open(tmp_output, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 14:
                    continue
                hsp = {
                    'qseqid': fields[0], 'sseqid': fields[1],
                    'pident': float(fields[2]), 'length': int(fields[3]),
                    'mismatch': int(fields[4]), 'gapopen': int(fields[5]),
                    'qstart': int(fields[6]), 'qend': int(fields[7]),
                    'sstart': int(fields[8]), 'send': int(fields[9]),
                    'evalue': float(fields[10]), 'bitscore': float(fields[11]),
                    'qseq': fields[12], 'sseq': fields[13]
                }
                hsps.append(hsp)

        print(f"BLASTN found {len(hsps)} HSP(s)")
        return hsps

    except subprocess.TimeoutExpired:
        print("BLASTN timed out (300s)")
        return []
    except FileNotFoundError:
        print(f"BLASTN executable not found: {blastn_exe}")
        return []
    except Exception as e:
        print(f"BLASTN error: {e}")
        import traceback
        traceback.print_exc()
        return []
    finally:
        if tmp_output and os.path.exists(tmp_output):
            try:
                os.remove(tmp_output)
            except Exception:
                pass


def _build_mapping_from_blastn(
    consensus_id,
    consensus_seq,
    sseqid,
    sstart,
    send,
    qseq,
    sseq,
    qstart=1,
    qend=None,
    mapping_status="unique",
    pident=None,
    alignment_length=None,
    evalue=None,
    bitscore=None,
    query_coverage=None,
):
    """Build a coordinate- and strand-aware BLASTN position mapping.

    Mapping keys are coordinates in the original consensus/query sequence.
    Reference bases are always reported in the reference FASTA forward
    orientation. ``allele_complement_required`` indicates whether alleles
    observed in the original sequencing orientation must be complemented for
    REF/ALT nomenclature on the reference-forward strand.
    """
    del consensus_id
    mapping = {}
    if not qseq or not sseq or len(qseq) != len(sseq):
        print("Invalid BLASTN qseq/sseq alignment")
        return {}, None

    query_residues = sum(base != "-" for base in qseq)
    subject_residues = sum(base != "-" for base in sseq)
    if qend is None:
        qend = qstart + max(query_residues - 1, 0)

    query_direction = 1 if qend >= qstart else -1
    subject_direction = 1 if send >= sstart else -1
    query_position = qstart
    subject_position = sstart
    complement = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
    complement_required = query_direction != subject_direction

    for query_symbol, subject_symbol in zip(qseq.upper(), sseq.upper()):
        current_query_position = None
        current_subject_position = None

        if query_symbol != "-":
            current_query_position = query_position
            query_position += query_direction
        if subject_symbol != "-":
            current_subject_position = subject_position
            subject_position += subject_direction

        if current_query_position is None or current_subject_position is None:
            continue

        reference_base = (
            subject_symbol.translate(complement)
            if subject_direction < 0
            else subject_symbol
        )
        original_query_base = (
            query_symbol.translate(complement)
            if query_direction < 0
            else query_symbol
        )
        mapping[current_query_position] = {
            "reference_sequence": sseqid,
            "genome_position": current_subject_position,
            "reference_base": reference_base,
            "strand": "+" if not complement_required else "-",
            "query_direction": query_direction,
            "subject_direction": subject_direction,
            "query_base": original_query_base,
            "mapping_status": mapping_status,
            "allele_complement_required": complement_required,
            "pident": pident,
            "alignment_length": alignment_length,
            "evalue": evalue,
            "bitscore": bitscore,
            "query_coverage": query_coverage,
        }

    print("BLASTN genome mapping built:")
    print(f"  Reference: {sseqid}")
    print(f"  Query range: {qstart}-{qend}")
    print(f"  Subject range: {sstart}-{send}")
    print(f"  Mapped query positions: {len(mapping)}")
    print(f"  Mapping status: {mapping_status}")
    print(f"  Complement observed alleles: {complement_required}")
    if query_residues == 0 or subject_residues == 0:
        print("  WARNING: BLAST HSP contains no alignable residues")
    return mapping, sseqid



def align_consensus_to_reference(
    consensus_file,
    ref_genome_file,
    output_aln,
    blastn_exe,
    minimum_identity=80.0,
    minimum_query_coverage=0.5,
):
    """Map a consensus to a reference and classify mapping confidence."""
    try:
        if not os.path.exists(consensus_file):
            print(f"Consensus file not found: {consensus_file}")
            return None, None, None
        if not os.path.exists(ref_genome_file):
            print(f"Reference genome file not found: {ref_genome_file}")
            return None, None, None

        records = list(SeqIO.parse(consensus_file, "fasta"))
        if not records:
            print("Consensus file is empty or invalid")
            return None, None, None
        consensus_record = records[0]
        consensus_sequence = str(consensus_record.seq)
        consensus_length = len(consensus_sequence)

        hsps = _align_with_blastn(consensus_file, ref_genome_file, blastn_exe)
        if not hsps:
            print("No BLASTN hits found")
            return None, None, None

        ranked_hsps = sorted(
            hsps,
            key=lambda hit: (
                -hit["bitscore"],
                hit["evalue"],
                -hit["pident"],
                -hit["length"],
                hit["sseqid"],
            ),
        )
        best = ranked_hsps[0]
        aligned_query_bases = sum(base != "-" for base in best["qseq"])
        query_coverage = (
            aligned_query_bases / consensus_length if consensus_length else 0.0
        )

        near_ties = [
            hit
            for hit in ranked_hsps
            if hit["bitscore"] >= best["bitscore"] * 0.99
            and hit["evalue"] <= max(best["evalue"] * 10.0, 1e-300)
        ]
        best_subject_interval = sorted((best["sstart"], best["send"]))

        def is_distinct_locus(hit):
            if hit["sseqid"] != best["sseqid"]:
                return True
            interval = sorted((hit["sstart"], hit["send"]))
            overlap = max(
                0,
                min(best_subject_interval[1], interval[1])
                - max(best_subject_interval[0], interval[0])
                + 1,
            )
            shorter = min(
                best_subject_interval[1] - best_subject_interval[0] + 1,
                interval[1] - interval[0] + 1,
            )
            return shorter > 0 and overlap / shorter < 0.5

        ambiguous_hits = [hit for hit in near_ties[1:] if is_distinct_locus(hit)]
        if ambiguous_hits:
            mapping_status = "ambiguous"
        elif best["pident"] < minimum_identity or query_coverage < minimum_query_coverage:
            mapping_status = "low_confidence"
        else:
            mapping_status = "unique"

        print("Best BLASTN hit:")
        print(f"  Subject: {best['sseqid']}")
        print(f"  Identity: {best['pident']:.2f}%")
        print(f"  Query coverage: {query_coverage:.4f}")
        print(f"  Query range: {best['qstart']}-{best['qend']}")
        print(f"  Subject range: {best['sstart']}-{best['send']}")
        print(f"  Mapping status: {mapping_status}")

        mapping, reference_name = _build_mapping_from_blastn(
            consensus_record.id,
            consensus_sequence,
            best["sseqid"],
            best["sstart"],
            best["send"],
            best["qseq"],
            best["sseq"],
            qstart=best["qstart"],
            qend=best["qend"],
            mapping_status=mapping_status,
            pident=best["pident"],
            alignment_length=best["length"],
            evalue=best["evalue"],
            bitscore=best["bitscore"],
            query_coverage=query_coverage,
        )

        with open(output_aln, "w", encoding="utf-8") as output_handle:
            output_handle.write("# BLASTN alignment result\n")
            output_handle.write(f"# Software: ABI-SNP Analyzer v{SOFTWARE_VERSION}\n")
            output_handle.write(f"# Query: {consensus_record.id}\n")
            output_handle.write(f"# Subject: {best['sseqid']}\n")
            output_handle.write(f"# Identity: {best['pident']:.6f}%\n")
            output_handle.write(f"# Query coverage: {query_coverage:.6f}\n")
            output_handle.write(f"# Query range: {best['qstart']}-{best['qend']}\n")
            output_handle.write(f"# Subject range: {best['sstart']}-{best['send']}\n")
            output_handle.write(f"# Mapping status: {mapping_status}\n")
            output_handle.write(f"# Bit score: {best['bitscore']}\n")
            output_handle.write(f"# E-value: {best['evalue']}\n")
            if ambiguous_hits:
                output_handle.write("# Ambiguous near-tied candidate loci:\n")
                for candidate in ambiguous_hits:
                    output_handle.write(
                        f"#   {candidate['sseqid']}:"
                        f"{candidate['sstart']}-{candidate['send']} "
                        f"bitscore={candidate['bitscore']} "
                        f"identity={candidate['pident']}\n"
                    )
            output_handle.write(f">Query_aln\n{best['qseq']}\n")
            output_handle.write(f">Subject_aln\n{best['sseq']}\n")
        print(f"BLASTN alignment information saved to: {output_aln}")
        return None, mapping, reference_name
    except Exception as exc:
        print(f"BLASTN reference mapping error: {exc}")
        import traceback

        traceback.print_exc()
        return None, None, None



def build_genome_mapping(alignment, ref_id, consensus_id, ref_seq):
    """Build consensus-to-genome position mapping from alignment (legacy, pre-v1.1).

    Deprecated: use _build_mapping_from_blastn() instead which provides
    structured mapping with REF base from the reference genome.
    """
    mapping = {}

    try:
        ref_idx = None
        cons_idx = None
        for i, record in enumerate(alignment):
            if record.id == ref_id:
                ref_idx = i
            elif record.id == consensus_id:
                cons_idx = i

        if ref_idx is None or cons_idx is None:
            print("Reference genome or consensus sequence not found in alignment")
            return mapping

        ref_aln_seq = str(alignment[ref_idx].seq)
        cons_aln_seq = str(alignment[cons_idx].seq)

        ref_start_in_aln = 0
        for i, base in enumerate(ref_aln_seq):
            if base != '-':
                ref_start_in_aln = i
                break

        consensus_pos = 0
        ref_genome_pos = ref_start_in_aln + 1

        for i in range(len(ref_aln_seq)):
            ref_base = ref_aln_seq[i]
            cons_base = cons_aln_seq[i]

            if cons_base != '-':
                consensus_pos += 1

            if ref_base != '-':
                if cons_base != '-':
                    mapping[consensus_pos] = ref_genome_pos
                ref_genome_pos += 1

        print(f"\nSNP genome position mapping built:")
        print(f"  Consensus length: {consensus_pos}")
        print(f"  Mapped positions: {len(mapping)}")

        if mapping:
            sample_keys = list(mapping.keys())[:5]
            print(f"  Sample mappings (first 5): {[(k, mapping[k]) for k in sample_keys]}")

    except Exception as e:
        print(f"Error building genome position mapping: {e}")
        import traceback
        traceback.print_exc()

    return mapping, ref_id


def get_genome_snp_name(genome_pos, alleles, ref_base, ref_seq_name=None):
    """Return a reference-forward SNP name, preserving every ALT allele."""
    if genome_pos in {None, ""}:
        return ""
    reference = str(ref_base).strip().upper()
    if reference not in "ACGT" or len(reference) != 1:
        return ""
    observed = sorted(
        {str(allele).strip().upper() for allele in alleles if str(allele).strip().upper() in "ACGT"},
        key="ACGT".index,
    )
    alternates = [allele for allele in observed if allele != reference]
    if not alternates:
        return ""
    suffix = f"{reference}>{','.join(alternates)}"
    if ref_seq_name:
        clean_name = str(ref_seq_name).split()[0].replace(" ", "_")
        return f"{clean_name}:{genome_pos}{suffix}"
    return f"{genome_pos}{suffix}"



def compute_consensus_position_mapping(alignment):
    """Map alignment column indices (0-based) to consensus sequence positions (1-based).

    Only alignment columns that contain at least one non-gap base count toward
    the consensus.  Returns: {alignment_col_index: consensus_1based_position}.
    """
    mapping = {}
    consensus_pos = 0
    for col in range(alignment.get_alignment_length()):
        has_base = any(
            str(rec.seq[col]).upper() not in ("-", ".")
            for rec in alignment
        )
        if has_base:
            consensus_pos += 1
            mapping[col] = consensus_pos
    return mapping


def compute_sample_position_maps(alignment, consensus_pos_map):
    """For each sample, map 1-based sample sequence position -> 1-based consensus position.

    Gaps introduced by MUSCLE are skipped.  Returns:
        {sample_id: {sample_1based_pos: consensus_1based_pos}}
    """
    sample_maps = {}
    for rec in alignment:
        sample_seq = str(rec.seq).upper()
        sample_pos = 0
        pos_map = {}
        for col, char in enumerate(sample_seq):
            if char not in ("-", "."):
                sample_pos += 1
                cons_pos = consensus_pos_map.get(col)
                if cons_pos is not None:
                    pos_map[sample_pos] = cons_pos
        sample_maps[rec.id] = pos_map
    return sample_maps


def create_genotype_matrix(
    alignment_file,
    output_excel,
    similarity_threshold,
    minor_allele_freq_threshold=0.10,
    min_observed_heterozygosity=0.1,
    genome_snp_mapping=None,
    ref_seq_name=None,
    minimum_call_rate=0.5,
):
    """Create a diploid genotype matrix from reliable decoded alleles.

    Candidate SNPs are alignment columns with at least two observed A/C/G/T
    alleles after decoding reliable diploid IUPAC symbols. V/H/D/B/N/gaps are
    missing calls. Filtering uses call rate, MAF, and the optional minimum
    observed heterozygosity (default disabled).
    """
    try:
        similar_sequences, excluded_sequences = calculate_sequence_similarity(
            alignment_file, similarity_threshold
        )
        del excluded_sequences

        if not similar_sequences:
            print("No sequences meet the similarity threshold")
            empty_dataframe = pd.DataFrame(columns=["Sample"])
            with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
                empty_dataframe.to_excel(
                    writer, sheet_name="Genotype_Matrix", index=False
                )
            return

        alignment = AlignIO.read(alignment_file, "fasta")
        retained_ids = {item["id"] for item in similar_sequences}
        retained_samples = sorted(retained_ids)
        snp_data = []
        snp_genotype_matrix = {}

        # Compute consensus position mapping: alignment column -> consensus 1-based position
        alignment_to_consensus = compute_consensus_position_mapping(alignment)

        for position in range(alignment.get_alignment_length()):
            column = _analyze_snp_column(alignment, retained_ids, position)
            observed_alleles = column["observed_alleles"]
            if len(observed_alleles) < 2:
                continue

            consensus_position = alignment_to_consensus.get(position, position + 1)
            reference_base = ""
            genome_position_name = ""
            reference_alleles = list(observed_alleles)
            mapping_status = "unmapped"
            mapping_strand = ""
            mapping = (
                genome_snp_mapping.get(consensus_position)
                if genome_snp_mapping
                else None
            )

            if isinstance(mapping, dict):
                mapping_status = mapping.get("mapping_status", "low_confidence")
                mapping_strand = mapping.get("strand", "")
                reference_base = str(mapping.get("reference_base", "")).upper()
                if mapping.get("allele_complement_required"):
                    complement = str.maketrans("ACGT", "TGCA")
                    reference_alleles = [
                        allele.translate(complement) for allele in observed_alleles
                    ]
                    reference_alleles.sort(key="ACGT".index)
                genome_position = mapping.get("genome_position", "")
                reference_name = mapping.get("reference_sequence", ref_seq_name)
                if mapping_status == "unique" and genome_position:
                    genome_position_name = get_genome_snp_name(
                        genome_position,
                        reference_alleles,
                        reference_base,
                        reference_name,
                    )

            snp_id = genome_position_name or f"SNP_{consensus_position}"
            if column["allele_counts"]:
                consensus_base = min(
                    column["allele_counts"],
                    key=lambda allele: (
                        -column["allele_counts"][allele],
                        "ACGT".index(allele),
                    ),
                )
            else:
                consensus_base = "N"

            pass_call_rate = column["call_rate"] >= minimum_call_rate
            pass_maf = column["maf"] >= minor_allele_freq_threshold
            pass_heterozygosity = (
                column["heterozygote_rate"] >= min_observed_heterozygosity
            )
            retained = pass_call_rate and pass_maf and pass_heterozygosity

            snp_genotype_matrix[snp_id] = column["sample_genotypes"]
            snp_data.append(
                {
                    "Position": consensus_position,
                    "SNP_ID": snp_id,
                    "Genome_Position": genome_position_name,
                    "Consensus_Base": consensus_base,
                    "Reference_Base": reference_base,
                    "Sequencing_Alleles": "/".join(observed_alleles),
                    "Reference_Alleles": "/".join(reference_alleles),
                    "Mapping_Status": mapping_status,
                    "Mapping_Strand": mapping_strand,
                    "Mapping_Identity": mapping.get("pident", "") if isinstance(mapping, dict) else "",
                    "Mapping_Query_Coverage": mapping.get("query_coverage", "") if isinstance(mapping, dict) else "",
                    "Mapping_Evalue": mapping.get("evalue", "") if isinstance(mapping, dict) else "",
                    "Mapping_Bitscore": mapping.get("bitscore", "") if isinstance(mapping, dict) else "",
                    "All_Bases": "/".join(reference_alleles),
                    "Total_Sequences": column["total_sequences"],
                    "Valid_Bases_Count": column["called_genotypes"],
                    "Missing_Genotypes_Count": column["missing_genotypes"],
                    "Call_Rate": column["call_rate"],
                    "Base_Counts": dict(column["symbol_counts"]),
                    "Allele_Counts": dict(column["allele_counts"]),
                    "Genotype_Counts": dict(column["genotype_counts"]),
                    "Minor_Allele_Frequency": column["maf"],
                    "Heterozygote_Rate": column["heterozygote_rate"],
                    "Pass_Call_Rate": pass_call_rate,
                    "Pass_MAF": pass_maf,
                    "Pass_Min_Heterozygosity": pass_heterozygosity,
                    "Retained": retained,
                }
            )

        retained_snps = [item["SNP_ID"] for item in snp_data if item["Retained"]]
        position_by_id = {item["SNP_ID"]: item["Position"] for item in snp_data}
        retained_snps.sort(key=lambda snp_id: position_by_id[snp_id])

        genotype_rows = []
        for sample in retained_samples:
            row = {"Sample": sample}
            for snp_id in retained_snps:
                row[snp_id] = snp_genotype_matrix[snp_id].get(sample, "--")
            genotype_rows.append(row)

        genotype_dataframe = pd.DataFrame(
            genotype_rows, columns=["Sample"] + retained_snps
        )
        with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
            genotype_dataframe.to_excel(
                writer, sheet_name="Genotype_Matrix", index=False
            )

        report_output = output_excel.replace(".xlsx", "_snp_analysis.xlsx")
        create_snp_analysis_report(snp_data, report_output)

        print(f"Genotype matrix saved to: {output_excel}")
        print(
            f"Matrix dimensions: {len(retained_samples)} samples x "
            f"{len(retained_snps)} SNPs"
        )
        print(
            f"Candidate SNPs: {len(snp_data)}; call-rate threshold: "
            f"{minimum_call_rate}; MAF threshold: {minor_allele_freq_threshold}; "
            f"minimum observed heterozygosity: {min_observed_heterozygosity}"
        )
    except Exception as exc:
        print(f"Error generating genotype matrix: {exc}")
        import traceback

        traceback.print_exc()
        raise






def _parse_population_genotype(value):
    """Parse a genotype-matrix value into a complete diploid call or None."""
    if value is None:
        return None
    text = str(value).strip().upper()
    if text in {"", "--", "NN", "N", "N/A", "NA", "NAN"}:
        return None
    if len(text) == 1 and text in "ACGT":
        return text, text
    if len(text) >= 2 and text[0] in "ACGT" and text[1] in "ACGT":
        return text[0], text[1]
    return None


def build_population_analysis_inputs(genotype_matrix_data, snp_names, populations):
    """Build consistent locus-level and sample-level statistical inputs.

    Returns ``(all_data, sample_genotypes, sample_to_population, unassigned)``.
    A sample assigned to more than one population raises ``ValueError``.
    """
    sample_to_population = {}
    for population, samples in populations.items():
        for sample in samples:
            previous = sample_to_population.get(sample)
            if previous is not None and previous != population:
                raise ValueError(
                    f"Sample {sample!r} is assigned to both {previous!r} and "
                    f"{population!r}"
                )
            sample_to_population[sample] = population

    all_known_samples = {
        sample
        for locus_calls in genotype_matrix_data.values()
        for sample in locus_calls
    }
    unassigned = sorted(all_known_samples - set(sample_to_population), key=str)

    sample_genotypes = {}
    for sample in sorted(sample_to_population, key=str):
        calls = {}
        for locus in snp_names:
            value = genotype_matrix_data.get(locus, {}).get(sample)
            calls[locus] = _parse_population_genotype(value)
        sample_genotypes[sample] = calls

    all_data = {}
    for locus in snp_names:
        by_population = {}
        for sample, population in sample_to_population.items():
            genotype = sample_genotypes.get(sample, {}).get(locus)
            if genotype is not None:
                by_population.setdefault(population, []).append(genotype)
        if len(by_population) >= 2:
            all_data[locus] = by_population

    return all_data, sample_genotypes, sample_to_population, unassigned



def _sha256_file(path, chunk_size=1024 * 1024):
    """Return a streaming SHA-256 checksum without loading a file into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_metadata(path):
    """Return privacy-preserving file metadata (basename only)."""
    if not path:
        return None
    result = {"file_name": os.path.basename(path)}
    try:
        result["size_bytes"] = os.path.getsize(path)
        result["sha256"] = _sha256_file(path)
    except OSError as exc:
        result["metadata_error"] = str(exc)
    return result


def _package_versions():
    packages = ["numpy", "pandas", "scipy", "matplotlib", "biopython", "openpyxl"]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _external_tool_metadata(executable, version_arguments=("-version",)):
    if not executable:
        return None
    metadata = _file_metadata(executable) or {"file_name": os.path.basename(executable)}
    try:
        result = subprocess.run(
            [executable, *version_arguments],
            capture_output=True,
            text=True,
            timeout=10,
        )
        version_output = (result.stdout + result.stderr).strip().splitlines()
        metadata["version"] = version_output[0] if version_output else "unknown"
        metadata["version_return_code"] = result.returncode
    except (OSError, subprocess.TimeoutExpired) as exc:
        metadata["version"] = "unavailable"
        metadata["version_error"] = str(exc)
    return metadata


def _collect_output_metadata(output_dir, excluded_names=None):
    excluded = set(excluded_names or ())
    outputs = []
    if not os.path.isdir(output_dir):
        return outputs
    for root_dir, _, file_names in os.walk(output_dir):
        for file_name in sorted(file_names):
            if file_name in excluded:
                continue
            path = os.path.join(root_dir, file_name)
            relative = os.path.relpath(path, output_dir).replace(os.sep, "/")
            try:
                outputs.append(
                    {
                        "relative_path": relative,
                        "size_bytes": os.path.getsize(path),
                        "sha256": _sha256_file(path),
                    }
                )
            except OSError as exc:
                outputs.append({"relative_path": relative, "metadata_error": str(exc)})
    return outputs


def write_analysis_metadata(output_dir, metadata, file_name="analysis_metadata.json"):
    """Finalize and atomically write a reproducibility manifest."""
    os.makedirs(output_dir, exist_ok=True)
    finalized = dict(metadata)
    finalized["generated_outputs"] = _collect_output_metadata(
        output_dir, excluded_names={file_name}
    )
    destination = os.path.join(output_dir, file_name)
    temporary = destination + ".tmp"
    with open(temporary, "w", encoding="utf-8") as file_handle:
        json.dump(finalized, file_handle, indent=2, ensure_ascii=False, sort_keys=True)
        file_handle.write("\n")
    os.replace(temporary, destination)
    return destination


# ========== GUI section ==========
class ABIToSNPApp:
    def __init__(self, root):
        self.root = root
        self.root.title(f"ABI-SNP Analyzer v{SOFTWARE_VERSION} (MIT License © 2026 Xia Chengyu et al. Xidian University)")
        self.root.geometry("1000x700")
        self.root.minsize(900, 600)
        self.root.resizable(True, True)
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        # Initialize variables
        self.abi_files = []
        self.output_dir = ""
        self.muscle_path = find_external_tool("muscle.exe") or os.path.join(APP_DIR, "muscle.exe")

        # Create main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)

        # ========== 1. File selection area ==========
        file_frame = ttk.LabelFrame(main_frame, text="File Selection", padding="10")
        file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        file_frame.columnconfigure(1, weight=1)

        ttk.Label(file_frame, text="ABI Files:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.abi_file_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.abi_file_var, state="readonly").grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_frame, text="Select Files/Folder", command=self.select_abi_files).grid(row=0, column=2, padx=5)

        ttk.Label(file_frame, text="MUSCLE Path:").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.muscle_path_var = tk.StringVar(value=self.muscle_path)
        ttk.Entry(file_frame, textvariable=self.muscle_path_var).grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_frame, text="Browse", command=self.select_muscle_path).grid(row=1, column=2, padx=5)

        ttk.Label(file_frame, text="Reference Genome (optional):").grid(row=2, column=0, sticky=tk.W, padx=5)
        self.ref_genome_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.ref_genome_var, width=50).grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_frame, text="Browse", command=self.select_ref_genome).grid(row=2, column=2, padx=5)

        ttk.Label(file_frame, text="BLASTN Path (optional):").grid(row=3, column=0, sticky=tk.W, padx=5)
        self.blastn_path_var = tk.StringVar(value=find_external_tool("blastn.exe"))
        ttk.Entry(file_frame, textvariable=self.blastn_path_var, width=50).grid(row=3, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_frame, text="Browse", command=self.select_blastn_path).grid(row=3, column=2, padx=5)

        ttk.Label(file_frame, text="Output Directory:").grid(row=4, column=0, sticky=tk.W, padx=5)
        self.output_dir_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.output_dir_var, state="readonly").grid(row=4, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(file_frame, text="Select", command=self.select_output_dir).grid(row=4, column=2, padx=5)

        # ========== 2. Filter threshold settings ==========
        param_frame = ttk.LabelFrame(main_frame, text="Filter Threshold Settings", padding="10")
        param_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        ttk.Label(param_frame, text="Sequence Similarity:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.similarity_threshold_var = tk.DoubleVar(value=0.9)
        ttk.Entry(param_frame, textvariable=self.similarity_threshold_var, width=10).grid(row=0, column=1, padx=5)
        ttk.Label(param_frame, text="(0-1, default 0.9)").grid(row=0, column=2, sticky=tk.W)

        ttk.Label(param_frame, text="MAF Threshold:").grid(row=0, column=3, sticky=tk.W, padx=5)
        self.maf_threshold_var = tk.DoubleVar(value=0.10)
        ttk.Entry(param_frame, textvariable=self.maf_threshold_var, width=10).grid(row=0, column=4, padx=5)
        ttk.Label(param_frame, text="(0-1, default 0.10)").grid(row=0, column=5, sticky=tk.W)

        # [FIX v1.1] Renamed: "Minimum observed heterozygosity (optional)" with default 0.1
        ttk.Label(param_frame, text="Min Observed Het (optional):").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.min_obs_het_var = tk.DoubleVar(value=0.1)  # [FIX v1.1] Default 0.1
        ttk.Entry(param_frame, textvariable=self.min_obs_het_var, width=10).grid(row=1, column=1, padx=5)
        ttk.Label(param_frame, text="(0-1, default 0.0)").grid(row=1, column=2, sticky=tk.W)

        # [FIX v1.1] Clarified signal ratio definition
        ttk.Label(param_frame, text="Signal Ratio (sec./primary):").grid(row=1, column=3, sticky=tk.W, padx=5)
        self.sigratio_var = tk.DoubleVar(value=0.3)
        ttk.Entry(param_frame, textvariable=self.sigratio_var, width=10).grid(row=1, column=4, padx=5)
        ttk.Label(param_frame, text="(0-1, default 0.3)").grid(row=1, column=5, sticky=tk.W)

        ttk.Label(param_frame, text="Consensus Quality:").grid(row=2, column=0, sticky=tk.W, padx=5)
        self.quality_threshold_var = tk.DoubleVar(value=0.2)
        ttk.Entry(param_frame, textvariable=self.quality_threshold_var, width=10).grid(row=2, column=1, padx=5)
        ttk.Label(param_frame, text="(0-1, default 0.2)").grid(row=2, column=2, sticky=tk.W)

        ttk.Label(param_frame, text="Min Base Quality (optional):").grid(row=2, column=3, sticky=tk.W, padx=5)
        self.min_base_quality_var = tk.IntVar(value=5)
        ttk.Entry(param_frame, textvariable=self.min_base_quality_var, width=10).grid(row=2, column=4, padx=5)
        ttk.Label(param_frame, text="(0-100, default 10)").grid(row=2, column=5, sticky=tk.W)

        # ========== 3. Run buttons ==========
        run_frame = ttk.Frame(main_frame)
        run_frame.grid(row=2, column=0, columnspan=3, pady=10)

        self.run_button = ttk.Button(run_frame, text="Start Analysis", command=self.run_analysis)
        self.run_button.grid(row=0, column=0, padx=5)

        self.stop_button = ttk.Button(run_frame, text="Stop", command=self.stop_analysis, state=tk.DISABLED)
        self.stop_button.grid(row=0, column=1, padx=5)

        self.load_data_button = ttk.Button(run_frame, text="Load Existing Data",
                                           command=self._load_genotype_data)
        self.load_data_button.grid(row=0, column=2, padx=5)

        self.pop_manager_button = ttk.Button(run_frame, text="Population Manager",
                                             command=self._open_population_manager, state=tk.DISABLED)
        self.pop_manager_button.grid(row=0, column=3, padx=5)

        self.stats_button = ttk.Button(run_frame, text="Statistical Analysis",
                                       command=self._run_statistical_analysis, state=tk.DISABLED)
        self.stats_button.grid(row=0, column=4, padx=5)

        self.export_formats_button = ttk.Button(run_frame, text="Export SNP Formats",
                                                command=self._export_snp_formats, state=tk.DISABLED)
        self.export_formats_button.grid(row=0, column=5, padx=5)

        # ========== 4. Log output area ==========
        log_frame = ttk.LabelFrame(main_frame, text="Run Log", padding="10")
        log_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, height=20)
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Thread-safe GUI/log queues must exist before stdout is redirected.
        self._gui_queue = queue.Queue()
        self._log_queue = queue.Queue()
        self.original_stdout = sys.stdout
        sys.stdout = self

        # Thread control
        self.analysis_thread = None
        self.stop_flag = False

        # Population genetics state
        self.genotype_matrix_data = None
        self.sample_names = []
        self.snp_names = []
        self.populations = {}
        self.last_output_dir = ""

        self._process_gui_queue()

    def _process_gui_queue(self):
        """Process worker callbacks and redirected stdout on the Tk main thread."""
        try:
            while True:
                text = self._log_queue.get_nowait()
                self.log_text.insert(tk.END, text)
                self.log_text.see(tk.END)
        except queue.Empty:
            pass

        try:
            while True:
                callback = self._gui_queue.get_nowait()
                try:
                    callback()
                except Exception as exc:
                    self.log_text.insert(tk.END, f"[GUI Queue] {exc}\n")
                    self.log_text.see(tk.END)
        except queue.Empty:
            pass

        try:
            self.root.after(100, self._process_gui_queue)
        except tk.TclError:
            return

    def _safe_gui_call(self, callback):
        """Schedule a callback to run on the Tk main thread."""
        self._gui_queue.put(callback)

    def write(self, text):
        """Queue redirected stdout without touching Tk from worker threads."""
        if text:
            self._log_queue.put(str(text))

    def flush(self):
        """flush implementation"""
        pass

    def _on_close(self):
        """Stop workers, restore stdout, and close the Tk application safely."""
        if self.analysis_thread is not None and self.analysis_thread.is_alive():
            should_close = messagebox.askyesno(
                "Analysis Running",
                "An analysis is still running. Stop it and close the application?",
                parent=self.root,
            )
            if not should_close:
                return
            self.stop_flag = True
        if sys.stdout is self:
            sys.stdout = self.original_stdout
        try:
            self.root.destroy()
        except tk.TclError:
            pass

    def select_abi_files(self):
        choice = messagebox.askyesno("Selection Type", "Select a folder?\n(Yes=folder, No=single/multiple files)")
        if choice:
            dir_path = filedialog.askdirectory(title="Select folder containing ABI files")
            if dir_path:
                self.abi_files = list_ab1_files(dir_path)
                self.abi_file_var.set(f"Folder: {dir_path} ({len(self.abi_files)} AB1 files)")
        else:
            file_paths = filedialog.askopenfilenames(
                title="Select ABI files",
                filetypes=[("AB1 files", "*.ab1"), ("All files", "*.*")]
            )
            if file_paths:
                self.abi_files = list(file_paths)
                self.abi_file_var.set(f"Files: {len(self.abi_files)}")

    def select_muscle_path(self):
        muscle_exe = filedialog.askopenfilename(
            title="Select MUSCLE executable",
            filetypes=[("Executable files", "*.exe"), ("All files", "*.*")],
            initialdir=os.path.dirname(self.muscle_path)
        )
        if muscle_exe:
            self.muscle_path = muscle_exe
            self.muscle_path_var.set(muscle_exe)

    def select_ref_genome(self):
        ref_genome = filedialog.askopenfilename(
            title="Select reference genome file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
        )
        if ref_genome:
            self.ref_genome_var.set(ref_genome)

    def select_blastn_path(self):
        blastn_exe = filedialog.askopenfilename(
            title="Select BLASTN executable (blastn.exe)",
            filetypes=[("Executable files", "*.exe"), ("All files", "*.*")]
        )
        if blastn_exe:
            self.blastn_path_var.set(blastn_exe)

    def select_output_dir(self):
        dir_path = filedialog.askdirectory(title="Select output directory")
        if dir_path:
            self.output_dir = dir_path
            self.output_dir_var.set(dir_path)

    def run_analysis(self):
        if not self.abi_files:
            messagebox.showerror("Error", "Please select ABI files first!")
            return
        if not self.output_dir:
            messagebox.showerror("Error", "Please select output directory first!")
            return
        if not os.path.exists(self.muscle_path):
            messagebox.showerror("Error", f"MUSCLE file not found: {self.muscle_path}")
            return

        ref_genome = self.ref_genome_var.get().strip()
        if ref_genome and not os.path.exists(ref_genome):
            messagebox.showerror("Error", f"Reference genome file not found: {ref_genome}")
            return

        blastn_path = self.blastn_path_var.get().strip()
        if ref_genome:
            if not blastn_path:
                messagebox.showerror("Error",
                    "BLASTN path is required when a reference genome is provided!\n"
                    "Please specify the path to blastn.exe (e.g., ncbi-blast-x64-win64/bin/blastn.exe)")
                return
            if not os.path.exists(blastn_path):
                messagebox.showerror("Error", f"BLASTN executable not found: {blastn_path}")
                return

        try:
            similarity_threshold = self.similarity_threshold_var.get()
            maf_threshold = self.maf_threshold_var.get()
            min_obs_het = self.min_obs_het_var.get()
            sigratio = self.sigratio_var.get()
            quality_threshold = self.quality_threshold_var.get()
            minimum_base_quality = self.min_base_quality_var.get()

            for name, val in [
                ("Similarity threshold", similarity_threshold),
                ("MAF threshold", maf_threshold),
                ("Minimum observed heterozygosity", min_obs_het),
                ("Signal ratio threshold", sigratio),
                ("Quality threshold", quality_threshold)
            ]:
                if not (0.0 <= val <= 1.0):
                    messagebox.showerror("Error", f"{name} must be between 0 and 1!")
                    return
            if not (0 <= minimum_base_quality <= 100):
                messagebox.showerror("Error", "Minimum base quality must be between 0 and 100!")
                return
        except Exception as exc:
            error_msg = str(exc)
            messagebox.showerror("Error", f"Parameter input error: {error_msg}")
            return

        self.run_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)
        self.stop_flag = False

        self.analysis_thread = threading.Thread(
            target=self._analysis_worker,
            args=(
                self.abi_files,
                self.output_dir,
                self.muscle_path,
                similarity_threshold,
                maf_threshold,
                min_obs_het,
                sigratio,
                quality_threshold,
                minimum_base_quality,
                ref_genome,
                blastn_path
            ),
            daemon=True
        )
        self.analysis_thread.start()

    def stop_analysis(self):
        self.stop_flag = True
        self.stop_button.config(state=tk.DISABLED)
        print("\nStopping analysis...")

    def _analysis_worker(
        self,
        abi_files,
        output_dir,
        muscle_path,
        similarity_threshold,
        maf_threshold,
        min_obs_het,
        sigratio,
        quality_threshold,
        minimum_base_quality,
        ref_genome=None,
        blastn_path=None,
    ):
        """Run the primary pipeline and always write a reproducibility manifest."""
        started_wall = datetime.now().astimezone()
        started_timer = time.perf_counter()
        genome_snp_mapping = None
        ref_seq_name = None
        run_status = "running"
        failure_message = None
        input_records = []
        for abi_file in abi_files:
            record = _file_metadata(abi_file) or {"file_name": os.path.basename(abi_file)}
            record.update({"status": "pending", "reason": ""})
            input_records.append(record)

        metadata = {
            "software": {"name": "ABI-SNP Analyzer", "version": SOFTWARE_VERSION},
            "run": {
                "started_at": started_wall.isoformat(),
                "status": run_status,
            },
            "runtime": {
                "python_version": platform.python_version(),
                "python_implementation": platform.python_implementation(),
                "operating_system": platform.system(),
                "operating_system_release": platform.release(),
                "machine": platform.machine(),
                "frozen_executable": bool(getattr(sys, "frozen", False)),
            },
            "dependencies": _package_versions(),
            "tools": {
                "muscle": _external_tool_metadata(muscle_path),
                "blastn": _external_tool_metadata(blastn_path) if blastn_path else None,
            },
            "parameters": {
                "sequence_similarity_threshold": similarity_threshold,
                "minimum_overlap_coverage": 0.8,
                "maf_threshold": maf_threshold,
                "minimum_observed_heterozygosity": min_obs_het,
                "signal_ratio_threshold": sigratio,
                "signal_ratio_formula": "secondary_peak_height / primary_peak_height",
                "consensus_valid_sequence_threshold": quality_threshold,
                "minimum_base_quality": minimum_base_quality,
                "trace_end_trim_bases": 20,
                "statistical_permutations_default": 999,
                "statistical_random_seed_default": 42,
            },
            "inputs": input_records,
            "reference_genome": _file_metadata(ref_genome) if ref_genome else None,
            "summary": {
                "input_file_count": len(abi_files),
                "successfully_parsed_count": 0,
                "failed_input_count": 0,
                "retained_genotype_sample_count": 0,
                "output_snp_count": 0,
            },
        }

        try:
            print(f"===== Analysis started - {started_wall.strftime('%Y-%m-%d %H:%M:%S %z')} =====")
            print(f"Software version: {SOFTWARE_VERSION}")
            print(f"ABI files: {len(abi_files)}")
            print("Parameter settings:")
            print(f"  Sequence similarity threshold: {similarity_threshold}")
            print("  Minimum overlap coverage: 0.8")
            print(f"  MAF threshold: {maf_threshold}")
            print(f"  Minimum observed heterozygosity: {min_obs_het}")
            print(f"  Signal ratio (secondary/primary): {sigratio}")
            print(f"  Consensus valid-sequence threshold: {quality_threshold}")
            print(f"  Minimum base quality: {minimum_base_quality} (default 5)")

            parser = AbifParser()
            total_files = len(abi_files)
            failed_files = []
            for idx, (abi_file, record) in enumerate(zip(abi_files, input_records)):
                if self.stop_flag:
                    run_status = "stopped_by_user"
                    record["status"] = "not_processed"
                    record["reason"] = "analysis stopped before this file"
                    continue

                # [FIX v1.2] Batch progress: only print every 50 files or on errors.
                # Per-file print with 400+ samples floods the Windows console buffer
                # and can cause PyInstaller EXE to hang.
                if (idx + 1) % 50 == 0 or idx == 0:
                    print(f"Processing file {idx+1}/{total_files}...")
                    sys.stdout.flush()

                try:
                    trace_data = parser.read_abi(abi_file)
                    if trace_data is None:
                        record["status"] = "failed"
                        record["reason"] = "invalid or incomplete ABIF data"
                        failed_files.append(os.path.basename(abi_file))
                        continue
                    heterozygotes = parser.detect_heterozygotes(
                        trace_data,
                        sigratio=sigratio,
                        minimum_base_quality=minimum_base_quality,
                    )
                    iupac_sequence = parser.generate_iupac_sequence(
                        trace_data, heterozygotes, trim_bp=20
                    )
                    save_fasta_file_with_iupac(
                        trace_data, iupac_sequence, abi_file, output_dir
                    )

                    # [FIX v1.0.0] Defer trace plot and report generation to Phase 2
                    # (after alignment so we can use consensus positions).
                    # Only save FASTA now to reduce memory pressure with many samples.

                    record["_heterozygotes"] = heterozygotes
                    record["_abi_file"] = abi_file

                    record["status"] = "parsed"
                    record["sequence_length"] = len(trace_data.basecalls1)
                    record["quality_score_count"] = len(trace_data.qual)
                    record["heterozygous_call_count"] = sum(
                        item.get("qc_status") == "ok" for item in heterozygotes
                    )
                    record["ambiguous_peak_count"] = sum(
                        item.get("qc_status") == "ambiguous_3plus_peaks"
                        for item in heterozygotes
                    )

                    # [MEMORY] Release large trace_data arrays immediately
                    del trace_data

                    # [MEMORY] Periodic garbage collection every 20 files
                    if (idx + 1) % 20 == 0:
                        gc.collect()

                except Exception as exc:
                    record["status"] = "failed"
                    record["reason"] = f"processing error: {exc}"
                    failed_files.append(f"{os.path.basename(abi_file)}: {exc}")

            # [FIX v1.2] Report failed files in one batch instead of per-file spam
            if failed_files:
                print(f"\n{len(failed_files)} file(s) failed to parse:")
                for f in failed_files[:20]:  # Cap at 20 lines
                    print(f"  - {f}")
                if len(failed_files) > 20:
                    print(f"  ... and {len(failed_files) - 20} more")
                sys.stdout.flush()

            # [MEMORY] Final GC after Phase 1
            gc.collect()
            metadata["summary"]["successfully_parsed_count"] = sum(
                record["status"] == "parsed" for record in input_records
            )
            metadata["summary"]["failed_input_count"] = sum(
                record["status"] == "failed" for record in input_records
            )
            het_sample_count = sum(
                1 for r in input_records
                if r.get("_heterozygotes") and any(
                    h.get("qc_status") == "ok" for h in r["_heterozygotes"]
                )
            )
            print(
                f"\nPhase 1 complete. "
                f"{metadata['summary']['successfully_parsed_count']}/{total_files} "
                f"files parsed, {het_sample_count} with heterozygous calls."
            )
            sys.stdout.flush()
            if self.stop_flag:
                run_status = "stopped_by_user"
                print("\nAnalysis stopped by user")
                return

            combined_fasta = combine_iupac_fastas(output_dir)
            if not combined_fasta:
                run_status = "failed_no_valid_input"
                failure_message = "No valid IUPAC FASTA files were generated"
                print(f"\n{failure_message}")
                return

            alignment_file = os.path.join(output_dir, "alignment.aln")
            if not run_muscle_alignment(combined_fasta, alignment_file, muscle_path):
                run_status = "failed_muscle"
                failure_message = "MUSCLE alignment failed"
                print(f"\n{failure_message}")
                return

            # [FIX v1.0.0] Generate trace plots and heterozygote reports using
            # consensus sequence positions (deferred from Phase 1 to reduce memory pressure).
            # [FIX v1.2] For large batches, ONLY generate plots for samples with
            # actual heterozygote calls. Reading 400 ABI files + creating 400 TkAgg
            # figures exhausts Windows GDI handles and causes silent crash.
            try:
                alignment = AlignIO.read(alignment_file, "fasta")
                cons_pos_map = compute_consensus_position_mapping(alignment)
                sample_pos_maps = compute_sample_position_maps(alignment, cons_pos_map)
                parsed_records = [r for r in input_records if r.get("status") == "parsed"]
                total_parsed = len(parsed_records)

                # Only process samples with actual heterozygote calls (QC ok)
                het_records = [
                    r for r in parsed_records
                    if r.get("_heterozygotes") and any(
                        h.get("qc_status") == "ok" for h in r["_heterozygotes"]
                    )
                ]
                non_het_count = len(parsed_records) - len(het_records)

                if non_het_count > 0:
                    print(
                        f"\nSkipping trace plots for {non_het_count}/{total_parsed} "
                        f"samples (no heterozygous calls) to reduce memory pressure."
                    )

                if het_records:
                    print(
                        f"\nGenerating trace plots for {len(het_records)} samples "
                        f"with heterozygous calls..."
                    )
                else:
                    print(
                        f"\nNo heterozygous calls found across {total_parsed} samples; "
                        f"skipping trace plot generation."
                    )

                for idx, record in enumerate(het_records):
                    abi_file = record.get("_abi_file")
                    heterozygotes = record.get("_heterozygotes")
                    if not abi_file or heterozygotes is None:
                        continue
                    try:
                        pos_map = None
                        for align_rec in alignment:
                            if os.path.splitext(os.path.basename(abi_file))[0] in align_rec.id:
                                pos_map = sample_pos_maps.get(align_rec.id)
                                break
                        if pos_map is not None:
                            trace_data = parser.read_abi(abi_file)
                            if trace_data is not None:
                                save_heterozygote_report(heterozygotes, abi_file, output_dir,
                                                         position_map=pos_map)
                                plot_raw_trace_data(trace_data, abi_file, heterozygotes,
                                                    output_dir, position_map=pos_map)
                                # [MEMORY] Agg backend: zero GDI handles per plot.
                                # Only need to release numpy arrays + Python objects.
                                del trace_data
                                # Deep GC every 100 plots + flush output to avoid
                                # Windows console buffer overflow.
                                batch_count = idx + 1
                                total_het = len(het_records)
                                if batch_count % 100 == 0:
                                    gc.collect()
                                    sys.stdout.flush()
                                if batch_count % 10 == 0 or batch_count == total_het:
                                    print(
                                        f"  [{batch_count}/{total_het}] trace plots "
                                        f"generated"
                                    )
                                    sys.stdout.flush()
                    except Exception as exc:
                        print(f"  -> Failed plot for {os.path.basename(abi_file)}: {exc}")

                # [MEMORY] Clear stored heterozygotes to free memory before downstream steps
                for record in input_records:
                    record.pop("_heterozygotes", None)
                    record.pop("_abi_file", None)
                gc.collect()
                if het_records:
                    print("Trace plot generation complete.")
            except Exception as exc:
                print(f"Warning: Could not generate trace plots: {exc}")

            if ref_genome and os.path.exists(ref_genome):
                consensus_file = os.path.join(output_dir, "consensus_sequence.fasta")
                consensus_sequence = create_consensus_sequence(
                    alignment_file, consensus_file, quality_threshold
                )
                if consensus_sequence:
                    reference_alignment = os.path.join(
                        output_dir, "consensus_vs_reference.aln"
                    )
                    _, genome_mapping, ref_seq_name = align_consensus_to_reference(
                        consensus_file, ref_genome, reference_alignment, blastn_path
                    )
                    if genome_mapping:
                        genome_snp_mapping = genome_mapping
                        statuses = sorted(
                            {
                                item.get("mapping_status", "unknown")
                                for item in genome_mapping.values()
                            }
                        )
                        metadata["reference_mapping"] = {
                            "reference_sequence": ref_seq_name,
                            "mapped_position_count": len(genome_mapping),
                            "mapping_statuses": statuses,
                        }
                    else:
                        metadata["reference_mapping"] = {
                            "mapped_position_count": 0,
                            "mapping_statuses": ["unmapped"],
                        }

            genotype_matrix_file = os.path.join(output_dir, "genotype_matrix.xlsx")
            create_genotype_matrix(
                alignment_file,
                genotype_matrix_file,
                similarity_threshold,
                maf_threshold,
                min_obs_het,
                genome_snp_mapping,
                ref_seq_name,
            )
            consensus_file = os.path.join(output_dir, "consensus_sequence.fasta")
            create_consensus_sequence(alignment_file, consensus_file, quality_threshold)

            try:
                genotype_frame = pd.read_excel(
                    genotype_matrix_file, sheet_name="Genotype_Matrix"
                )
                self.sample_names = list(genotype_frame.iloc[:, 0].astype(str))
                self.snp_names = list(genotype_frame.columns[1:])
                self.genotype_matrix_data = {
                    snp: {
                        sample: str(genotype_frame.iloc[index][snp])
                        for index, sample in enumerate(self.sample_names)
                    }
                    for snp in self.snp_names
                }
                self.last_output_dir = output_dir
                # v1.1.0: keep the BLASTN mapping for the VCF export (CHROM /
                # POS / reference-forward REF>ALT per locus)
                self.genome_snp_mapping = genome_snp_mapping
                self.reference_sequence_name = ref_seq_name
                metadata["summary"]["retained_genotype_sample_count"] = len(
                    self.sample_names
                )
                metadata["summary"]["output_snp_count"] = len(self.snp_names)
                self._safe_gui_call(
                    lambda: self.pop_manager_button.config(state=tk.NORMAL)
                )
                self._safe_gui_call(
                    lambda: self.export_formats_button.config(state=tk.NORMAL)
                )
            except Exception as exc:
                print(f"Warning: Could not load genotype matrix: {exc}")
                metadata["genotype_matrix_load_warning"] = str(exc)

            run_status = "completed"
            print(
                f"\n===== Analysis complete - "
                f"{datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')} ====="
            )

            def _show_complete():
                messagebox.showinfo(
                    "Complete",
                    "Analysis completed successfully!\n\n"
                    "A reproducibility manifest was saved as analysis_metadata.json.",
                )

            self._safe_gui_call(_show_complete)
        except Exception as exc:
            run_status = "failed_exception"
            failure_message = str(exc)
            print(f"\nError during analysis: {exc}")
            import traceback

            traceback.print_exc()
            self._safe_gui_call(
                lambda message=failure_message: messagebox.showerror("Error", message)
            )
        finally:
            ended = datetime.now().astimezone()
            metadata["run"].update(
                {
                    "ended_at": ended.isoformat(),
                    "duration_seconds": round(time.perf_counter() - started_timer, 6),
                    "status": run_status,
                    "failure_message": failure_message,
                }
            )
            try:
                manifest = write_analysis_metadata(output_dir, metadata)
                print(f"Analysis metadata saved: {manifest}")
            except Exception as metadata_error:
                print(f"WARNING: Could not save analysis metadata: {metadata_error}")
            self._safe_gui_call(lambda: self.run_button.config(state=tk.NORMAL))
            self._safe_gui_call(lambda: self.stop_button.config(state=tk.DISABLED))


    # ===========================================================
    # v1.0 Population Genetics – Load Existing Genotype Data
    # ===========================================================

    def _load_genotype_data(self):
        """Load SNP genotypes from any supported genotype file.

        Accepts the native ``genotype_matrix.xlsx`` *and* every format written
        by "Export SNP Formats" (GENEPOP, Arlequin .arp, hierfstat, GenAlEx,
        FSTAT, POPGENE32, STRUCTURE).  The format is auto-detected; if the
        content is ambiguous the user is asked to pick one explicitly.
        """
        if not SNP_IO_AVAILABLE:
            messagebox.showerror("Error",
                "The SNP genotype import module (snp_format_io.py) is missing.\n"
                "Please keep it next to the program / reinstall the application.",
                parent=self.root)
            return

        path = filedialog.askopenfilename(
            title="Select SNP genotype file",
            filetypes=[
                ("All supported SNP genotype files",
                 "*.xlsx *.xls *.arp *.txt *.dat *.pop *.csv"),
                ("Genotype matrix (Excel, native)", "*.xlsx *.xls"),
                ("Arlequin project", "*.arp"),
                ("GENEPOP / hierfstat (text)", "*.txt"),
                ("GenAlEx (CSV)", "*.csv"),
                ("FSTAT data", "*.dat"),
                ("POPGENE32", "*.pop"),
                ("All files", "*.*")],
            initialdir=self.last_output_dir or ".")
        if not path:
            return
        self._import_genotype_file(path)

    def _import_genotype_file(self, path, forced_format=None):
        """Parse ``path`` into the in-memory genotype table and preview it."""
        filename = os.path.basename(path)
        try:
            try:
                data = snpio.read_genotype_file(path, forced_format)
            except snpio.AmbiguousFormatError:
                chosen = self._ask_import_format(filename)
                if not chosen:
                    return
                data = snpio.read_genotype_file(path, chosen)
        except Exception as e:
            messagebox.showerror("Error",
                f"Failed to read genotype file:\n{e}\n\n"
                "Supported: genotype_matrix.xlsx, GENEPOP, Arlequin (.arp), "
                "hierfstat, GenAlEx (CSV), FSTAT (.dat), POPGENE32 (.pop), "
                "STRUCTURE.",
                parent=self.root)
            return

        self.sample_names = list(data.samples)
        self.snp_names = list(data.loci)
        self.genotype_matrix_data = data.to_matrix_dict()
        self.populations = {pop: set(members)
                            for pop, members in data.populations.items()}
        self.last_output_dir = os.path.dirname(path) or self.last_output_dir
        # v1.1.0: an imported matrix carries no BLASTN mapping; the VCF writer
        # then falls back to the reference-forward SNP names or amplicon
        # coordinates.
        self.genome_snp_mapping = None
        self.reference_sequence_name = None

        print(f"[Population Genetics] Loaded genotype data from: {filename}")
        print(f"  Detected format : {snpio.FORMAT_LABELS.get(data.source_format, data.source_format)}")
        print(f"  Samples         : {len(self.sample_names)}")
        print(f"  SNPs            : {len(self.snp_names)}")
        if data.populations:
            print(f"  Populations     : {len(data.populations)} "
                  f"({', '.join(sorted(data.populations))})")
        for warning in data.warnings:
            print(f"  Note            : {warning}")

        self.pop_manager_button.config(state=tk.NORMAL)
        self.export_formats_button.config(state=tk.NORMAL)
        # A file that already carries >= 2 populations can go straight to the
        # statistical analysis; otherwise the Population Manager stays as the
        # required next step.
        self.stats_button.config(
            state=tk.NORMAL if len(self.populations) >= 2 else tk.DISABLED)

        self._show_data_preview(
            filename,
            source_format=data.source_format,
            warnings=data.warnings)

    def _ask_import_format(self, filename):
        """Let the user pick the genotype format when auto-detection fails."""
        choice = {"value": None}
        win = tk.Toplevel(self.root)
        win.title("Select Genotype File Format")
        win.geometry("520x420")
        win.minsize(420, 320)
        win.transient(self.root)

        ttk.Label(win, text=f"Could not identify the format of:\n{filename}",
                  wraplength=460, justify=tk.LEFT).pack(pady=(12, 4))
        ttk.Label(win, text="Please select the format to import:",
                  foreground="gray").pack(pady=(0, 8))

        var = tk.StringVar(value="genepop")
        frame = ttk.LabelFrame(win, text="Format", padding="10")
        frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=5)
        for fmt_id, label in snpio.FORMAT_LABELS.items():
            ttk.Radiobutton(frame, text=label, variable=var,
                            value=fmt_id).pack(anchor=tk.W, pady=2)

        def accept():
            choice["value"] = var.get()
            win.destroy()

        bottom = ttk.Frame(win)
        bottom.pack(fill=tk.X, pady=10, padx=15)
        ttk.Button(bottom, text="Import", command=accept).pack(side=tk.RIGHT, padx=5)
        ttk.Button(bottom, text="Cancel",
                   command=lambda: (choice.__setitem__("value", None),
                                    win.destroy())).pack(side=tk.RIGHT, padx=5)

        win.grab_set()
        win.wait_window()
        return choice["value"]

    def _show_data_preview(self, filename, source_format="", warnings=None):
        preview_win = tk.Toplevel(self.root)
        preview_win.title("Genotype Data Preview - Confirm Data")
        preview_win.geometry("900x550")
        preview_win.minsize(700, 400)
        preview_win.resizable(True, True)
        preview_win.transient(self.root)

        info_frame = ttk.Frame(preview_win)
        info_frame.pack(fill=tk.X, padx=15, pady=(15, 5))

        warnings = list(warnings or [])
        format_label = ""
        if source_format and SNP_IO_AVAILABLE:
            format_label = snpio.FORMAT_LABELS.get(source_format, source_format)

        ttk.Label(info_frame, text=f"Source file: {filename}",
                  font=("Microsoft YaHei", 10, "bold")).pack(anchor=tk.W)
        if format_label:
            ttk.Label(info_frame, text=f"Detected format: {format_label}",
                      font=("Microsoft YaHei", 9),
                      foreground="#1a6fb5").pack(anchor=tk.W, pady=(2, 0))
        pop_text = (f"   |   Populations: {len(self.populations)}"
                    if self.populations else "   |   Populations: none (use Population Manager)")
        ttk.Label(info_frame,
                  text=(f"Samples: {len(self.sample_names)}   |   "
                        f"SNPs/Loci: {len(self.snp_names)}{pop_text}"),
                  font=("Microsoft YaHei", 9)).pack(anchor=tk.W, pady=(3, 0))
        for warning in warnings:
            ttk.Label(info_frame, text=f"⚠ {warning}",
                      font=("Microsoft YaHei", 9), foreground="darkorange",
                      wraplength=820, justify=tk.LEFT).pack(anchor=tk.W, pady=(2, 0))

        table_frame = ttk.Frame(preview_win)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        n_show_samples = min(len(self.sample_names), 15)
        n_show_snps = min(len(self.snp_names), 10)

        cols = ["Sample"] + self.snp_names[:n_show_snps]
        tree = ttk.Treeview(table_frame, columns=cols, show="headings")
        tree.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.E, tk.W))

        vsb = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=tree.yview)
        vsb.grid(row=0, column=1, sticky=(tk.N, tk.S))
        hsb = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=tree.xview)
        hsb.grid(row=1, column=0, sticky=(tk.E, tk.W))
        tree.config(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tree.column("Sample", width=180, minwidth=120)
        for snp in self.snp_names[:n_show_snps]:
            tree.column(snp, width=70, minwidth=50, anchor=tk.CENTER)
            tree.heading(snp, text=snp)

        for sample in self.sample_names[:n_show_samples]:
            vals = [sample]
            for snp in self.snp_names[:n_show_snps]:
                vals.append(self.genotype_matrix_data.get(snp, {}).get(sample, ""))
            tree.insert("", tk.END, values=vals)

        note = ""
        if len(self.sample_names) > n_show_samples or len(self.snp_names) > n_show_snps:
            note = f"(Showing first {n_show_samples} of {len(self.sample_names)} samples "
            if len(self.snp_names) > n_show_snps:
                note += f"and first {n_show_snps} of {len(self.snp_names)} loci)"
            else:
                note += ")"
            ttk.Label(table_frame, text=note, foreground="gray").grid(
                row=2, column=0, sticky=tk.W, pady=(5, 0))

        btn_frame = ttk.Frame(preview_win)
        btn_frame.pack(fill=tk.X, padx=15, pady=(5, 15))

        confirmed = tk.BooleanVar(value=False)

        def accept():
            confirmed.set(True)
            preview_win.destroy()

        ttk.Button(btn_frame, text="Confirm — Data is Correct", command=accept).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel — Load Different File",
                   command=preview_win.destroy).pack(side=tk.RIGHT, padx=5)

        ttk.Label(btn_frame, text="⚠ Please verify the genotype data above before proceeding.",
                  foreground="darkorange").pack(side=tk.LEFT)

        preview_win.wait_window()

        if not confirmed.get():
            self.genotype_matrix_data = None
            self.sample_names = []
            self.snp_names = []
            self.populations = {}
            self.pop_manager_button.config(state=tk.DISABLED)
            self.export_formats_button.config(state=tk.DISABLED)
            self.stats_button.config(state=tk.DISABLED)
            print("[Population Genetics] Data preview was cancelled — data unloaded.")
        else:
            print("[Population Genetics] Data preview confirmed.")

    # ===========================================================
    # v1.0 Population Genetics – Export SNP Genotypes
    # ===========================================================

    def _export_snp_formats(self):
        if not self.genotype_matrix_data:
            result = messagebox.askyesno("No Data",
                "No genotype data loaded.\n\n"
                "Yes - Load an existing genotype file\n"
                "No  - Cancel",
                parent=self.root)
            if result:
                self._load_genotype_data()
                if not self.genotype_matrix_data:
                    return
            else:
                return

        fmt_win = tk.Toplevel(self.root)
        fmt_win.title("Export SNP Genotypes - Format Selection")
        fmt_win.geometry("760x600")
        fmt_win.minsize(600, 400)
        fmt_win.resizable(True, True)
        fmt_win.transient(self.root)

        ttk.Label(fmt_win, text="Select export format for SNP genotype data",
                  font=("Microsoft YaHei", 11, "bold")).pack(pady=10)

        fmt_var = tk.StringVar(value="genepop")

        formats_frame = ttk.LabelFrame(fmt_win, text="Available Formats", padding="10")
        formats_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)

        fmts = [
            ("Genepop Format (*.txt)", "genepop",
             "Classic population genetics format.",
             "Compatible with Fstat, Genepop, and most analysis software."),
            ("Arlequin Format (*.arp)", "arlequin",
             "Arlequin 3.5 project file.",
             "Direct import into Arlequin for AMOVA, Fst, etc."),
            ("hierfstat Format (*.txt)", "hierfstat",
             "R table read by read.table(header=TRUE).",
             "Column 1 = population; one column per locus (A=1,T=2,C=3,G=4; 00 = missing)."),
            ("GenAlEx Format (*.csv)", "genalex",
             "GenAlEx 6.5 diploid codominant layout.",
             "3 header rows + 2 allele columns per locus; opens directly in Excel."),
            ("STRUCTURE Format (*.txt)", "structure",
             "Two-column allele integer format.",
             "Each SNP = 2 columns (allele1, allele2). Missing = -9."),
            ("FSTAT Format (*.dat)", "fstat",
             "FSTAT native format.",
             "Includes locus names and population structure headers."),
            ("POPGENE32 Format (*.pop)", "popgene",
             "POPGENE 1.32 diploid native format.",
             "2-letter genotype codes (AA, AT, TT...). Missing = .."),
            ("VCF Format (*.vcf)", "vcf",
             "Variant Call Format 4.2, one sample column per individual.",
             "CHROM/POS from the unique BLASTN mapping; 0/0 0/1 1/1 genotypes, "
             "'./.' = missing. Readable by bcftools, vcftools, tabix."),
        ]

        for i, (label, value, desc1, desc2) in enumerate(fmts):
            ttk.Radiobutton(formats_frame, text=label, variable=fmt_var,
                           value=value).grid(row=i, column=0, sticky=tk.W, pady=(8, 2))
            full_desc = desc1 + " " + desc2
            desc_lbl = ttk.Label(formats_frame, text=full_desc, foreground="gray",
                                wraplength=550)
            desc_lbl.grid(row=i, column=1, sticky=tk.W, padx=(20, 0))

        formats_frame.columnconfigure(1, weight=1)

        include_pop = tk.BooleanVar(value=True)
        ttk.Checkbutton(fmt_win, text="Include population assignment info (if available)",
                       variable=include_pop).pack(pady=5, padx=20, anchor=tk.W)

        ttk.Label(fmt_win, text=f"Samples: {len(self.sample_names)}  |  SNPs: {len(self.snp_names)}",
                  foreground="gray").pack(pady=5)

        bottom_btns = ttk.Frame(fmt_win)
        bottom_btns.pack(fill=tk.X, pady=10, padx=20)

        def do_export():
            fmt = fmt_var.get()
            ext_map = {"popgene": ".pop", "fstat": ".dat", "arlequin": ".arp",
                       "genepop": ".txt", "structure": ".txt",
                       "hierfstat": ".txt", "genalex": ".csv", "vcf": ".vcf"}
            expected_ext = ext_map.get(fmt, ".txt")
            # The Windows save dialog appends the extension of the *selected*
            # filter to a bare filename, so the filter list must lead with the
            # format's own extension — otherwise every export is saved .txt.
            format_labels = {
                "genepop": "GENEPOP files",
                "arlequin": "Arlequin project files",
                "hierfstat": "hierfstat files",
                "genalex": "GenAlEx CSV files",
                "vcf": "VCF files",
                "structure": "STRUCTURE files",
                "fstat": "FSTAT data files",
                "popgene": "POPGENE files",
            }
            own_label = format_labels.get(fmt, f"{fmt.upper()} files")
            path = filedialog.asksaveasfilename(
                defaultextension=expected_ext,
                initialfile=f"SNP_genotypes{expected_ext}",
                filetypes=[(f"{own_label} (*{expected_ext})", f"*{expected_ext}"),
                           ("All files", "*.*")],
                title=f"Save as {fmt.upper()} Format",
                initialdir=self.last_output_dir or ".")
            if not path:
                return
            # Enforce correct extension (asksaveasfilename defaultextension is unreliable)
            if not path.lower().endswith(expected_ext):
                path += expected_ext
            try:
                use_pop = include_pop.get() and bool(self.populations) and len(self.populations) >= 1
                if fmt == "genepop":
                    n = self._write_genepop_snp(path, use_pop)
                elif fmt == "arlequin":
                    n = self._write_arlequin_snp(path, use_pop)
                elif fmt == "hierfstat":
                    n = self._write_hierfstat_snp(path, use_pop)
                elif fmt == "genalex":
                    n = self._write_genalex_snp(path, use_pop)
                elif fmt == "structure":
                    n = self._write_structure_snp(path, use_pop)
                elif fmt == "fstat":
                    n = self._write_fstat_snp(path, use_pop)
                elif fmt == "popgene":
                    n = self._write_popgene_snp(path, use_pop)
                elif fmt == "vcf":
                    n = self._write_vcf_snp(path, use_pop)
                else:
                    return
                messagebox.showinfo("Exported",
                    f"Successfully exported {n} samples.\nSaved to:\n{path}", parent=fmt_win)
                fmt_win.destroy()
            except Exception as e:
                error_msg = str(e)
                messagebox.showerror("Error", f"Export failed:\n{error_msg}", parent=fmt_win)

        ttk.Button(bottom_btns, text="Export", command=do_export).pack(side=tk.RIGHT, padx=5)
        ttk.Button(bottom_btns, text="Cancel", command=fmt_win.destroy).pack(side=tk.RIGHT, padx=5)

    # ---- Format writers for SNP genotypes ----
    # Every writer delegates to ``snp_format_io`` so that export and import
    # share one implementation and the round trip stays loss-less (verified by
    # tests/test_snp_format_io.py).

    def _export_populations(self, use_pop):
        """Population mapping handed to a writer (``None`` = single group)."""
        if use_pop and self.populations:
            return {str(p): sorted(s) for p, s in self.populations.items()}
        return None

    def _require_snp_io(self):
        if not SNP_IO_AVAILABLE:
            raise RuntimeError(
                "The SNP genotype export module (snp_format_io.py) is missing.")

    @staticmethod
    def _snp_allele_code(allele):
        mapping = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 'a': 1, 't': 2, 'c': 3, 'g': 4}
        return mapping.get(allele, 0)

    @staticmethod
    def _parse_diploid(gt):
        if not gt or gt in ("--", "NN", "N/A", "nan"):
            return (None, None)
        gt_clean = gt.strip().upper()
        if not gt_clean or gt_clean in ("--", "NN", "N/A", "NAN"):
            return (None, None)
        if len(gt_clean) == 1 and gt_clean in "ACGT":
            return (gt_clean, gt_clean)
        if len(gt_clean) >= 2:
            a1 = gt_clean[0] if gt_clean[0] in "ACGT" else None
            a2 = gt_clean[1] if gt_clean[1] in "ACGT" else None
            return (a1, a2)
        return (None, None)

    def _write_genepop_snp(self, path, use_pop):
        """Export SNP genotypes in GENEPOP 4.x native format.

        ``snp_format_io`` emits a title line carrying the allele coding and the
        population names so the file can be re-imported byte-for-byte.
        """
        self._require_snp_io()
        return snpio.write_genepop(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop),
            title=f"Genotype data - ABI-SNP Analyzer v{SOFTWARE_VERSION}")

    def _write_arlequin_snp(self, path, use_pop):
        """Export SNP genotypes in Arlequin 3.5 project (.arp) format."""
        self._require_snp_io()
        return snpio.write_arlequin(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop),
            title=f"ABI-SNP Analyzer v{SOFTWARE_VERSION} - SNP genotypes")

    def _write_hierfstat_snp(self, path, use_pop):
        """Export SNP genotypes as a hierfstat / R genotype table.

        ``read.table(file, header=TRUE)`` in R yields one column per locus plus
        a leading ``pop`` column; allele codes are A=1, T=2, C=3, G=4 with
        ``00`` marking a missing call, which is what ``hierfstat::basic.stats()``
        and ``hierfstat::wc()`` expect.
        """
        self._require_snp_io()
        return snpio.write_hierfstat(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop))

    def _write_genalex_snp(self, path, use_pop):
        """Export SNP genotypes in GenAlEx 6.5 diploid codominant layout (CSV)."""
        self._require_snp_io()
        return snpio.write_genalex(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop),
            title=f"ABI-SNP Analyzer v{SOFTWARE_VERSION} - SNP genotypes")

    def _write_structure_snp(self, path, use_pop):
        self._require_snp_io()
        return snpio.write_structure(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop))

    def _write_fstat_snp(self, path, use_pop):
        """Export SNP genotypes in FSTAT 2.9.3/2.9.4 native format (``.dat``)."""
        self._require_snp_io()
        return snpio.write_fstat(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop))

    def _write_popgene_snp(self, path, use_pop):
        self._require_snp_io()
        return snpio.write_popgene(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            self._export_populations(use_pop),
            title=f"SNP Genotype Data - ABI-SNP Analyzer v{SOFTWARE_VERSION}")

    def _write_vcf_snp(self, path, use_pop):
        """Export SNP genotypes in VCF 4.2 format (v1.1.0).

        CHROM/POS/REF/ALT come from the unique BLASTN mapping when the
        reference genome was analysed (alleles are complemented onto the
        reference-forward strand); reference-forward SNP names such as
        ``NC_031982.2:26437556C>T`` are parsed directly; loci without any
        mapping fall back to the virtual ``amplicon`` contig with POS = the
        consensus coordinate (declared in the VCF header).  Population
        structure is not part of the VCF specification and is ignored.
        """
        self._require_snp_io()
        return snpio.write_vcf(
            path, self.sample_names, self.snp_names, self.genotype_matrix_data,
            locus_metadata=snpio.vcf_locus_metadata_from_mapping(
                self.snp_names, self.genotype_matrix_data,
                getattr(self, "genome_snp_mapping", None)),
            reference_name=getattr(self, "reference_sequence_name", None),
            source=f"ABI-SNP_Analyzer_v{SOFTWARE_VERSION}")

    # ===========================================================
    # v1.0 Population Genetics – Population Manager
    # ===========================================================

    def _open_population_manager(self):
        if not self.genotype_matrix_data:
            result = messagebox.askyesno("No Data",
                "No genotype data loaded yet.\n\n"
                "Yes - Load an existing genotype file\n"
                "No  - Cancel",
                parent=self.root)
            if result:
                self._load_genotype_data()
                if not self.genotype_matrix_data:
                    return
            else:
                return

        win = tk.Toplevel(self.root)
        win.title("Population Manager")
        win.geometry("1200x750")
        win.minsize(800, 500)
        win.resizable(True, True)
        win.transient(self.root)

        ttk.Label(win, text="Assign samples to populations for statistical analysis.",
                  font=("Microsoft YaHei", 10)).pack(pady=5)

        main = ttk.Frame(win, padding="10")
        main.pack(fill=tk.BOTH, expand=True)

        left_frame = ttk.LabelFrame(main, text="Available Samples", padding="5")
        left_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(0, weight=1)
        left_frame.rowconfigure(1, weight=0)

        avail_listbox = tk.Listbox(left_frame, selectmode=tk.EXTENDED, width=45, height=18)
        avail_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        avail_scroll_v = ttk.Scrollbar(left_frame, orient=tk.VERTICAL, command=avail_listbox.yview)
        avail_scroll_v.grid(row=0, column=1, sticky=(tk.N, tk.S))
        avail_listbox.config(yscrollcommand=avail_scroll_v.set)
        avail_scroll_h = ttk.Scrollbar(left_frame, orient=tk.HORIZONTAL, command=avail_listbox.xview)
        avail_scroll_h.grid(row=1, column=0, sticky=(tk.W, tk.E))
        avail_listbox.config(xscrollcommand=avail_scroll_h.set)
        for s in self.sample_names:
            avail_listbox.insert(tk.END, s)

        mid_frame = ttk.Frame(main)
        mid_frame.grid(row=0, column=1, padx=10)

        ttk.Label(mid_frame, text="Population:").grid(row=0, column=0, pady=5)
        pop_var = tk.StringVar(value="Pop1")
        ttk.Entry(mid_frame, textvariable=pop_var, width=15).grid(row=0, column=1, pady=5)

        def refresh_pop_list():
            pop_listbox.delete(0, tk.END)
            for pname in sorted(self.populations.keys()):
                for s in sorted(self.populations[pname]):
                    pop_listbox.insert(tk.END, f"{pname} -> {s}")

        def add_to_pop():
            sel = avail_listbox.curselection()
            if not sel:
                return
            pname = pop_var.get().strip()
            if not pname:
                messagebox.showwarning("Warning", "Please enter a population name.", parent=win)
                return
            if pname not in self.populations:
                self.populations[pname] = set()
            for idx in sel:
                self.populations[pname].add(avail_listbox.get(idx))
            refresh_pop_list()
            avail_listbox.selection_clear(0, tk.END)

        def remove_from_pop():
            sel = pop_listbox.curselection()
            if not sel:
                return
            for idx in reversed(sel):
                parts = pop_listbox.get(idx).split(" -> ")
                if len(parts) == 2:
                    pname, sname = parts
                    if pname in self.populations:
                        self.populations[pname].discard(sname)
                        if not self.populations[pname]:
                            del self.populations[pname]
            refresh_pop_list()

        def import_csv():
            path = filedialog.askopenfilename(
                title="Import Population Assignment CSV",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
                parent=win)
            if not path:
                return
            try:
                imported = 0
                with open(path, newline="", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        sample = row.get("Sample", "").strip()
                        pop = row.get("Population", "").strip()
                        if sample and pop:
                            self.populations.setdefault(pop, set()).add(sample)
                            imported += 1
                refresh_pop_list()
                messagebox.showinfo("Imported", f"{imported} rows imported.", parent=win)
            except Exception as e:
                error_msg = str(e)
                messagebox.showerror("Error", f"Import failed: {error_msg}", parent=win)

        ttk.Button(mid_frame, text="Add ->", command=add_to_pop).grid(
            row=1, column=0, columnspan=2, pady=5, sticky=tk.EW)
        ttk.Button(mid_frame, text="<- Remove", command=remove_from_pop).grid(
            row=2, column=0, columnspan=2, pady=5, sticky=tk.EW)
        ttk.Separator(mid_frame, orient=tk.HORIZONTAL).grid(
            row=3, column=0, columnspan=2, sticky=tk.EW, pady=10)
        ttk.Button(mid_frame, text="Import CSV", command=import_csv).grid(
            row=4, column=0, columnspan=2, pady=5, sticky=tk.EW)
        ttk.Button(mid_frame, text="Clear All",
                   command=lambda: (self.populations.clear(), refresh_pop_list())).grid(
            row=5, column=0, columnspan=2, pady=5, sticky=tk.EW)

        right_frame = ttk.LabelFrame(main, text="Current Assignments", padding="5")
        right_frame.grid(row=0, column=2, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)
        right_frame.rowconfigure(1, weight=0)

        pop_listbox = tk.Listbox(right_frame, selectmode=tk.EXTENDED, width=50, height=18)
        pop_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        pop_scroll_v = ttk.Scrollbar(right_frame, orient=tk.VERTICAL, command=pop_listbox.yview)
        pop_scroll_v.grid(row=0, column=1, sticky=(tk.N, tk.S))
        pop_listbox.config(yscrollcommand=pop_scroll_v.set)
        pop_scroll_h = ttk.Scrollbar(right_frame, orient=tk.HORIZONTAL, command=pop_listbox.xview)
        pop_scroll_h.grid(row=1, column=0, sticky=(tk.W, tk.E))
        pop_listbox.config(xscrollcommand=pop_scroll_h.set)

        refresh_pop_list()

        main.columnconfigure(0, weight=1)
        main.columnconfigure(1, weight=0)
        main.columnconfigure(2, weight=1)
        main.rowconfigure(0, weight=1)

        bottom = ttk.Frame(win)
        bottom.pack(fill=tk.X, pady=5)

        def done():
            total = sum(len(v) for v in self.populations.values())
            if total < 2:
                messagebox.showwarning("Warning",
                    "At least 2 samples must be assigned to populations.", parent=win)
                return
            if len(self.populations) < 2:
                messagebox.showwarning("Warning",
                    "At least 2 populations are required for statistical analysis.", parent=win)
                return
            self.stats_button.config(state=tk.NORMAL)
            messagebox.showinfo("Success",
                f"{total} samples assigned to {len(self.populations)} populations.", parent=win)
            win.destroy()

        def export_assignments():
            if not self.populations:
                messagebox.showwarning("Warning", "No populations assigned yet.", parent=win)
                return
            path = filedialog.asksaveasfilename(defaultextension=".csv",
                filetypes=[("CSV", "*.csv")], title="Export Population Assignments",
                initialdir=self.last_output_dir or ".")
            if not path:
                return
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["Sample", "Population"])
                for pname in sorted(self.populations.keys()):
                    for sname in sorted(self.populations[pname]):
                        w.writerow([sname, pname])
            messagebox.showinfo("Saved", f"Assignments saved to: {path}", parent=win)

        ttk.Button(bottom, text="Done", command=done).pack(side=tk.RIGHT, padx=10)
        ttk.Button(bottom, text="Export CSV", command=export_assignments).pack(side=tk.RIGHT, padx=5)
        ttk.Label(bottom, text=f"Total samples: {len(self.sample_names)}").pack(side=tk.LEFT, padx=10)

    # ===========================================================
    # v1.0 Population Genetics – Run Statistical Analysis
    # ===========================================================

    def _run_statistical_analysis(self):
        if not POPGEN_AVAILABLE:
            messagebox.showerror("Error",
                "population_genetics module not available.\n"
                "Please ensure population_genetics.py is in the same folder as this program.")
            return
        if not self.genotype_matrix_data:
            result = messagebox.askyesno("No Data",
                "No genotype data loaded yet.\n\n"
                "Yes - Load an existing genotype file\n"
                "No  - Cancel",
                parent=self.root)
            if result:
                self._load_genotype_data()
                if not self.genotype_matrix_data:
                    return
            else:
                return
        if len(self.populations) < 2:
            messagebox.showerror("Error",
                "Please assign samples to at least 2 populations using Population Manager first.")
            return

        try:
            (
                all_data,
                sample_genotypes,
                sample_to_population,
                unassigned_samples,
            ) = build_population_analysis_inputs(
                self.genotype_matrix_data,
                self.snp_names,
                self.populations,
            )
        except ValueError as exc:
            messagebox.showerror("Population Assignment Error", str(exc))
            return

        if unassigned_samples:
            print(
                f"[Population Genetics] {len(unassigned_samples)} unassigned "
                "sample(s) excluded from statistical analysis."
            )

        if not all_data:
            messagebox.showerror("Error",
                "No valid SNP data with at least 2 populations found.")
            return

        # ---- how should linked SNPs (same amplicon) be handled? ---------
        if not LINKAGE_AVAILABLE:
            messagebox.showwarning(
                "Linkage module missing",
                "linkage_analysis.py was not found next to the program, so "
                "the linked-SNP handling modes are unavailable.\n"
                "The analysis will run in 'All SNPs' mode.",
                parent=self.root)
            mode_cfg = {"mode": "all", "r2_threshold": 0.8,
                        "ld_permutations": 999, "write_report": False}
        else:
            mode_cfg = self._ask_analysis_mode(len(all_data))
        if not mode_cfg:
            return

        prog_win = tk.Toplevel(self.root)
        prog_win.title("Running Statistical Analysis...")
        prog_win.geometry("420x160")
        prog_win.transient(self.root)
        ttk.Label(prog_win, text="Running population genetics analysis...").pack(pady=15)
        prog_bar = ttk.Progressbar(prog_win, mode="determinate", maximum=100)
        prog_bar.pack(fill=tk.X, padx=20, pady=8)
        prog_lbl = ttk.Label(prog_win, text="Initializing...")
        prog_lbl.pack(pady=5)
        prog_win.update()

        def _set_progress(value, text):
            if prog_win.winfo_exists():
                prog_bar["value"] = value
                prog_lbl.config(text=text)

        def _worker():
            try:
                n_permutations = 999
                random_seed = 42
                loci_all = list(all_data.keys())

                def _subset_genotypes(locus_list):
                    return {
                        sample: {locus: calls[locus]
                                 for locus in locus_list if locus in calls}
                        for sample, calls in sample_genotypes.items()
                    }

                def _run_mode(locus_list, p0, p1, label):
                    """Run the full statistical pipeline on a locus subset."""
                    subset = {locus: all_data[locus]
                              for locus in locus_list if locus in all_data}
                    total = len(subset)
                    locus_summaries = {}
                    for i, snp_name in enumerate(subset):
                        progress = p0 + (i / max(total, 1)) * (p1 - p0) * 0.75
                        label_i = f"{label}: {snp_name} ({i + 1}/{total})..."
                        self._safe_gui_call(
                            lambda v=progress, t=label_i: _set_progress(v, t)
                        )
                        locus_summaries[snp_name] = PopulationGenetics.analyze_locus(
                            subset[snp_name],
                            snp_name,
                            n_permutations=n_permutations,
                            random_seed=random_seed,
                        )

                    self._safe_gui_call(
                        lambda: _set_progress(
                            p0 + (p1 - p0) * 0.82,
                            f"{label}: pairwise-difference AMOVA (Excoffier 1992)...")
                    )
                    amova_results = []
                    amova_warning = None
                    try:
                        amova_results = PopulationGenetics.amova(
                            _subset_genotypes(list(subset)),
                            sample_to_population,
                            n_permutations=n_permutations,
                            random_seed=random_seed,
                            min_common_loci=1,
                        )
                    except ValueError as exc:
                        amova_warning = str(exc)
                        print(f"[AMOVA] Not calculated: {amova_warning}")

                    self._safe_gui_call(
                        lambda: _set_progress(
                            p0 + (p1 - p0) * 0.95,
                            f"{label}: genotype frequency tests...")
                    )
                    gft_results = PopulationGenetics.test_all_loci_genotype_freq(subset)
                    return (subset, locus_summaries, amova_results,
                            gft_results, amova_warning)

                def _estimate_linkage(progress_span):
                    """LD report + one representative SNP per linkage block."""
                    if not LINKAGE_AVAILABLE:
                        raise RuntimeError(
                            "The linkage module (linkage_analysis.py) is missing.")
                    p_start, p_end = progress_span

                    def ld_progress(fraction, text):
                        self._safe_gui_call(
                            lambda: _set_progress(
                                p_start + fraction * (p_end - p_start),
                                f"Linkage disequilibrium: {text}"))

                    pairs = linkage.compute_pairwise_ld(
                        sample_genotypes, loci_all, sample_to_population,
                        n_permutations=mode_cfg["ld_permutations"],
                        random_seed=random_seed,
                        progress_cb=ld_progress)
                    summaries = linkage.summarise_ld_pairs(pairs)
                    blocks = linkage.cluster_linkage_blocks(
                        loci_all, summaries, mode_cfg["r2_threshold"])
                    linkage.choose_representatives(
                        sample_genotypes, sorted(sample_genotypes), blocks)
                    kept = [locus for locus in loci_all
                            if locus in {b.representative for b in blocks}]

                    report_path = None
                    if mode_cfg["write_report"]:
                        base_dir = self.last_output_dir or os.getcwd()
                        report_path = os.path.join(base_dir, "linkage_report.xlsx")
                        try:
                            linkage.write_linkage_report(
                                report_path, pairs, summaries, blocks, kept,
                                r2_threshold=mode_cfg["r2_threshold"],
                                n_permutations=mode_cfg["ld_permutations"],
                                random_seed=random_seed)
                            print(f"[Linkage] LD report saved: {report_path}")
                        except Exception as report_error:
                            print(f"[Linkage] Could not save LD report: "
                                  f"{report_error}")
                            report_path = None

                    n_multi = sum(1 for b in blocks if len(b.loci) > 1)
                    print(f"[Linkage] {len(loci_all)} SNPs -> {len(blocks)} "
                          f"linkage blocks ({n_multi} with >1 SNP, "
                          f"r2 >= {mode_cfg['r2_threshold']}); "
                          f"{len(kept)} independent loci kept.")
                    return {
                        "pairs": pairs,
                        "summaries": summaries,
                        "blocks": blocks,
                        "kept_loci": kept,
                        "all_loci": loci_all,
                        "r2_threshold": mode_cfg["r2_threshold"],
                        "n_permutations": mode_cfg["ld_permutations"],
                        "report_path": report_path,
                        "available": True,
                    }

                mode_label = {
                    "all": "All SNPs (independent-loci assumption)",
                    "representative": "One representative SNP per linkage block",
                }.get(mode_cfg["mode"], "")

                linkage_info = None
                comparison = None

                if mode_cfg["mode"] == "all":
                    main_results = _run_mode(loci_all, 0, 92, "All SNPs")
                else:
                    linkage_info = _estimate_linkage((2, 10))
                    kept = linkage_info["kept_loci"] or loci_all
                    if mode_cfg["mode"] == "representative":
                        main_results = _run_mode(kept, 10, 92, "Representative SNPs")
                    else:                                   # compare both modes
                        results_all = _run_mode(loci_all, 10, 50, "Mode A (all SNPs)")
                        results_kept = _run_mode(kept, 52, 92, "Mode B (representatives)")
                        main_results = results_kept
                        comparison = {
                            "Mode A - all SNPs":
                                self._mode_summary(results_all),
                            "Mode B - one SNP per block":
                                self._mode_summary(results_kept),
                        }

                results_data = main_results

                def _finish(results=results_data):
                    _set_progress(100, "Done!")
                    if prog_win.winfo_exists():
                        prog_win.destroy()
                    self._show_analysis_results(
                        *results,
                        mode_label=mode_label,
                        linkage_info=linkage_info,
                        mode_comparison=comparison)

                self._safe_gui_call(_finish)
            except Exception as exc:
                import traceback

                traceback.print_exc()
                error_message = str(exc)

                def _show_error(message=error_message):
                    if prog_win.winfo_exists():
                        prog_win.destroy()
                    messagebox.showerror("Error", f"Analysis failed: {message}")

                self._safe_gui_call(_show_error)

        threading.Thread(target=_worker, daemon=True).start()

    @staticmethod
    def _mode_summary(results):
        """Compact per-mode metrics used by the mode-comparison table."""
        _all_data, locus_summaries, amova_results, _gft, _warning = results
        fsts = [s.fst for s in locus_summaries.values() if s.fst is not None]

        def _mean(values):
            values = [v for v in values if v is not None]
            return sum(values) / len(values) if values else None

        phi_st = None
        phi_p = None
        if amova_results:
            phi_st = amova_results[0].phi_stat
            phi_p = amova_results[0].p_value
        return {
            "n_loci": len(locus_summaries),
            "mean_ho": _mean([s.overall_ho for s in locus_summaries.values()]),
            "mean_he": _mean([s.overall_he for s in locus_summaries.values()]),
            "mean_pic": _mean([s.overall_pic for s in locus_summaries.values()]),
            "mean_fst": _mean(fsts),
            "amova_phi_st": phi_st,
            "amova_p_value": phi_p,
        }

    def _ask_analysis_mode(self, n_loci):
        """Ask how SNPs that share an amplicon (i.e. are linked) are handled.

        Returns ``None`` on cancel, otherwise a dict with the keys ``mode``
        (``all`` | ``representative`` | ``both``), ``r2_threshold``,
        ``ld_permutations`` and ``write_report``.
        """
        if n_loci < 2:
            return {"mode": "all", "r2_threshold": 0.8,
                    "ld_permutations": 999, "write_report": False}

        choice = {"value": None}
        win = tk.Toplevel(self.root)
        win.title("Linked SNPs - Analysis Mode")
        win.geometry("680x440")
        win.minsize(560, 380)
        win.transient(self.root)

        ttk.Label(win, text=(
            f"{n_loci} SNPs will be analysed.\n"
            "SNPs called from the same amplicon are physically linked and can "
            "be in strong linkage disequilibrium (LD), which inflates the "
            "degrees of freedom of multi-locus tests and violates the "
            "independent-loci assumption of software such as STRUCTURE.\n"
            "Please choose how they should be handled:"),
            justify=tk.LEFT, wraplength=620).pack(padx=15, pady=(12, 8))

        mode_var = tk.StringVar(value="representative")
        frame = ttk.LabelFrame(win, text="Analysis mode", padding="10")
        frame.pack(fill=tk.X, padx=15, pady=5)
        options = [
            ("all", "All SNPs (independent-loci assumption)",
             "Every SNP is kept. Over-states the information content when SNPs share an amplicon."),
            ("representative", "One representative SNP per linkage block (recommended)",
             "SNPs in LD (r2 >= threshold) are clustered and only the most reliable SNP of each block is analysed."),
            ("both", "Run both modes and compare",
             "Produces the LD report plus a side-by-side comparison of the two modes."),
        ]
        for value, label, desc in options:
            ttk.Radiobutton(frame, text=label, variable=mode_var,
                            value=value).pack(anchor=tk.W, pady=(4, 0))
            ttk.Label(frame, text=desc, foreground="gray",
                      wraplength=600).pack(anchor=tk.W, padx=(24, 0))

        params = ttk.Frame(win)
        params.pack(fill=tk.X, padx=15, pady=(8, 0))
        ttk.Label(params, text="r2 threshold:").grid(row=0, column=0, sticky=tk.W)
        r2_var = tk.StringVar(value="0.80")
        ttk.Spinbox(params, from_=0.5, to=1.0, increment=0.05, width=6,
                    textvariable=r2_var).grid(row=0, column=1, padx=(4, 20))
        ttk.Label(params, text="LD permutations:").grid(row=0, column=2, sticky=tk.W)
        perm_var = tk.StringVar(value="999")
        ttk.Spinbox(params, from_=99, to=9999, increment=100, width=7,
                    textvariable=perm_var).grid(row=0, column=3, padx=4)

        report_var = tk.BooleanVar(value=True)
        report_check = ttk.Checkbutton(
            win, text="Write LD / linkage-block report (linkage_report.xlsx)",
            variable=report_var)
        report_check.pack(anchor=tk.W, padx=15, pady=(10, 0))

        def accept():
            try:
                threshold = float(r2_var.get())
                if not 0.5 <= threshold <= 1.0:
                    raise ValueError
            except ValueError:
                messagebox.showerror("Error", "r2 threshold must be between 0.50 and 1.00.",
                                     parent=win)
                return
            try:
                permutations = int(perm_var.get())
                if permutations < 99:
                    raise ValueError
            except ValueError:
                messagebox.showerror("Error", "LD permutations must be an integer >= 99.",
                                     parent=win)
                return
            choice["value"] = {"mode": mode_var.get(), "r2_threshold": threshold,
                               "ld_permutations": permutations,
                               "write_report": report_var.get()}
            win.destroy()

        bottom = ttk.Frame(win)
        bottom.pack(fill=tk.X, pady=10, padx=15)
        ttk.Button(bottom, text="Run Analysis", command=accept).pack(side=tk.RIGHT, padx=5)
        ttk.Button(bottom, text="Cancel", command=win.destroy).pack(side=tk.RIGHT, padx=5)

        win.grab_set()
        win.wait_window()
        return choice["value"]

    # ===========================================================
    # v1.0 Population Genetics – Show Results Window
    # ===========================================================

    def _show_analysis_results(
        self, all_data, locus_summaries, amova_results, gft_results,
        amova_warning=None, mode_label=None, linkage_info=None,
        mode_comparison=None,
    ):
        win = tk.Toplevel(self.root)
        win.title("Population Genetics Statistical Results")
        win.geometry("1200x800")
        win.minsize(900, 600)
        win.resizable(True, True)
        win.transient(self.root)

        notebook = ttk.Notebook(win)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # ---------- Tab 1: Data Overview ----------
        tab1 = ttk.Frame(notebook)
        notebook.add(tab1, text="Data Overview")
        ttk.Label(tab1, text="Genotype Data Overview",
                  font=("Microsoft YaHei", 12, "bold")).pack(pady=5)
        ov_tree = ttk.Treeview(tab1,
            columns=("Population", "N_Samples", "N_SNPs"), show="headings")
        for c in ("Population", "N_Samples", "N_SNPs"):
            ov_tree.heading(c, text=c)
            ov_tree.column(c, width=150)
        ov_tree.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        for pname in sorted(self.populations.keys()):
            ov_tree.insert("", tk.END,
                values=(pname, len(self.populations[pname]), len(all_data)))

        # ---------- Tab 2: Summary Statistics ----------
        tab2 = ttk.Frame(notebook)
        notebook.add(tab2, text="Summary Statistics")
        ttk.Label(tab2, text="Summary Statistics per SNP-Population",
                  font=("Microsoft YaHei", 12, "bold")).pack(pady=5)

        cols2 = ("SNP", "Population", "N", "Na", "Ho", "He", "PIC",
                 "HWE_p", "Fis", "Null_Freq")
        t2_container = ttk.Frame(tab2)
        t2_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        t2 = ttk.Treeview(t2_container, columns=cols2, show="headings")
        for c in cols2:
            t2.heading(c, text=c)
            t2.column(c, width=90)
        sb2 = ttk.Scrollbar(t2_container, orient=tk.VERTICAL, command=t2.yview)
        t2.configure(yscrollcommand=sb2.set)
        t2.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb2.pack(side=tk.RIGHT, fill=tk.Y)

        pop_metrics = {}
        for snp_name, summary in locus_summaries.items():
            for pname, stats in summary.per_population.items():
                t2.insert("", tk.END, values=(
                    snp_name, pname,
                    stats.n_individuals,
                    stats.n_alleles_observed,
                    f"{stats.ho:.4f}",
                    f"{stats.he:.4f}",
                    f"{stats.pic:.4f}",
                    f"{stats.hwe_p_value:.4f}" if stats.hwe_p_value is not None else "N/A",
                    f"{stats.fis:.4f}" if stats.fis is not None else "N/A",
                    f"{stats.null_allele_freq:.4f}" if stats.null_allele_freq is not None else "N/A",
                ))
                if pname not in pop_metrics:
                    pop_metrics[pname] = {"Ho": [], "He": [], "PIC": [],
                        "HWE_p": [], "Fis": [], "Null_Freq": []}
                pop_metrics[pname]["Ho"].append(stats.ho)
                pop_metrics[pname]["He"].append(stats.he)
                pop_metrics[pname]["PIC"].append(stats.pic)
                if stats.hwe_p_value is not None:
                    pop_metrics[pname]["HWE_p"].append(stats.hwe_p_value)
                if stats.fis is not None:
                    pop_metrics[pname]["Fis"].append(stats.fis)
                if stats.null_allele_freq is not None:
                    pop_metrics[pname]["Null_Freq"].append(stats.null_allele_freq)

        if pop_metrics:
            sum_frame2 = ttk.LabelFrame(tab2, text="Average +/- SD by Population (across all SNPs)")
            sum_frame2.pack(fill=tk.X, padx=5, pady=(0, 5))
            sum_cols2 = ("Pop", "N_SNPs", "Ho", "He", "PIC", "HWE_p", "Fis", "Null_Freq")
            sum_h = min(len(pop_metrics) + 1, 7)
            t2s = ttk.Treeview(sum_frame2, columns=sum_cols2, show="headings", height=sum_h)
            t2s.heading("Pop", text="Population")
            t2s.column("Pop", width=100)
            t2s.heading("N_SNPs", text="N_SNPs")
            t2s.column("N_SNPs", width=60)
            for c in ["Ho", "He", "PIC", "HWE_p", "Fis", "Null_Freq"]:
                t2s.heading(c, text=c)
                t2s.column(c, width=125)
            t2s.pack(fill=tk.BOTH, expand=True)
            for pname in sorted(pop_metrics.keys()):
                m = pop_metrics[pname]
                def fmt_mean_sd(vals):
                    if len(vals) > 1:
                        return f"{np.mean(vals):.4f}+/-{np.std(vals, ddof=1):.4f}"
                    elif len(vals) == 1:
                        return f"{vals[0]:.4f} (n=1)"
                    else:
                        return "N/A"
                t2s.insert("", tk.END, values=(
                    pname, len(m["Ho"]),
                    fmt_mean_sd(m["Ho"]),
                    fmt_mean_sd(m["He"]),
                    fmt_mean_sd(m["PIC"]),
                    fmt_mean_sd(m["HWE_p"]),
                    fmt_mean_sd(m["Fis"]),
                    fmt_mean_sd(m["Null_Freq"]),
                ))

        def export_summary():
            path = filedialog.asksaveasfilename(defaultextension=".csv",
                filetypes=[("CSV", "*.csv")], title="Save Summary Statistics",
                initialdir=self.last_output_dir or ".")
            if not path:
                return
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(cols2)
                for snp_name, summary in locus_summaries.items():
                    for pname, stats in summary.per_population.items():
                        w.writerow([
                            snp_name, pname,
                            stats.n_individuals, stats.n_alleles_observed,
                            f"{stats.ho:.4f}", f"{stats.he:.4f}", f"{stats.pic:.4f}",
                            f"{stats.hwe_p_value:.4f}" if stats.hwe_p_value is not None else "",
                            f"{stats.fis:.4f}" if stats.fis is not None else "",
                            f"{stats.null_allele_freq:.4f}" if stats.null_allele_freq is not None else "",
                        ])
            messagebox.showinfo("Saved", f"Saved to: {path}", parent=win)

        ttk.Button(win, text="Save Summary Statistics as CSV",
                   command=export_summary).pack(in_=tab2, side=tk.BOTTOM, pady=5)

        # ---------- Tab 3: Locus-by-Locus ----------
        tab3 = ttk.Frame(notebook)
        notebook.add(tab3, text="Locus-by-Locus")
        ttk.Label(tab3, text="Locus-by-Locus F-statistics",
                  font=("Microsoft YaHei", 12, "bold")).pack(pady=5)

        cols3 = ("SNP", "Overall_Ho", "Overall_He", "Overall_PIC",
                 "N_Alleles", "Fst", "Fit", "Fst_p")
        t3_container = ttk.Frame(tab3)
        t3_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 5))
        t3 = ttk.Treeview(t3_container, columns=cols3, show="headings")
        for c in cols3:
            t3.heading(c, text=c)
            t3.column(c, width=110)
        sb3 = ttk.Scrollbar(t3_container, orient=tk.VERTICAL, command=t3.yview)
        t3.configure(yscrollcommand=sb3.set)
        t3.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb3.pack(side=tk.RIGHT, fill=tk.Y)

        locus_metric_vals = {c: [] for c in cols3[1:]}
        for snp_name, summary in locus_summaries.items():
            t3.insert("", tk.END, values=(
                snp_name,
                f"{summary.overall_ho:.4f}",
                f"{summary.overall_he:.4f}",
                f"{summary.overall_pic:.4f}",
                summary.overall_n_alleles,
                f"{summary.fst:.4f}" if summary.fst is not None else "N/A",
                f"{summary.fit:.4f}" if summary.fit is not None else "N/A",
                f"{summary.fst_p_value:.4f}" if summary.fst_p_value is not None else "N/A",
            ))
            locus_metric_vals["Overall_Ho"].append(summary.overall_ho)
            locus_metric_vals["Overall_He"].append(summary.overall_he)
            locus_metric_vals["Overall_PIC"].append(summary.overall_pic)
            locus_metric_vals["N_Alleles"].append(summary.overall_n_alleles)
            if summary.fst is not None:
                locus_metric_vals["Fst"].append(summary.fst)
            if summary.fit is not None:
                locus_metric_vals["Fit"].append(summary.fit)
            if summary.fst_p_value is not None:
                locus_metric_vals["Fst_p"].append(summary.fst_p_value)

        sum_frame3 = ttk.LabelFrame(tab3, text="Average +/- SD across all Loci")
        sum_frame3.pack(fill=tk.X, padx=5, pady=(0, 5))
        sum_cols3 = ("Measure", "Mean", "SD", "N")
        t3s = ttk.Treeview(sum_frame3, columns=sum_cols3, show="headings",
                           height=len(cols3))
        for c in sum_cols3:
            t3s.heading(c, text=c)
            t3s.column(c, width=105)
        t3s.pack(fill=tk.BOTH, expand=True)
        for metric in cols3[1:]:
            vals = locus_metric_vals[metric]
            n_total = len(locus_summaries)
            if len(vals) > 1:
                t3s.insert("", tk.END, values=(
                    metric,
                    f"{np.mean(vals):.4f}",
                    f"{np.std(vals, ddof=1):.4f}",
                    f"{len(vals)}/{n_total}",
                ))
            elif len(vals) == 1:
                t3s.insert("", tk.END, values=(
                    metric, f"{vals[0]:.4f}", "n=1", f"1/{n_total}",
                ))
            else:
                t3s.insert("", tk.END, values=(
                    metric, "N/A", "N/A", f"0/{n_total}",
                ))

        def export_locus():
            path = filedialog.asksaveasfilename(defaultextension=".csv",
                filetypes=[("CSV", "*.csv")], title="Save Locus-by-Locus F-statistics",
                initialdir=self.last_output_dir or ".")
            if not path:
                return
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(cols3)
                for snp_name, summary in locus_summaries.items():
                    w.writerow([
                        snp_name,
                        f"{summary.overall_ho:.4f}",
                        f"{summary.overall_he:.4f}",
                        f"{summary.overall_pic:.4f}",
                        summary.overall_n_alleles,
                        f"{summary.fst:.4f}" if summary.fst is not None else "",
                        f"{summary.fit:.4f}" if summary.fit is not None else "",
                        f"{summary.fst_p_value:.4f}" if summary.fst_p_value is not None else "",
                    ])
            messagebox.showinfo("Saved", f"Saved to: {path}", parent=win)

        ttk.Button(win, text="Save Locus-by-Locus as CSV",
                   command=export_locus).pack(in_=tab3, side=tk.BOTTOM, pady=5)

        # ---------- Tab 4: AMOVA ----------
        tab4 = ttk.Frame(notebook)
        notebook.add(tab4, text="AMOVA")
        ttk.Label(
            tab4,
            text="Pairwise-difference AMOVA (Excoffier et al. 1992) — Arlequin-compatible",
            font=("Microsoft YaHei", 12, "bold"),
        ).pack(pady=5)

        if amova_results:
            first_amova = amova_results[0]
            ttk.Label(
                tab4,
                text=(
                    f"Method: {first_amova.method}  |  "
                    f"Individuals: {first_amova.n_individuals}  |  "
                    f"Loci: {first_amova.n_loci}  |  "
                    f"Permutations: {first_amova.n_permutations}  |  "
                    f"Seed: {first_amova.random_seed}"
                ),
                foreground="gray",
            ).pack(pady=(0, 5))

            cols4 = ("Source", "df", "SS", "MS", "Var_Comp", "%", "Phi", "p_value")
            t4 = ttk.Treeview(tab4, columns=cols4, show="headings")
            for column in cols4:
                t4.heading(column, text=column)
                t4.column(column, width=120)
            t4.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            for result in amova_results:
                t4.insert(
                    "",
                    tk.END,
                    values=(
                        result.source,
                        result.df,
                        f"{result.ss:.4f}",
                        f"{result.ms:.4f}",
                        f"{result.variance_component:.6f}",
                        f"{result.percentage:.2f}%",
                        f"{result.phi_stat:.4f}" if result.phi_stat is not None else "",
                        f"{result.p_value:.4f}" if result.p_value is not None else "",
                    ),
                )

            def export_amova():
                path = filedialog.asksaveasfilename(
                    defaultextension=".csv",
                    filetypes=[("CSV", "*.csv")],
                    title="Save AMOVA Results",
                    initialdir=self.last_output_dir or ".",
                )
                if not path:
                    return
                with open(path, "w", newline="", encoding="utf-8-sig") as file_handle:
                    writer = csv.writer(file_handle)
                    writer.writerow(
                        cols4
                        + (
                            "Method",
                            "N_Individuals",
                            "N_Loci",
                            "N_Permutations",
                            "Random_Seed",
                        )
                    )
                    for result in amova_results:
                        writer.writerow(
                            [
                                result.source,
                                result.df,
                                result.ss,
                                result.ms,
                                result.variance_component,
                                result.percentage,
                                result.phi_stat,
                                result.p_value,
                                result.method,
                                result.n_individuals,
                                result.n_loci,
                                result.n_permutations,
                                result.random_seed,
                            ]
                        )
                messagebox.showinfo("Saved", f"Saved to: {path}", parent=win)

            ttk.Button(
                tab4, text="Save AMOVA Results as CSV", command=export_amova
            ).pack(side=tk.BOTTOM, pady=5)
        else:
            warning_text = amova_warning or "AMOVA was not calculated."
            ttk.Label(
                tab4,
                text=f"AMOVA unavailable:\n{warning_text}",
                foreground="darkorange",
                wraplength=900,
                justify=tk.LEFT,
            ).pack(pady=20, padx=20)

        # ---------- Tab 5: Genotype Freq. Test ----------
        tab5 = ttk.Frame(notebook)
        notebook.add(tab5, text="Genotype Freq. Test")
        ttk.Label(
            tab5,
            text="Genotype Frequency Differences with Multiple-testing Correction",
            font=("Microsoft YaHei", 12, "bold"),
        ).pack(pady=5)

        cols5 = (
            "SNP",
            "Chi-square",
            "df",
            "Raw_P",
            "Bonferroni_P",
            "FDR_Q",
            "Significant_FDR",
            "Note",
        )
        widths5 = {
            "SNP": 130,
            "Chi-square": 90,
            "df": 50,
            "Raw_P": 100,
            "Bonferroni_P": 110,
            "FDR_Q": 100,
            "Significant_FDR": 110,
            "Note": 320,
        }
        t5_frame = ttk.Frame(tab5)
        t5_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        t5 = ttk.Treeview(t5_frame, columns=cols5, show="headings")
        for column in cols5:
            t5.heading(column, text=column)
            t5.column(column, width=widths5[column], minwidth=50)
        sb5v = ttk.Scrollbar(t5_frame, orient=tk.VERTICAL, command=t5.yview)
        sb5h = ttk.Scrollbar(t5_frame, orient=tk.HORIZONTAL, command=t5.xview)
        t5.configure(yscrollcommand=sb5v.set, xscrollcommand=sb5h.set)
        t5.grid(row=0, column=0, sticky="nsew")
        sb5v.grid(row=0, column=1, sticky="ns")
        sb5h.grid(row=1, column=0, sticky="ew")
        t5_frame.rowconfigure(0, weight=1)
        t5_frame.columnconfigure(0, weight=1)

        for snp_name, result in gft_results.items():
            if result is None:
                values = (snp_name, "N/A", "", "", "", "", "", "Insufficient data")
            else:
                values = (
                    snp_name,
                    f"{result.chi2_stat:.4f}",
                    result.degrees_of_freedom,
                    f"{result.raw_p_value:.6g}"
                    if result.raw_p_value is not None
                    else "",
                    f"{result.bonferroni_p_value:.6g}"
                    if result.bonferroni_p_value is not None
                    else "",
                    f"{result.fdr_q_value:.6g}"
                    if result.fdr_q_value is not None
                    else "",
                    "Yes" if result.significant_fdr else "No",
                    result.warning or "",
                )
            t5.insert("", tk.END, values=values)

        def export_gft():
            path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=[("CSV", "*.csv")],
                title="Save Genotype Frequency Tests",
                initialdir=self.last_output_dir or ".",
            )
            if not path:
                return
            with open(path, "w", newline="", encoding="utf-8-sig") as file_handle:
                writer = csv.writer(file_handle)
                writer.writerow(cols5)
                for snp_name, result in gft_results.items():
                    if result is None:
                        writer.writerow(
                            [snp_name, "", "", "", "", "", "", "Insufficient data"]
                        )
                    else:
                        writer.writerow(
                            [
                                snp_name,
                                result.chi2_stat,
                                result.degrees_of_freedom,
                                result.raw_p_value,
                                result.bonferroni_p_value,
                                result.fdr_q_value,
                                result.significant_fdr,
                                result.warning or "",
                            ]
                        )
            messagebox.showinfo("Saved", f"Saved to: {path}", parent=win)

        ttk.Button(
            tab5, text="Save Genotype Frequency Tests as CSV", command=export_gft
        ).pack(side=tk.BOTTOM, pady=5)

        # ---------- Tab 6: Frequency Plots ----------
        tab6 = ttk.Frame(notebook)
        notebook.add(tab6, text="Frequency Plots")
        ttk.Label(tab6, text="Genotype Frequency Bar Plots",
                  font=("Microsoft YaHei", 12, "bold")).pack(pady=5)

        ctrl = ttk.Frame(tab6)
        ctrl.pack(fill=tk.X, padx=5, pady=5)

        valid_snps = [s for s in all_data if all_data[s]]
        ttk.Label(ctrl, text="Select SNP:").pack(side=tk.LEFT, padx=5)
        plot_snp_var = tk.StringVar()
        plot_combo = ttk.Combobox(ctrl, textvariable=plot_snp_var,
                                  values=valid_snps, state="readonly", width=30)
        if valid_snps:
            plot_combo.current(0)
            plot_combo.pack(side=tk.LEFT, padx=5)
        else:
            ttk.Label(ctrl, text="No valid SNPs", foreground="gray").pack(side=tk.LEFT)

        fig_holder = [None]
        canvas_holder = [None]

        def save_figure():
            if fig_holder[0] is None:
                messagebox.showwarning("Warning",
                    "No figure generated yet. Please select a SNP first.", parent=win)
                return
            path = filedialog.asksaveasfilename(defaultextension=".png",
                filetypes=[("PNG 300 DPI", "*.png"), ("PDF", "*.pdf")],
                title="Save Figure", initialdir=self.last_output_dir or ".")
            if not path:
                return
            try:
                fig_holder[0].savefig(path, dpi=300, bbox_inches="tight")
                messagebox.showinfo("Saved", f"Figure saved to: {path}", parent=win)
            except Exception as e:
                error_msg = str(e)
                messagebox.showerror("Error", f"Save failed: {error_msg}", parent=win)

        ttk.Button(ctrl, text="Save Figure (300 DPI)",
                   command=save_figure).pack(side=tk.LEFT, padx=10)

        plot_frame = ttk.Frame(tab6)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        def _draw_freq_plot(*_):
            loc = plot_snp_var.get()
            if not loc or loc not in all_data:
                return
            try:
                pop_names = sorted(all_data[loc].keys())
                all_gts = set()
                for p in pop_names:
                    for gt in all_data[loc][p]:
                        if gt is not None:
                            all_gts.add(gt)
                if not all_gts:
                    return
                all_gts = sorted(all_gts, key=lambda g: (str(g[0]), str(g[1])))

                counts = {p: Counter(gt for gt in all_data[loc][p] if gt)
                          for p in pop_names}
                n_pops = len(pop_names)
                n_gts = len(all_gts)
                freq_matrix = np.zeros((n_pops, n_gts))
                for i, p in enumerate(pop_names):
                    total = sum(counts[p].values())
                    if total > 0:
                        for j, gt in enumerate(all_gts):
                            freq_matrix[i, j] = counts[p].get(gt, 0) / total

                if fig_holder[0] is not None:
                    plt.close(fig_holder[0])
                fig_holder[0] = Figure(figsize=(10, 5), dpi=100)
                ax = fig_holder[0].add_subplot(111)

                x = np.arange(n_gts)
                bw = max(0.6 / n_pops, 0.15)
                offsets = np.linspace(-(n_pops - 1) * bw / 2,
                                      (n_pops - 1) * bw / 2, n_pops)
                colors = plt.cm.Set3(np.linspace(0, 1, n_pops))
                for i, p in enumerate(pop_names):
                    ax.bar(x + offsets[i], freq_matrix[i],
                           width=bw, label=p, color=colors[i], alpha=0.85)

                ax.set_xlabel("Genotype")
                ax.set_ylabel("Frequency")
                ax.set_title(f"Genotype Frequency – {loc}")
                ax.set_xticks(x)
                ax.set_xticklabels([f"{g[0]}/{g[1]}" for g in all_gts],
                                   rotation=45, ha="right", fontsize=8)
                ax.legend(title="Population")
                ax.grid(axis="y", linestyle="--", alpha=0.4)
                fig_holder[0].tight_layout()

                if canvas_holder[0] is not None:
                    canvas_holder[0].get_tk_widget().destroy()
                canvas_holder[0] = FigureCanvasTkAgg(fig_holder[0], master=plot_frame)
                canvas_holder[0].draw()
                canvas_holder[0].get_tk_widget().pack(fill=tk.BOTH, expand=True)
            except Exception as _e:
                import traceback
                traceback.print_exc()
                # Surface the error on the plot area instead of leaving the
                # tab silently blank (a blank canvas looks like a bug).
                try:
                    if fig_holder[0] is not None:
                        plt.close(fig_holder[0])
                    fig_holder[0] = Figure(figsize=(10, 5), dpi=100)
                    ax = fig_holder[0].add_subplot(111)
                    ax.axis("off")
                    ax.text(0.02, 0.95,
                            "Plot failed:\n" + traceback.format_exc()[-1500:],
                            transform=ax.transAxes, va="top", fontsize=8,
                            family="monospace", wrap=True)
                    if canvas_holder[0] is not None:
                        canvas_holder[0].get_tk_widget().destroy()
                    canvas_holder[0] = FigureCanvasTkAgg(fig_holder[0],
                                                         master=plot_frame)
                    canvas_holder[0].draw()
                    canvas_holder[0].get_tk_widget().pack(fill=tk.BOTH,
                                                          expand=True)
                except Exception:
                    traceback.print_exc()

        if valid_snps:
            plot_combo.bind("<<ComboboxSelected>>", _draw_freq_plot)
            tab6.after(300, _draw_freq_plot)

        # ---------- Tab 7: Linkage / mode comparison ---------------------
        if linkage_info or mode_comparison:
            tab7 = ttk.Frame(notebook)
            notebook.add(tab7, text="Linkage / Modes")
            ttk.Label(tab7, text="Linked SNPs and analysis modes",
                      font=("Microsoft YaHei", 12, "bold")).pack(pady=5)

            info_frame = ttk.Frame(tab7)
            info_frame.pack(fill=tk.X, padx=10)
            if mode_label:
                ttk.Label(info_frame, text=f"Mode analysed: {mode_label}",
                          font=("Microsoft YaHei", 10, "bold")).pack(anchor=tk.W)
            if linkage_info:
                n_multi = sum(1 for b in linkage_info["blocks"]
                              if len(b.loci) > 1)
                ttk.Label(info_frame, text=(
                    f"{len(linkage_info['all_loci'])} SNPs clustered into "
                    f"{len(linkage_info['blocks'])} linkage blocks "
                    f"({n_multi} with >1 SNP) at r2 >= "
                    f"{linkage_info['r2_threshold']:.2f}; "
                    f"{len(linkage_info['kept_loci'])} independent loci kept."
                    ), wraplength=1000, justify=tk.LEFT).pack(anchor=tk.W)
                if linkage_info.get("report_path"):
                    ttk.Label(info_frame,
                              text=f"LD report: {linkage_info['report_path']}",
                              foreground="#1a6fb5",
                              wraplength=1000,
                              justify=tk.LEFT).pack(anchor=tk.W)

            if mode_comparison:
                cmp_frame = ttk.LabelFrame(tab7, text="Mode comparison",
                                           padding="8")
                cmp_frame.pack(fill=tk.X, padx=10, pady=(10, 5))
                metric_keys = [
                    ("n_loci", "Loci analysed", "{:.0f}"),
                    ("mean_ho", "Mean Ho", "{:.4f}"),
                    ("mean_he", "Mean He", "{:.4f}"),
                    ("mean_pic", "Mean PIC", "{:.4f}"),
                    ("mean_fst", "Mean Fst (WC84)", "{:.4f}"),
                    ("amova_phi_st", "AMOVA Phi-ST", "{:.4f}"),
                    ("amova_p_value", "AMOVA p (permutation)", "{:.4f}"),
                ]
                header = ("Metric",) + tuple(mode_comparison)
                tree = ttk.Treeview(cmp_frame, columns=header, show="headings",
                                    height=len(metric_keys))
                for col in header:
                    tree.heading(col, text=col)
                    tree.column(col, width=230 if col == "Metric" else 200,
                                anchor=tk.CENTER)
                for key, label, fmt in metric_keys:
                    values = [label]
                    for mode in mode_comparison:
                        value = mode_comparison[mode].get(key)
                        values.append("n/a" if value is None else fmt.format(value))
                    tree.insert("", tk.END, values=values)
                tree.pack(fill=tk.X)
                ttk.Label(cmp_frame,
                          text=("Mode B keeps only SNPs that are independent "
                                "(r2 below the threshold with every other kept "
                                "SNP), so multi-locus estimates no longer "
                                "over-count the information shared inside an "
                                "amplicon."),
                          foreground="gray", wraplength=1000,
                          justify=tk.LEFT).pack(anchor=tk.W, pady=(4, 0))

            if linkage_info:
                blocks_frame = ttk.LabelFrame(
                    tab7, text="Linkage blocks and representatives",
                    padding="8")
                blocks_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
                b_cols = ("Block", "SNPs", "N", "Representative", "Reason")
                b_tree = ttk.Treeview(blocks_frame, columns=b_cols,
                                      show="headings", height=8)
                for col in b_cols:
                    b_tree.heading(col, text=col)
                b_tree.column("Block", width=80, anchor=tk.CENTER)
                b_tree.column("SNPs", width=420)
                b_tree.column("N", width=40, anchor=tk.CENTER)
                b_tree.column("Representative", width=200)
                b_tree.column("Reason", width=280)
                b_scroll = ttk.Scrollbar(blocks_frame, orient=tk.VERTICAL,
                                         command=b_tree.yview)
                b_tree.config(yscrollcommand=b_scroll.set)
                b_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
                b_scroll.pack(side=tk.RIGHT, fill=tk.Y)
                for block in linkage_info["blocks"]:
                    b_tree.insert("", tk.END, values=(
                        f"Block{block.block_id}", ", ".join(block.loci),
                        len(block.loci), block.representative or "",
                        block.representative_reason))

                pairs_frame = ttk.LabelFrame(
                    tab7, text="Pairwise LD (within populations)", padding="8")
                pairs_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
                p_cols = ("Locus1", "Locus2", "Population", "N", "r2", "D'",
                          "p (perm)")
                p_tree = ttk.Treeview(pairs_frame, columns=p_cols,
                                      show="headings", height=10)
                for col in p_cols:
                    p_tree.heading(col, text=col)
                p_tree.column("Locus1", width=180)
                p_tree.column("Locus2", width=180)
                p_tree.column("Population", width=110, anchor=tk.CENTER)
                p_tree.column("N", width=50, anchor=tk.CENTER)
                for col in ("r2", "D'", "p (perm)"):
                    p_tree.column(col, width=80, anchor=tk.CENTER)
                p_scroll = ttk.Scrollbar(pairs_frame, orient=tk.VERTICAL,
                                         command=p_tree.yview)
                p_tree.config(yscrollcommand=p_scroll.set)
                p_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
                p_scroll.pack(side=tk.RIGHT, fill=tk.Y)
                for pair in linkage_info["pairs"]:
                    def _fmt(value):
                        if value is None:
                            return ""
                        try:
                            if not np.isfinite(value):
                                return ""
                        except TypeError:
                            return ""
                        return f"{value:.4f}"
                    p_tree.insert("", tk.END, values=(
                        pair.locus1, pair.locus2, pair.population,
                        pair.n_individuals, _fmt(pair.r2), _fmt(pair.d_prime),
                        _fmt(pair.p_value)))


# Main entry point
if __name__ == "__main__":
    application = None
    try:
        root = tk.Tk()
        application = ABIToSNPApp(root)
        root.mainloop()
    finally:
        if application is not None and sys.stdout is application:
            sys.stdout = application.original_stdout

