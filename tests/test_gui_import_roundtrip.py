# -*- coding: utf-8 -*-
"""End-to-end check of the "Export SNP Formats" -> "Load Existing Data" path.

The GUI class is exercised without creating a Tk window: the heavy optional
dependencies (matplotlib, Biopython, Pillow, scipy) are stubbed, so this test
only needs tkinter + pandas.  It asserts that a file written by the export
dialog is re-imported into the same genotype table and can be handed straight
to ``build_population_analysis_inputs``.
"""

import importlib.util
import os
import sys
import tempfile
import types

import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN_SCRIPT = os.path.join(PROJECT_ROOT, "abi_to_snp_analyzer_v1.1.0.py")


class _StubModule(types.ModuleType):
    """Anything-goes stub used to satisfy optional imports."""

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        module = _StubModule(name)
        sys.modules.setdefault(self.__name__ + "." + name, module)
        return module

    def __call__(self, *args, **kwargs):
        return None


def _load_app():
    """Import the GUI module with optional dependencies stubbed.

    The stubs are removed from ``sys.modules`` afterwards so that they do not
    leak into unrelated tests (a stale ``scipy`` stub makes every later
    chi-square test silently fall back to its no-scipy branch).
    """
    import tempfile
    os.environ.setdefault(
        "MPLCONFIGDIR",
        os.path.join(tempfile.gettempdir(), "abi_snp_mplconfig"))
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    added = {}
    for name in ("Bio", "Bio.Align", "Bio.SeqIO", "matplotlib",
                 "matplotlib.pyplot", "matplotlib.figure",
                 "matplotlib.backends", "matplotlib.backends.backend_tkagg",
                 "matplotlib.backends.backend_agg", "PIL", "PIL.Image",
                 "PIL.ImageTk", "scipy", "scipy.stats"):
        if name not in sys.modules:
            added[name] = _StubModule(name)
            sys.modules[name] = added[name]
    try:
        if PROJECT_ROOT not in sys.path:
            sys.path.insert(0, PROJECT_ROOT)
        spec = importlib.util.spec_from_file_location("abi_snp_app", MAIN_SCRIPT)
        app = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(app)
    finally:
        # Remove every stub left behind, including nested ones that the stub
        # package's __getattr__ auto-created (e.g. "Bio.AlignIO").  A leftover
        # stub for a submodule breaks the real import later in the session.
        for name, module in list(sys.modules.items()):
            if isinstance(module, _StubModule):
                del sys.modules[name]
    return app


APP = _load_app()
GUI = next(obj for _, obj in vars(APP).items()
           if isinstance(obj, type) and "_import_genotype_file" in obj.__dict__)

LOCI = ["SNP01_A>C", "SNP02", "SNP03_G>T"]
SAMPLES = ["A1", "A2", "B1", "B2", "C1", "C2"]
_GENOTYPES = ["AA", "AT", "TT", "CC", "GG", "--"]
CALLS = {locus: {sample: _GENOTYPES[(i + j) % 6]
                 for i, sample in enumerate(SAMPLES)}
         for j, locus in enumerate(LOCI)}
POPULATIONS = {"PopA": ["A1", "A2"], "PopB": ["B1", "B2"], "PopC": ["C1", "C2"]}

EXPORTERS = [
    ("genepop", "_write_genepop_snp", "g.txt"),
    ("arlequin", "_write_arlequin_snp", "g.arp"),
    ("hierfstat", "_write_hierfstat_snp", "g.hierfstat.txt"),
    ("genalex", "_write_genalex_snp", "g.csv"),
    ("fstat", "_write_fstat_snp", "g.dat"),
    ("popgene", "_write_popgene_snp", "g.pop"),
    ("structure", "_write_structure_snp", "g.str.txt"),
]

# STRUCTURE has no header line, so locus names cannot be carried.
LOCUS_NAME_LOSSY = {"structure"}


class _DummyButton(object):
    def __init__(self):
        self.state = "disabled"

    def config(self, **kwargs):
        self.state = kwargs.get("state", self.state)


def _make_gui():
    gui = object.__new__(GUI)
    gui.sample_names = list(SAMPLES)
    gui.snp_names = list(LOCI)
    gui.genotype_matrix_data = {k: dict(v) for k, v in CALLS.items()}
    gui.populations = {k: set(v) for k, v in POPULATIONS.items()}
    gui.last_output_dir = None
    gui.pop_manager_button = _DummyButton()
    gui.export_formats_button = _DummyButton()
    gui.stats_button = _DummyButton()
    gui.previews = []
    gui._show_data_preview = lambda filename, **kw: gui.previews.append((filename, kw))
    return gui


@pytest.mark.parametrize("fmt,method,filename", EXPORTERS)
def test_export_then_load_feeds_population_analysis(fmt, method, filename):
    gui = _make_gui()
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, filename)
        assert getattr(gui, method)(path, True) == len(SAMPLES)

        fresh = _make_gui()
        fresh.sample_names, fresh.snp_names = [], []
        fresh.genotype_matrix_data, fresh.populations = {}, {}
        fresh._import_genotype_file(path)

        assert fresh.previews, "data preview was not shown"
        assert len(fresh.sample_names) == len(SAMPLES)
        if fmt not in LOCUS_NAME_LOSSY:
            assert fresh.snp_names == LOCI
        assert len(fresh.populations) == len(POPULATIONS)

        # The imported table must be usable by the statistical module.
        all_data, _geno, _map, _un = APP.build_population_analysis_inputs(
            fresh.genotype_matrix_data, fresh.snp_names, fresh.populations)
        assert len(all_data) == len(LOCI)
        assert all(len(pops) >= 2 for pops in all_data.values())

        # >= 2 populations -> "Statistical Analysis" is available immediately.
        assert fresh.stats_button.state == "normal"


@pytest.mark.parametrize("fmt", ["genepop", "arlequin", "genalex"])
def test_export_then_load_keeps_genotypes_and_ids(fmt):
    """The three formats that carry individual IDs must round-trip exactly."""
    method, filename = next((m, f) for k, m, f in EXPORTERS if k == fmt)
    gui = _make_gui()
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, filename)
        getattr(gui, method)(path, True)

        fresh = _make_gui()
        fresh._import_genotype_file(path)
        assert set(fresh.sample_names) == set(SAMPLES)
        for locus in LOCI:
            for sample in SAMPLES:
                assert fresh.genotype_matrix_data[locus][sample] == CALLS[locus][sample]


# ---------------------------------------------------------------------------
# Linked-SNP dual mode: one representative per linkage block vs all SNPs
# ---------------------------------------------------------------------------

import numpy as np  # noqa: E402

_DUAL_LOCI = ["SNP1", "SNP2", "SNP3", "SNP4"]
_DUAL_POPS = {"PopA": 0.9, "PopB": 0.5, "PopC": 0.1}   # minor-allele freqs


def _build_dual_mode_matrix(seed=5):
    """SNP1/SNP2 share an amplicon (perfect LD); SNP3/SNP4 are independent."""
    rng = np.random.default_rng(seed)
    matrix = {locus: {} for locus in _DUAL_LOCI}
    populations = {}
    index = 0
    for pop, minor_freq in _DUAL_POPS.items():
        members = []
        for _ in range(20):
            sample = f"{pop}_S{index:02d}"
            index += 1
            # SNP1/SNP2: haplotypes A-C (major/major) and T-G (minor/minor)
            dosage = int(rng.choice([0, 1, 2],
                                    p=[(1 - minor_freq) ** 2,
                                       2 * minor_freq * (1 - minor_freq),
                                       minor_freq ** 2]))
            matrix["SNP1"][sample] = ["AA", "AT", "TT"][dosage]
            matrix["SNP2"][sample] = ["CC", "CG", "GG"][dosage]
            for locus in ("SNP3", "SNP4"):
                d = int(rng.choice([0, 1, 2],
                                   p=[(1 - minor_freq) ** 2,
                                      2 * minor_freq * (1 - minor_freq),
                                      minor_freq ** 2]))
                matrix[locus][sample] = {"SNP3": ["GG", "GT", "TT"],
                                         "SNP4": ["CC", "CA", "AA"]}[locus][d]
            members.append(sample)
        populations[pop] = set(members)
    return matrix, populations


def test_dual_mode_representative_subset_is_smaller_and_analysable():
    matrix, populations = _build_dual_mode_matrix()
    all_data, sample_genotypes, sample_to_population, _un = \
        APP.build_population_analysis_inputs(matrix, _DUAL_LOCI, populations)

    assert APP.LINKAGE_AVAILABLE
    kept, blocks, pairs, summaries = APP.linkage.select_independent_loci(
        sample_genotypes, _DUAL_LOCI, sample_to_population,
        r2_threshold=0.8, n_permutations=99)

    assert len(kept) == 3, f"expected 3 independent loci, got {kept}"
    assert "SNP3" in kept and "SNP4" in kept
    assert kept.count("SNP1") + kept.count("SNP2") == 1
    linked_block = next(b for b in blocks if len(b.loci) > 1)
    assert set(linked_block.loci) == {"SNP1", "SNP2"}

    # the linked pair must be flagged as such, the independent pairs must not
    linked = next(s for s in summaries
                  if {s.locus1, s.locus2} == {"SNP1", "SNP2"})
    assert linked.r2_max > 0.9
    for summary in summaries:
        if {summary.locus1, summary.locus2} == {"SNP1", "SNP2"}:
            continue
        assert summary.r2_max is None or summary.r2_max < 0.5

    # both modes must produce analysable summaries
    summaries_all = {}
    for locus in _DUAL_LOCI:
        summaries_all[locus] = APP.PopulationGenetics.analyze_locus(
            all_data[locus], locus, n_permutations=99, random_seed=1)
    summaries_kept = {}
    for locus in kept:
        summaries_kept[locus] = APP.PopulationGenetics.analyze_locus(
            all_data[locus], locus, n_permutations=99, random_seed=1)

    mode_all = GUI._mode_summary((all_data, summaries_all, [], [], None))
    mode_kept = GUI._mode_summary((all_data, summaries_kept, [], [], None))
    assert mode_all["n_loci"] == 4
    assert mode_kept["n_loci"] == 3
    assert mode_all["mean_fst"] is not None
    assert mode_kept["mean_fst"] is not None
