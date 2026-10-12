# -*- mode: python ; coding: utf-8 -*-
# Slim build: PyPI (OpenBLAS) stack + aggressive exclusion of unrelated
# packages that the anaconda environment pulled into the 480 MB build.
from pathlib import Path

_base = Path('.')

# ---------- module excludes (never imported by the application) ----------
EXCLUDES = [
    # notebook / viz stacks dragged in from the anaconda env
    'panel', 'bokeh', 'botocore', 'boto3', 'llvmlite', 'numba',
    'notebook', 'jupyterlab', 'jupyter_server', 'jupyter_client',
    'jupyter_core', 'jupyter', 'nbconvert', 'nbformat', 'nbclient',
    'nbclassic', 'widgetsnbextension', 'ipykernel', 'ipywidgets',
    'ipython', 'jedi', 'parso', 'tornado', 'zmq', 'pygments',
    'plotly', 'pyarrow', 'sklearn', 'skimage', 'babel', 'grpc',
    'google.protobuf', 'lxml', 'cryptography', 'astropy',
    'astropy_iers_data', 'lief', 'jsonschema', 'prometheus_client',
    'send2trash', 'requests', 'urllib3', 'certifi', 'charset_normalizer',
    'idna', 'setuptools', 'pkg_resources', 'dask', 'distributed',
    'numexpr', 'blosc', 'tables', 'pandas.tests', 'numpy.tests',
    'numpy.distutils', 'numpy.doc', 'scipy.tests',
    'matplotlib.tests', 'matplotlib.testing', 'pywin32', 'win32com',
    'pythoncom',
]
# NOTE: numpy.f2py must NOT be excluded.  numpy's lazy __getattr__ imports
# submodules by dotted name at runtime, which PyInstaller's static analysis
# cannot see; scipy's array_api_compat clones numpy via `from numpy import *`
# and therefore needs numpy.f2py present or every scipy.stats import fails
# in the frozen app ("No module named 'numpy.f2py'").

# ---------- binary / data filters (name substrings) ----------
DROP_SUBSTR = (
    'mkl_avx512', 'mkl_mc3', 'mkl_vml_avx512', 'mkl_vml_mc3',
    'mkl_tbb_thread', 'mkl_scalapack', 'tbb12.dll',
    'bokeh', 'panel', 'notebook', 'jupyter', 'plotly', 'pyarrow',
    'arrow.dll', 'arrow_compute', 'sklearn', 'babel', 'grpc',
    'libprotobuf', 'protobuf', 'lief', '_lief', 'astropy_iers_data',
    'cryptography', '_rust.pyd', 'mfc140u', 'win32',
    'vtk.js', 'mapbox-gl', 'maplibre-gl', 'widgetbundle',
)

# Keep the full DejaVu family plus matplotlib's hard-coded fallback font.
# LastResortHE-Regular.ttf is REQUIRED: matplotlib's font fallback chain
# references it during text layout (e.g. tight_layout of rotated tick
# labels); if missing, drawing raises FileNotFoundError and the Frequency
# Plots tab renders silently blank.
FONT_KEEP = (
    'mpl-data/fonts/ttf/dejavu',
    'mpl-data/fonts/ttf/lastresorthe-regular.ttf',
)


def _low(name):
    return name.replace('\\', '/').lower()


def keep_entry(name):
    low = _low(name)
    # Tcl/Tk runtime data: keep unconditionally (e.g. ttk/notebook.tcl must
    # never be caught by the 'notebook' drop rule).
    if low.startswith('_tcl_data/') or low.startswith('_tk_data/'):
        return True
    for d in DROP_SUBSTR:
        if d in low:
            return False
    if '/tests/' in low or low.endswith('/tests') or '/test/' in low:
        return False
    if 'mpl-data/fonts/' in low:
        if not any(f in low for f in FONT_KEEP):
            return False
    if 'mpl-data/images' in low or 'mpl-data/sample_data' in low:
        return False
    return True


a = Analysis(
    [str(_base / 'abi_to_snp_analyzer_v1.1.0.py')],
    pathex=[str(_base)],
    binaries=[],
    datas=[
        (str(_base / 'population_genetics.py'), '.'),
        (str(_base / 'chromatogram.py'), '.'),
        (str(_base / 'snp_format_io.py'), '.'),
        (str(_base / 'linkage_analysis.py'), '.'),
    ],
    hiddenimports=[
        'population_genetics', 'chromatogram', 'snp_format_io',
        'linkage_analysis',
        'scipy', 'scipy.stats', 'scipy.spatial', 'scipy.spatial.distance',
        'matplotlib', 'matplotlib.backends.backend_tkagg',
        'tkinter', 'tkinter.ttk', 'tkinter.messagebox', 'tkinter.filedialog',
        'tkinter.scrolledtext',
        'PIL', 'PIL.Image', 'PIL.ImageTk',
        'Bio', 'Bio.Align',
        'openpyxl', 'pandas', 'numpy',
    ],
    excludes=EXCLUDES,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    noarchive=False,
    optimize=0,
)

# drop leftovers that Analysis still collected
a.binaries = [x for x in a.binaries if keep_entry(x[0])]
a.datas = [x for x in a.datas if keep_entry(x[0])]

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='ABI_SNP_Analyzer_v1.1.0',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    version=str(_base / 'version_info.txt'),
)
