# -*- coding: utf-8 -*-
"""Load the main application module (dotted file name) for unit tests.

The main script is named ``abi_to_snp_analyzer_v1.1.0.py``; its dotted name
cannot be imported with a plain ``import`` statement, so it is loaded by file
location.  Older test files referenced ``abi_to_snp_analyzer_v1_1``, which no
longer exists.
"""

import importlib.util
import os
import sys
import tempfile

# Matplotlib builds a font cache in ~/.matplotlib on first import, which is
# not writable in restricted environments; point the cache at a temp dir
# BEFORE the application module imports matplotlib.
os.environ.setdefault(
    "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "abi_snp_mplconfig"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_TESTS_DIR)
_MAIN_SCRIPT = os.path.join(_PROJECT_ROOT, "abi_to_snp_analyzer_v1.1.0.py")
_MODULE_NAME = "abi_snp_app_full"


def load_app():
    """Import and cache the main application module (with real dependencies)."""
    module = sys.modules.get(_MODULE_NAME)
    if module is not None:
        return module
    if _PROJECT_ROOT not in sys.path:
        sys.path.insert(0, _PROJECT_ROOT)
    spec = importlib.util.spec_from_file_location(_MODULE_NAME, _MAIN_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[_MODULE_NAME] = module
    spec.loader.exec_module(module)
    return module
