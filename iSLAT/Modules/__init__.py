"""
Modules for iSLAT (interactive Spectral Line Analysis Tool).

This module contains core components for spectral analysis, data processing,
file handling, GUI widgets, and debugging utilities.
"""

import sys
import warnings
from importlib.util import find_spec as _find_spec

# Package metadata
from iSLAT import __version__
__author__ = "iSLAT Development Team"
__description__ = "interactive Spectral Line Analysis Tool (iSLAT) Components"

# Check Python version
if sys.version_info < (3, 8):
    raise RuntimeError("iSLAT COMPONENTS requires Python 3.8 or higher")

# Lightweight dependency checks (does NOT import the packages)
if _find_spec("numpy") is None:
    raise ImportError("NumPy is required for iSLAT COMPONENTS module")

if _find_spec("pandas") is None:
    raise ImportError("Pandas is required for iSLAT COMPONENTS module")

if _find_spec("scipy") is None:
    warnings.warn("SciPy not found. Some fitting features may be limited.")

if _find_spec("matplotlib") is None:
    warnings.warn("Matplotlib not found. Plotting features will be limited.")

if _find_spec("astroquery") is None:
    warnings.warn("Astroquery not found. HITRAN data download will be limited.")

# Define public API (lazy -- consumers import from the submodule directly)
__all__ = [
    'download_hitran_data',
    'get_Hitran_data'
]

# Configuration constants
DEFAULT_SAVE_PATH = "DATAFILES/SAVES"
DEFAULT_CONFIG_PATH = "DATAFILES/CONFIG"
DEFAULT_THEME_PATH = "DATAFILES/CONFIG/GUIThemes"

def __getattr__(name):
    """Lazy access to Hitran_data functions."""
    if name in ("download_hitran_data", "get_Hitran_data"):
        from .Hitran_data import download_hitran_data, get_Hitran_data
        globals()["download_hitran_data"] = download_hitran_data
        globals()["get_Hitran_data"] = get_Hitran_data
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")