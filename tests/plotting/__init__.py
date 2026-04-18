# -*- coding: utf-8 -*-
"""Shared fixtures and helpers for the plotting test suite.

All matplotlib work is done with the non-interactive 'Agg' backend
so the suite is safe for headless CI environments.
"""

from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

from iSLAT.Modules.Plotting import SpectralPanel


# ======================================================================
# Synthetic-data generators
# ======================================================================

def make_wave_flux(
    wmin: float = 10.0,
    wmax: float = 20.0,
    n: int = 500,
    gap: Optional[Tuple[float, float]] = None,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate synthetic wave, flux, and error arrays.

    If *gap* is given as ``(lo, hi)`` the flux values in that range
    are set to NaN to simulate a data gap.
    """
    rng = np.random.default_rng(seed)
    wave = np.linspace(wmin, wmax, n)
    flux = 0.05 + 0.1 * np.sin(2 * np.pi * (wave - wmin) / (wmax - wmin))
    flux += rng.normal(0, 0.005, n)
    err = np.full(n, 0.005)
    if gap is not None:
        mask = (wave >= gap[0]) & (wave <= gap[1])
        flux[mask] = np.nan
    return wave, flux, err


def make_atomic_lines(waves: List[float]) -> pd.DataFrame:
    return pd.DataFrame({
        "wave": waves,
        "species": [f"Ion{i}" for i in range(len(waves))],
        "line": [f"L{i}" for i in range(len(waves))],
    })


def make_line_list(waves: List[float]) -> pd.DataFrame:
    return pd.DataFrame({
        "wave": waves,
        "species": [f"Sp{i}" for i in range(len(waves))],
        "line": [f"ML{i}" for i in range(len(waves))],
    })


# ======================================================================
# Concrete SpectralPanel stub (the ABC cannot be instantiated directly)
# ======================================================================

class ConcreteSpectralPanel(SpectralPanel):
    """Minimal concrete implementation for testing the ABC."""

    def generate_plot(self, **kwargs) -> None:
        ax = self._resolve_axes()
        wave, flux, err = self.get_panel_data()
        ax.plot(wave, flux)
        ax.set_xlim(*self.xlim)
