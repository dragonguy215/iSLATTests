"""Shared spectral resampling utilities.

Contains the flux-conserving ``spectres``-style resampling function used
by both :class:`~iSLAT.Modules.DataTypes.Molecule.Molecule` and
:class:`~iSLAT.Modules.DataTypes.Spectrum.Spectrum`.  Having a single
canonical implementation avoids duplication and makes updates easier.
"""

from __future__ import annotations

from typing import Tuple
import warnings

import numpy as np

import iSLAT.Constants as c

def make_bins(wavs: np.ndarray):
    """Given a series of wavelength points, find the edges and widths
    of corresponding wavelength bins.

    Parameters
    ----------
    wavs : np.ndarray
        1-D array of wavelength centers (must be sorted, length ≥ 2).

    Returns
    -------
    edges : np.ndarray
        Bin edges (length ``len(wavs) + 1``).
    widths : np.ndarray
        Bin widths (same length as *wavs*).
    """
    edges = np.zeros(wavs.shape[0] + 1)
    widths = np.zeros(wavs.shape[0])
    edges[0] = wavs[0] - (wavs[1] - wavs[0]) / 2
    widths[-1] = (wavs[-1] - wavs[-2])
    edges[-1] = wavs[-1] + (wavs[-1] - wavs[-2]) / 2
    edges[1:-1] = (wavs[1:] + wavs[:-1]) / 2
    widths[:-1] = edges[1:-1] - edges[:-2]
    return edges, widths


def spectres(
    new_wavs: np.ndarray,
    spec_wavs: np.ndarray,
    spec_fluxes: np.ndarray,
    fill: float = 0.0,
    verbose: bool = False,
) -> np.ndarray:
    """Flux-conserving spectral resampling onto a new wavelength basis.

    Vectorized implementation using ``np.searchsorted`` to map new bin edges
    onto the old wavelength grid in *O(n log m)* time.

    Parameters
    ----------
    new_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired.
    spec_wavs : numpy.ndarray
        1-D array containing the current wavelength sampling.
    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths in *spec_wavs*.
    fill : float, optional
        Value to use where *new_wavs* extends outside *spec_wavs* range.
    verbose : bool, optional
        If *True*, warn when fill values are used.

    Returns
    -------
    new_fluxes : numpy.ndarray
        Array of resampled flux values with same length as *new_wavs*.
    """
    old_edges, old_widths = make_bins(spec_wavs)
    new_edges, _ = make_bins(new_wavs)

    n_new = new_wavs.shape[0]
    new_fluxes = np.full(n_new, fill, dtype=np.float64)

    # Identify new bins that fall entirely within the old grid
    valid = (new_edges[:-1] >= old_edges[0]) & (new_edges[1:] <= old_edges[-1])
    if not np.any(valid):
        if verbose:
            warnings.warn(
                "spectres: new_wavs contains values outside the range "
                "in spec_wavs, new_fluxes will be filled with fill value.",
                category=RuntimeWarning,
            )
        return new_fluxes

    valid_idx = np.where(valid)[0]

    # Map new bin edges onto old bin edges via searchsorted.
    left_edges = new_edges[valid_idx]
    right_edges = new_edges[valid_idx + 1]

    start_idx = np.searchsorted(old_edges, left_edges, side="right") - 1
    stop_idx = np.searchsorted(old_edges, right_edges, side="left") - 1

    # Clip to valid bin range
    n_old = spec_wavs.shape[0]
    np.clip(start_idx, 0, n_old - 1, out=start_idx)
    np.clip(stop_idx, 0, n_old - 1, out=stop_idx)

    # Fast path: bins where start == stop (new bin is entirely inside one old bin)
    same = start_idx == stop_idx
    if np.any(same):
        new_fluxes[valid_idx[same]] = spec_fluxes[start_idx[same]]

    # Handle bins that span multiple old bins
    diff_mask = ~same
    if np.any(diff_mask):
        diff_idx = valid_idx[diff_mask]
        d_start = start_idx[diff_mask]
        d_stop = stop_idx[diff_mask]
        d_left = left_edges[diff_mask]
        d_right = right_edges[diff_mask]

        start_factor = (
            (old_edges[d_start + 1] - d_left)
            / (old_edges[d_start + 1] - old_edges[d_start])
        )
        end_factor = (
            (d_right - old_edges[d_stop])
            / (old_edges[d_stop + 1] - old_edges[d_stop])
        )

        for k in range(len(diff_idx)):
            s = d_start[k]
            e = d_stop[k]
            sl = slice(s, e + 1)
            w = old_widths[sl].copy()
            w[0] *= start_factor[k]
            w[-1] *= end_factor[k]
            fw = w * spec_fluxes[sl]
            new_fluxes[diff_idx[k]] = fw.sum() / w.sum()

    return new_fluxes

def flux_integral(lam, flux, err, lam_min, lam_max) -> Tuple[float, float]:
    wavelength_mask = (lam >= lam_min) & (lam <= lam_max)

    if not np.any (wavelength_mask):
        return 0.0, 0.0

    lam_range = lam[wavelength_mask]
    flux_range = flux[wavelength_mask]

    if len (lam_range) < 2:
        return 0.0, 0.0

    # Convert to frequency space for proper integration
    freq_range = c.SPEED_OF_LIGHT_MICRONS / lam_range[::-1]

    # Integrate in frequency space (reverse order for proper frequency ordering)
    line_flux_meas = np.trapezoid(flux_range[::-1], x=freq_range[::-1])
    line_flux_meas = -line_flux_meas * 1e-23  # Convert Jy*Hz to erg/s/cm^2

    # Calculate error propagation if error data provided
    if err is not None:
        err_range = err[wavelength_mask]
        line_err_meas = np.trapezoid(err_range[::-1], x=freq_range[::-1])
        line_err_meas = -line_err_meas * 1e-23
    else:
        line_err_meas = 0.0

    return line_flux_meas, line_err_meas