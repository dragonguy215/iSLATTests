"""
SpectrumPanel -- Concrete :class:`SpectralPanel` for observed spectrum display.

Renders observed data, individual molecule models, a summed model
spectrum, and optionally line-list annotations and atomic-line markers
onto a single axes.  This is the panel type used by
:class:`FullSpectrumPlot` (one per row) and the upper sub-panel in each
:class:`ResidualSpectrumPlot` row.
"""

from typing import Optional, Tuple, List, Dict, Any, TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from .SpectralPanel import SpectralPanel

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict


class SpectrumPanel(SpectralPanel):
    """
    Single-axes panel showing observed spectrum + molecule overlays.

    Parameters
    ----------
    wave_data : np.ndarray
        Full observed wavelength array.
    flux_data : np.ndarray
        Full observed flux array.
    xmin, xmax : float
        Wavelength bounds for this panel.
    error_data : np.ndarray, optional
        Flux uncertainties.
    molecules : MoleculeDict, optional
        Molecule collection -- visible ones are plotted.
    mol_cache : list[tuple], optional
        Pre-computed molecule cache from
        :meth:`FullSpectrumPlot._build_mol_cache`.  When provided the
        panel slices into this cache instead of recomputing molecule
        spectra itself.
    summed_wave : np.ndarray, optional
        Wavelength array for the summed model spectrum.
    summed_flux : np.ndarray, optional
        Flux array for the summed model spectrum.
    line_list : pd.DataFrame, optional
        Saved-line annotations.
    atomic_lines : pd.DataFrame, optional
        Atomic-line annotations.
    wave_data_obs : np.ndarray, optional
        Observer-frame wavelengths (for matched spectral sampling).
    ax : Axes, optional
        Pre-existing axes to draw into.
    **kwargs
        Forwarded to :class:`SpectralPanel` / :class:`BasePlot`.
    """

    def __init__(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        xmin: float,
        xmax: float,
        error_data: Optional[np.ndarray] = None,
        molecules: Optional["MoleculeDict"] = None,
        mol_cache: Optional[List[tuple]] = None,
        summed_wave: Optional[np.ndarray] = None,
        summed_flux: Optional[np.ndarray] = None,
        line_list: Optional[pd.DataFrame] = None,
        atomic_lines: Optional[pd.DataFrame] = None,
        wave_data_obs: Optional[np.ndarray] = None,
        ax: Optional[Axes] = None,
        **kwargs,
    ):
        super().__init__(
            wave_data=wave_data,
            flux_data=flux_data,
            xmin=xmin,
            xmax=xmax,
            error_data=error_data,
            molecules=molecules,
            ax=ax,
            **kwargs,
        )
        self.mol_cache: List[tuple] = mol_cache or []
        self.summed_wave = summed_wave
        self.summed_flux = summed_flux
        self.line_list = line_list
        self.atomic_lines = atomic_lines
        self.wave_data_obs: Optional[np.ndarray] = (
            np.asarray(wave_data_obs) if wave_data_obs is not None else None
        )

    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Render the observed spectrum, molecule overlays, and annotations."""
        ax = self._resolve_axes()
        xr = self.xlim

        # -- Observed spectrum ------------------------------------------
        panel_wave, panel_flux, panel_err = self.get_panel_data()
        self._plot_observed_spectrum(
            ax, panel_wave, panel_flux, panel_err, deduplicate=True,
        )

        # -- Molecule models (slice pre-computed cache) -----------------
        for m_lam, m_flux, m_color, m_label, m_name in self.mol_cache:
            m_mask = (m_lam >= xr[0]) & (m_lam <= xr[1])
            if np.any(m_mask):
                (line,) = ax.plot(
                    m_lam[m_mask],
                    m_flux[m_mask],
                    linestyle="--",
                    color=m_color,
                    alpha=self._get_theme_value(
                        "full_spectrum_model_alpha", 0.8
                    ),
                    linewidth=self._get_theme_value(
                        "full_spectrum_model_linewidth", 0.8
                    ),
                    label=m_label,
                    zorder=self._get_theme_value("zorder_model", 3),
                )
                line._molecule_name = m_name

        # -- Summed spectrum fill ---------------------------------------
        if self.summed_wave is not None and self.summed_flux is not None:
            s_mask = (self.summed_wave >= xr[0]) & (self.summed_wave <= xr[1])
            if np.any(s_mask):
                self._plot_summed_spectrum(
                    ax, self.summed_wave[s_mask], self.summed_flux[s_mask],
                )

        # -- Axes configuration -----------------------------------------
        ax.set_xlim(*xr)

        # NOTE: Line annotations and atomic lines are NOT drawn here.
        # They are deferred to _post_render_cell() (or explicit calls to
        # plot_saved_lines / plot_atomic_lines) so that the labels are
        # positioned relative to the *final* y-limits set by the parent
        # stacked-panel composer.
