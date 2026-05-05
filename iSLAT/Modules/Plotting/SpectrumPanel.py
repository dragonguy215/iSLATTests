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
from matplotlib.lines import Line2D

from .SpectralPanel import SpectralPanel, GapMode

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
        # Use gap-aware data when gap_mode is SKIP.
        panel_wave, panel_flux, panel_err = self.get_panel_data_with_gaps()
        self._plot_observed_spectrum(
            ax, panel_wave, panel_flux, panel_err, deduplicate=True,
        )

        # NOTE: Gap indicators are drawn later by _post_render_cell()
        # after y-limits are finalised, so that the break-mark lines
        # and text annotations are positioned correctly.

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

    # ------------------------------------------------------------------
    # Spectrum rendering helpers
    # ------------------------------------------------------------------
    def _plot_observed_spectrum(
        self,
        ax: Axes,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        error_data: Optional[np.ndarray] = None,
        color: Optional[str] = None,
        label: str = "Data",
        deduplicate: bool = False,
    ) -> None:
        """Plot observed spectrum data on *ax*.

        Parameters
        ----------
        ax : Axes
            Target axes.
        wave_data, flux_data : np.ndarray
            Observed wavelength / flux arrays.
        error_data : np.ndarray, optional
            Error bars.
        color : str, optional
            Line / marker colour.  Defaults to the theme foreground.
        label : str
            Legend label.
        deduplicate : bool
            When *True*, remove any existing artists tagged with
            ``_islat_observed`` before plotting new ones, and tag the
            newly created artists.  Useful in GUI contexts where the
            method may be called repeatedly on the same axes.
        """
        from .BasePlot import BasePlot
        if flux_data is None or len(flux_data) == 0:
            return
        if color is None:
            color = self._get_theme_value("foreground", "black")
        if deduplicate:
            BasePlot._clear_tagged_artists(ax, "_islat_observed")
        if error_data is not None and len(error_data) == len(flux_data):
            container = ax.errorbar(
                wave_data, flux_data, yerr=error_data,
                fmt="-", color=color, linewidth=1, label=label,
                zorder=self._get_theme_value("zorder_observed", 2),
                elinewidth=0.5, capsize=0,
            )
            if deduplicate:
                container[0]._islat_observed = True
                for cap in container[1]:
                    cap._islat_observed = True
                for bar_col in container[2]:
                    bar_col._islat_observed = True
        else:
            (line,) = ax.plot(
                wave_data, flux_data,
                color=color, linewidth=1, label=label,
                zorder=self._get_theme_value("zorder_observed", 2),
            )
            if deduplicate:
                line._islat_observed = True

    def _plot_summed_spectrum(
        self,
        ax: Axes,
        wave_data: np.ndarray,
        summed_flux: np.ndarray,
        color: Optional[str] = None,
        label: str = "Sum",
        deduplicate: bool = False,
    ) -> None:
        """Plot the summed model spectrum as a filled area on *ax*.

        Parameters
        ----------
        deduplicate : bool
            When *True*, remove any existing ``_islat_summed``-tagged
            collections before plotting a new fill.
        """
        from .BasePlot import BasePlot
        if summed_flux is None or len(summed_flux) == 0:
            return
        if not np.any(summed_flux > 0):
            return
        if deduplicate:
            BasePlot._clear_tagged_artists(ax, "_islat_summed", lines=False)
        fill_color = color or self._get_theme_value("summed_spectra_color", "lightgray")
        fill = ax.fill_between(
            wave_data, 0, summed_flux,
            color=fill_color, alpha=1.0, label=label,
            zorder=self._get_theme_value("zorder_summed", 1),
        )
        fill._islat_summed = True

    def _plot_molecule_spectrum(
        self,
        ax: Axes,
        molecule: "Molecule",
        wave_data: Optional[np.ndarray] = None,
        linewidth: Optional[float] = None,
        alpha: Optional[float] = None,
        linestyle: str = "--",
        interpolate_to_input: bool = False,
        target_wavelengths: Optional[np.ndarray] = None,
    ) -> Optional[Line2D]:
        """Plot a single molecule's model spectrum on *ax*."""
        if linewidth is None:
            linewidth = self._get_theme_value("model_linewidth", 0.8)
        if alpha is None:
            alpha = self._get_theme_value("model_alpha", 0.8)
        plot_lam, plot_flux = self.get_molecule_spectrum_data(
            molecule, wave_data,
            interpolate_to_input=interpolate_to_input,
            target_wavelengths=target_wavelengths,
        )
        if plot_lam is None or plot_flux is None or len(plot_flux) == 0:
            return None
        (line,) = ax.plot(
            plot_lam, plot_flux,
            linestyle=linestyle,
            color=self.get_molecule_color(molecule),
            alpha=alpha, linewidth=linewidth,
            label=self.get_molecule_display_name(molecule),
            zorder=self._get_theme_value("zorder_model", 3),
        )
        line._molecule_name = getattr(molecule, "name", "unknown")
        return line

    # ------------------------------------------------------------------
    @staticmethod
    def plot_gaussian_fit(
        ax: "Axes",
        gauss_fit: Any,
        fitted_wave: np.ndarray,
        fitted_flux: np.ndarray,
        color: str = "lime",
        linewidth: float = 2,
        zorder: int = 10,
        uncertainty_sigma: float = 3.0,
        fill_alpha: float = 0.3,
    ) -> None:
        """Plot a Gaussian fit result with uncertainty band on *ax*.

        Parameters
        ----------
        ax : Axes
            Target matplotlib Axes.
        gauss_fit : lmfit.model.ModelResult
            The fitted model result (must support ``eval_uncertainty``).
        fitted_wave, fitted_flux : np.ndarray
            Wavelength / flux arrays produced by the fit.
        color : str
            Line and fill colour.
        linewidth : float
            Width of the fit curve.
        zorder : int
            Drawing order for the fit curve.
        uncertainty_sigma : float
            Number of sigma for the uncertainty envelope.
        fill_alpha : float
            Transparency of the uncertainty band.
        """
        if gauss_fit is None or fitted_wave is None or fitted_flux is None:
            return
        ax.plot(
            fitted_wave,
            fitted_flux,
            color=color,
            linewidth=linewidth,
            zorder=zorder,
            linestyle="--",
        )
        try:
            dely = gauss_fit.eval_uncertainty(sigma=uncertainty_sigma)
            ax.fill_between(
                fitted_wave,
                fitted_flux - dely,
                fitted_flux + dely,
                color=color,
                alpha=fill_alpha,
            )
        except Exception:
            pass  # Uncertainty evaluation may fail for some fits