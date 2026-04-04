"""
ResidualSpectrumPlot -- Multi-panel spectrum overview with per-panel residuals.

Extends :class:`FullSpectrumPlot` by adding a residual sub-panel below each
wavelength-range row, showing ``(observed - model)`` with error bars.
Optionally annotates each row with its reduced chi^2 and displays global
goodness-of-fit statistics (chi^2, degrees of freedom, reduced chi^2) at the
bottom of the figure.

Designed for best-fit comparison plots, works with any model flux array.
"""

from typing import Optional, Tuple, List, Dict, Any, Union, TYPE_CHECKING
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec

from .FullSpectrumPlot import FullSpectrumPlot
from .BasePlot import BasePlot

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

class ResidualSpectrumPlot(FullSpectrumPlot):
    """
    Multi-panel spectrum + residual plot.

    Inherits the observed-spectrum, molecule-overlay, summed-model, and
    annotation rendering from :class:`FullSpectrumPlot` and adds:

    * A **residual sub-panel** below each spectrum panel showing
      ``(data - model)`` with error bars.
    * **Per-panel** and **global** reduced chi^2 statistics.
    * Support for **nuisance parameters** (noise floor, linear continuum).
    * **Excluded wavelength ranges** and **atomic-line exclusion windows**.
    * Explicit **model component** overlays (as an alternative to
      ``MoleculeDict``).

    Parameters
    ----------
    wave_data : np.ndarray
        Observed wavelength array.
    flux_data : np.ndarray
        Observed flux array.
    model_flux : np.ndarray
        Model flux array on the same grid as *wave_data* / *flux_data*.
    error_data : np.ndarray, optional
        Flux uncertainties (used for residual error bars and chi^2).
    model_components : list[dict], optional
        Individual model components for overlay.  Each dict should
        contain ``"wave"``, ``"flux"``, ``"color"``, and ``"label"``
        keys.
    molecules : MoleculeDict, optional
        If provided (instead of *model_components*), visible molecules
        are plotted as in :class:`FullSpectrumPlot`.
    model_wave : np.ndarray, optional
        Wavelength grid of the combined model (for the filled-area
        "summed" display).  Defaults to *wave_data*.
    model_flux_hires : np.ndarray, optional
        High-resolution combined model flux corresponding to
        *model_wave*.  When given, the filled area uses this array
        while residuals still use *model_flux* (which is on the data
        grid).
    line_list : pd.DataFrame, optional
        Saved-line annotations.
    atomic_lines : pd.DataFrame, optional
        Atomic-line annotations.
    n_panels : int
        Target number of wavelength rows.  Default 5.
    step : float, optional
        Fixed wavelength width per row (overrides *n_panels*).
    xlim_range : tuple[float, float], optional
        Wavelength sub-range to display.
    ymax_factor : float
        Fractional headroom above the flux peak per panel (0.2 = 20 %).
    residual_height_ratio : float
        Height of the residual sub-panel relative to the spectrum
        sub-panel.  Default 0.3 (residual is 30 % of spectrum height).
    uniform_ylim_spectra : bool
        When *True* all spectrum sub-panels share the same vertical
        scale (global min/max across panels).  Default *False*.
    uniform_ylim_residuals : bool
        When *True* all residual sub-panels share the same vertical
        scale (the worst-case symmetric range across panels).
        Default *True*.
    show_chi2_per_panel : bool
        Annotate each row with its per-panel reduced chi^2.  Default True.
    show_total_chi2 : bool
        Show global chi^2, dof, and reduced chi^2 at the bottom of the
        figure.  Default True.
    n_free_params : int
        Number of free model parameters (used for dof calculation).
        Default 0.
    excluded_ranges : list[tuple[float, float]], optional
        Wavelength ranges excluded from the fit.  These are shaded
        in the plot and excluded from the chi^2 calculation.
    exclude_lines_half_width : float, optional
        Half-width (in the same units as *wave_data*) of exclusion
        windows centred on each line in *atomic_lines*.  When set
        together with *atomic_lines*, data points within +/- this
        distance of any atomic line are excluded from chi^2 statistics
        (and shaded on the residual panels).  Default ``None``
        (no extra exclusion beyond *excluded_ranges*).
    noise_floor : float, optional
        Systematic noise floor *s* (Jy).  When given, the effective
        per-pixel uncertainty becomes ``sigma_eff = sqrt(sigma^2 + s^2)`` for
        residual error bars and adjusted chi^2 statistics.  Both the
        raw and adjusted chi^2 are shown.  Default ``None``.
    continuum_c0 : float, optional
        Continuum offset constant (Jy).  When given (together with
        *continuum_c1* and *lam_ref*), the additive continuum
        ``c_0 + c_1*(lam - lam_ref)`` is overlaid on spectrum panels and
        added to the model for residual / chi^2 computation.
    continuum_c1 : float, optional
        Continuum slope (Jy / um).  Default ``None``.
    lam_ref : float, optional
        Reference wavelength for the continuum model (um).
    figsize : tuple, optional
        Figure size.  Auto-scaled if *None*.
    """

    def __init__(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        model_flux: np.ndarray,
        error_data: Optional[np.ndarray] = None,
        model_components: Optional[List[Dict[str, Any]]] = None,
        molecules: Optional["MoleculeDict"] = None,
        model_wave: Optional[np.ndarray] = None,
        model_flux_hires: Optional[np.ndarray] = None,
        line_list: Optional[pd.DataFrame] = None,
        atomic_lines: Optional[pd.DataFrame] = None,
        n_panels: int = 5,
        step: Optional[float] = None,
        xlim_range: Optional[Tuple[float, float]] = None,
        ymax_factor: float = 0.2,
        residual_height_ratio: float = 0.3,
        uniform_ylim_spectra: bool = False,
        uniform_ylim_residuals: bool = True,
        show_chi2_per_panel: bool = True,
        show_total_chi2: bool = True,
        n_free_params: int = 0,
        excluded_ranges: Optional[List[Tuple[float, float]]] = None,
        exclude_lines_half_width: Optional[float] = None,
        noise_floor: Optional[float] = None,
        continuum_c0: Optional[float] = None,
        continuum_c1: Optional[float] = None,
        lam_ref: Optional[float] = None,
        **kwargs,
    ):
        # -- Initialise the FullSpectrumPlot base ---------------------
        # Pass figsize=None so the parent stores None; we override
        # self._figsize below to account for the residual sub-panels.
        caller_figsize = kwargs.pop("figsize", None)
        super().__init__(
            wave_data=wave_data,
            flux_data=flux_data,
            molecules=molecules,
            error_data=error_data,
            line_list=line_list,
            atomic_lines=atomic_lines,
            n_panels=n_panels,
            step=step,
            xlim_range=xlim_range,
            ymax_factor=ymax_factor,
            uniform_ylim=uniform_ylim_spectra,
            figsize=caller_figsize,
            **kwargs,
        )

        # Uniform y-limit flags (spectra flag stored on parent as
        # self.uniform_ylim; residuals flag is RSP-specific)
        self.uniform_ylim_residuals = uniform_ylim_residuals

        # -- Model arrays (unique to ResidualSpectrumPlot) ------------
        self.model_flux = np.asarray(model_flux)
        self.model_components = model_components
        self.model_wave = (
            np.asarray(model_wave) if model_wave is not None else self.wave_data
        )
        self.model_flux_hires = (
            np.asarray(model_flux_hires)
            if model_flux_hires is not None
            else None
        )

        # -- Residual / chi^2 settings --------------------------------
        self.residual_height_ratio = residual_height_ratio
        self.show_chi2_per_panel = show_chi2_per_panel
        self.show_total_chi2 = show_total_chi2
        self.n_free_params = n_free_params
        self.excluded_ranges = excluded_ranges or []
        self.exclude_lines_half_width = exclude_lines_half_width

        # -- Nuisance parameters --------------------------------------
        self.noise_floor = noise_floor
        self.continuum_c0 = continuum_c0
        self.continuum_c1 = continuum_c1
        self.lam_ref = lam_ref
        self._has_continuum = (
            continuum_c0 is not None
            and continuum_c1 is not None
            and lam_ref is not None
        )
        self._has_noise_floor = noise_floor is not None and noise_floor > 0.0

        # -- Pre-compute adjusted arrays ------------------------------
        if self._has_continuum:
            cont = continuum_c0 + continuum_c1 * (self.wave_data - lam_ref)
            self._model_flux_adj = self.model_flux + cont
            if self.model_flux_hires is not None:
                cont_hr = continuum_c0 + continuum_c1 * (self.model_wave - lam_ref)
                self._model_flux_hires_adj = self.model_flux_hires + cont_hr
            else:
                self._model_flux_hires_adj = None
        else:
            self._model_flux_adj = self.model_flux
            self._model_flux_hires_adj = self.model_flux_hires

        if self._has_noise_floor and self.error_data is not None:
            self._error_adj = np.sqrt(self.error_data ** 2 + noise_floor ** 2)
        else:
            self._error_adj = self.error_data

        # -- Override figsize for the taller residual layout ----------
        if caller_figsize is None:
            n = len(self._panel_edges)
            row_height = 2.0 + 2.0 * self.residual_height_ratio
            self._figsize = (
                14,
                row_height * n + (1.0 if self.show_total_chi2 else 0.0),
            )

        # Subplot storage override: {panel_idx: (spectrum_ax, residual_ax)}
        self.subplots: Dict[int, Tuple[Axes, Axes]] = {}

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _is_excluded(self, wave_val: float) -> bool:
        """Return True if a wavelength falls in an excluded range."""
        for lo, hi in self.excluded_ranges:
            if lo <= wave_val <= hi:
                return True
        return False

    def _get_fit_mask(self, wave_arr: np.ndarray) -> np.ndarray:
        """Boolean mask: True for data points NOT in any excluded range
        and NOT within *exclude_lines_half_width* of any atomic line."""
        mask = np.ones(len(wave_arr), dtype=bool)
        for lo, hi in self.excluded_ranges:
            mask &= ~((wave_arr >= lo) & (wave_arr <= hi))
        # Also exclude points near atomic / ionic lines
        if (
            self.exclude_lines_half_width is not None
            and self.atomic_lines is not None
            and len(self.atomic_lines) > 0
        ):
            hw = self.exclude_lines_half_width
            for wl in self.atomic_lines["wave"].values:
                mask &= ~(np.abs(wave_arr - wl) <= hw)
        return mask

    def _compute_chi2(
        self,
        flux: np.ndarray,
        model_raw: np.ndarray,
        model_adj: np.ndarray,
        error_raw: Optional[np.ndarray],
        error_adj: Optional[np.ndarray],
        fit_mask: np.ndarray,
    ) -> Tuple[float, float, int]:
        """Compute raw and adjusted chi^2 over the given fit mask.

        Parameters
        ----------
        flux : np.ndarray
            Observed flux.
        model_raw : np.ndarray
            Raw model flux (before nuisance continuum).
        model_adj : np.ndarray
            Adjusted model flux (with continuum added).
        error_raw : np.ndarray or None
            Original pipeline errors.
        error_adj : np.ndarray or None
            Effective errors (with noise floor added in quadrature).
        fit_mask : np.ndarray[bool]
            Boolean mask selecting points to include.

        Returns
        -------
        chi2_raw : float
            chi^2 computed with *model_raw* and *error_raw*.
        chi2_adj : float
            chi^2 computed with *model_adj* and *error_adj*.
        n_fit : int
            Number of True entries in *fit_mask*.
        """
        n_fit = int(np.sum(fit_mask))
        if n_fit == 0 or error_adj is None:
            return 0.0, 0.0, n_fit

        chi2_raw = (
            float(np.sum(
                ((flux[fit_mask] - model_raw[fit_mask])
                 / error_raw[fit_mask]) ** 2
            ))
            if error_raw is not None
            else 0.0
        )
        chi2_adj = float(np.sum(
            ((flux[fit_mask] - model_adj[fit_mask])
             / error_adj[fit_mask]) ** 2
        ))
        return chi2_raw, chi2_adj, n_fit

    def _format_chi2_annotation(
        self,
        chi2_raw: float,
        chi2_adj: float,
        dof: int,
        has_nuisance: bool,
    ) -> str:
        """Build a LaTeX annotation string for reduced chi^2 values."""
        if has_nuisance:
            r_raw = chi2_raw / dof
            r_adj = chi2_adj / dof
            return (
                f"$\\chi^2_\\nu = {r_raw:.2f}$\n"
                f"$\\chi^2_{{\\nu,\\mathrm{{adj}}}} = {r_adj:.2f}$"
            )
        return f"$\\chi^2_\\nu = {chi2_adj / dof:.2f}$"

    # ------------------------------------------------------------------
    def _residual_ylim_fn(self, mask: np.ndarray) -> Tuple[float, float]:
        """Per-panel y-limit function for residual sub-panels.

        Computes a symmetric ``(-pad, +pad)`` range based on the
        maximum absolute residual and (if available) the maximum
        adjusted error within the panel.

        Parameters
        ----------
        mask : np.ndarray
            Boolean mask into :attr:`wave_data` selecting the points
            that fall within the current panel.

        Returns
        -------
        tuple[float, float]
            ``(-pad, +pad)`` for the residual panel.
        """
        panel_flux = self.flux_data[mask]
        panel_model = self._model_flux_adj[mask]
        panel_err = (
            self._error_adj[mask] if self._error_adj is not None else None
        )
        if len(panel_flux) > 0:
            residuals = panel_flux - panel_model
            res_abs_max = float(np.nanmax(np.abs(residuals)))
            if panel_err is not None and len(panel_err) > 0:
                err_max = float(np.nanmax(panel_err))
                res_abs_max = max(res_abs_max, err_max)
            res_pad = res_abs_max * 1.3 if res_abs_max > 0 else 0.01
            return (-res_pad, res_pad)
        return (-0.01, 0.01)

    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Build the multi-panel figure with residual rows."""
        n = len(self._panel_edges)
        self._ensure_figure()
        self.fig.clf()
        self.subplots.clear()

        # Disable constrained_layout so every panel pair has exactly
        # the same height -- the engine would otherwise redistribute
        # vertical space per-row based on tick labels and legends.
        # Use None (not "none") to fully clear the engine -- "none"
        # installs a PlaceHolderLayoutEngine that still triggers
        # warnings from subplots_adjust().
        self.fig.set_layout_engine(None)

        fg = self._get_theme_value("foreground", "black")

        # --- Nested GridSpec layout ------------------------------------
        # Outer grid: n rows (one per wavelength range), matching
        # FullSpectrumPlot's panel positions so that the TOP edge of
        # each spectrum sub-panel lines up exactly.
        # Inner grid (per row): 2 sub-rows [spectrum, residual] with
        # hspace=0 so they are flush against each other.
        from matplotlib.gridspec import GridSpecFromSubplotSpec

        gs_outer = GridSpec(
            nrows=n, ncols=1,
            figure=self.fig,
            hspace=0.15,
        )

        # Fixed margins: top leaves room for the molecule legend,
        # bottom leaves room for the global chi-squared box.
        self.fig.subplots_adjust(
            left=0.06, right=0.92, top=0.93, bottom=0.06,
        )

        # --- Pre-compute molecule data via inherited helper ------------
        mol_cache, mol_labels, mol_colors = self._build_mol_cache()

        # --- Pre-compute per-panel y-limits (pass 1) ------------------
        # Spectrum y-limits -- delegated to the inherited helper
        spec_ylims = self._compute_panel_ylims()

        # Residual y-limits -- same helper, different per-panel function
        res_ylims = self._compute_panel_ylims(
            uniform=self.uniform_ylim_residuals,
            ylim_fn=self._residual_ylim_fn,
        )

        # --- Global chi-squared accumulators ---------------------------
        total_chi2_raw = 0.0
        total_chi2_adj = 0.0
        total_n_points = 0
        _has_nuisance = self._has_continuum or self._has_noise_floor

        # --- Render each panel -----------------------------------------
        for idx, xlim_start in enumerate(self._panel_edges):
            is_last = idx == n - 1
            panel_end = xlim_start + self._step
            xr = (xlim_start, panel_end)

            gs_inner = GridSpecFromSubplotSpec(
                2, 1,
                subplot_spec=gs_outer[idx, 0],
                height_ratios=[1.0, self.residual_height_ratio],
                hspace=0.0,
            )
            ax_spec = self.fig.add_subplot(gs_inner[0, 0])
            ax_res = self.fig.add_subplot(gs_inner[1, 0], sharex=ax_spec)
            self.subplots[idx] = (ax_spec, ax_res)

            # --- Data mask for this panel ------------------------------
            mask = (self.wave_data >= xr[0]) & (self.wave_data <= xr[1])
            panel_wave = self.wave_data[mask]
            panel_flux = self.flux_data[mask]
            panel_model_raw = self.model_flux[mask]
            panel_model = self._model_flux_adj[mask]
            panel_err_raw = (
                self.error_data[mask] if self.error_data is not None else None
            )
            panel_err = (
                self._error_adj[mask] if self._error_adj is not None else None
            )

            # --- SPECTRUM sub-panel ------------------------------------
            self._plot_observed_spectrum(
                ax_spec, panel_wave, panel_flux, panel_err_raw, deduplicate=True
            )

            # Molecule models (from MoleculeDict)
            if mol_cache:
                for m_lam, m_flux, m_color, m_label, m_name in mol_cache:
                    m_mask = (m_lam >= xr[0]) & (m_lam <= xr[1])
                    if np.any(m_mask):
                        line, = ax_spec.plot(
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

            # Explicit model components (from generated / manual)
            if self.model_components:
                for comp in self.model_components:
                    c_wave = np.asarray(comp["wave"])
                    c_flux = np.asarray(comp["flux"])
                    c_mask = (c_wave >= xr[0]) & (c_wave <= xr[1])
                    if np.any(c_mask):
                        ax_spec.plot(
                            c_wave[c_mask],
                            c_flux[c_mask],
                            color=comp.get("color", "blue"),
                            linewidth=comp.get("linewidth", 0.7),
                            alpha=comp.get("alpha", 0.6),
                            label=comp.get("label", ""),
                            zorder=self._get_theme_value("zorder_model", 3),
                        )

            # Combined model as filled area (use adjusted flux)
            _hires_adj = self._model_flux_hires_adj
            if _hires_adj is not None:
                mw_mask = (self.model_wave >= xr[0]) & (
                    self.model_wave <= xr[1]
                )
                if np.any(mw_mask):
                    self._plot_summed_spectrum(
                        ax_spec,
                        self.model_wave[mw_mask],
                        _hires_adj[mw_mask],
                        label="Combined model",
                    )
            else:
                # Use the data-grid model flux for the fill
                if len(panel_wave) > 0:
                    self._plot_summed_spectrum(
                        ax_spec, panel_wave, panel_model, label="Combined model"
                    )

            # Continuum offset overlay
            if self._has_continuum and len(panel_wave) > 0:
                cont_panel = (
                    self.continuum_c0
                    + self.continuum_c1 * (panel_wave - self.lam_ref)
                )
                ax_spec.plot(
                    panel_wave, cont_panel,
                    color="dimgray", ls="-.", lw=1.0, alpha=0.7,
                    label="Continuum offset" if idx == 0 else "",
                    zorder=self._get_theme_value("zorder_model", 3) - 1,
                )

            # y-limits for spectrum panel (pre-computed in pass 1)
            ymin, ymax = spec_ylims[idx]

            ax_spec.set_xlim(*xr)
            ax_spec.set_ylim(ymin, ymax)
            ax_spec.tick_params(axis="x", labelbottom=False, labelsize=7)
            ax_spec.tick_params(axis="y", labelsize=7)
            ax_spec.xaxis.set_major_locator(
                MaxNLocator(nbins=8, prune="both")
            )
            ax_spec.yaxis.set_major_locator(
                MaxNLocator(nbins=6, prune="both")
            )

            # --- Line annotations on spectrum panel --------------------
            if self.line_list is not None and len(self.line_list) > 0:
                self._plot_line_annotations(
                    ax_spec, self.line_list, xr, ymin, ymax
                )
            if self.atomic_lines is not None and len(self.atomic_lines) > 0:
                self._plot_atomic_lines(ax_spec, self.atomic_lines, xr=xr)

            # --- RESIDUAL sub-panel (uses adjusted model + errors) ----
            residuals = panel_flux - panel_model
            if panel_err is not None and len(panel_err) == len(residuals):
                ax_res.errorbar(
                    panel_wave,
                    residuals,
                    yerr=panel_err,
                    fmt=".",
                    ms=2,
                    color=fg,
                    ecolor="lightgray",
                    elinewidth=0.5,
                    zorder=2,
                )
            else:
                ax_res.plot(
                    panel_wave,
                    residuals,
                    ".",
                    ms=2,
                    color=fg,
                    zorder=2,
                )
            ax_res.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.7)

            # Error-envelope bands on residual panel
            if self._has_noise_floor and panel_err is not None and len(panel_err) > 0:
                # Show the ORIGINAL +/-1 sigma band (pipeline errors only)
                ax_res.fill_between(
                    panel_wave, -panel_err_raw, panel_err_raw,
                    color="salmon", alpha=0.13, zorder=1,
                    label=r"$\pm 1\,\sigma_{\rm pipe}$" if idx == 0 else "",
                )
                # Show the EXPANDED +/-1 sigma_eff band (with noise floor)
                ax_res.fill_between(
                    panel_wave, -panel_err, panel_err,
                    color="dodgerblue", alpha=0.13, zorder=1,
                    label=r"$\pm 1\,\sigma_{\rm eff}$" if idx == 0 else "",
                )

            # Residual y-limits (pre-computed in pass 1)
            ax_res.set_ylim(*res_ylims[idx])

            ax_res.set_xlim(*xr)
            ax_res.tick_params(axis="x", labelbottom=True, labelsize=7)
            ax_res.tick_params(axis="y", labelsize=7)
            ax_res.xaxis.set_major_locator(
                MaxNLocator(nbins=8, prune="both")
            )
            ax_res.yaxis.set_major_locator(
                MaxNLocator(nbins=4, prune="both")
            )
            # Only label the x-axis on the very bottom residual panel
            if is_last:
                ax_res.set_xlabel("Wavelength (μm)", fontsize=8, color=fg)
            else:
                ax_res.set_xlabel("")

            # --- Shade excluded ranges ---------------------------------
            for exc_lo, exc_hi in self.excluded_ranges:
                if exc_hi >= xr[0] and exc_lo <= xr[1]:
                    for _ax in (ax_spec, ax_res):
                        _ax.axvspan(
                            max(exc_lo, xr[0]),
                            min(exc_hi, xr[1]),
                            color="lightcoral",
                            alpha=0.15,
                        )

            # Shade atomic-line exclusion windows
            if (
                self.exclude_lines_half_width is not None
                and self.atomic_lines is not None
                and len(self.atomic_lines) > 0
            ):
                hw = self.exclude_lines_half_width
                for _, arow in self.atomic_lines.iterrows():
                    wl = arow["wave"]
                    a_lo, a_hi = wl - hw, wl + hw
                    if a_hi >= xr[0] and a_lo <= xr[1]:
                        for _ax in (ax_spec, ax_res):
                            _ax.axvspan(
                                max(a_lo, xr[0]),
                                min(a_hi, xr[1]),
                                color="lightsalmon",
                                alpha=0.12,
                            )

            # --- Per-panel chi-squared ---------------------------------
            if (
                self.show_chi2_per_panel
                and panel_err is not None
                and len(panel_err) > 0
            ):
                fit_mask = self._get_fit_mask(panel_wave)
                p_raw, p_adj, n_fit = self._compute_chi2(
                    panel_flux, panel_model_raw, panel_model,
                    panel_err_raw, panel_err, fit_mask,
                )
                if n_fit > 0:
                    panel_dof = max(n_fit - 1, 1)
                    total_chi2_raw += p_raw
                    total_chi2_adj += p_adj
                    total_n_points += n_fit
                    ann = self._format_chi2_annotation(
                        p_raw, p_adj, panel_dof, _has_nuisance,
                    )
                    ax_spec.text(
                        1.01, 0.5, ann,
                        transform=ax_spec.transAxes,
                        fontsize=7, va="center", ha="left",
                        color=fg,
                    )

        # --- Global labels & legend ------------------------------------
        self.fig.supylabel("Flux Density (Jy)", fontsize=10, color=fg, x=0.01)

        # Molecule colour legend on the first spectrum panel
        if 0 in self.subplots and (mol_labels or self.model_components):
            first_ax = self.subplots[0][0]
            if mol_labels:
                BasePlot.build_molecule_legend(
                    first_ax, mol_labels, mol_colors,
                    fontsize=7,
                    ncols=4,
                    bbox_to_anchor=(0.5, 1.5),
                    use_figure_transform=False,
                )
            else:
                # Build legend from model_components
                comp_labels = [c.get("label", "") for c in self.model_components if c.get("label")]
                comp_colors = [c.get("color", "blue") for c in self.model_components if c.get("label")]
                if comp_labels:
                    BasePlot.build_molecule_legend(
                        first_ax, comp_labels, comp_colors,
                        fontsize=7,
                        ncols=4,
                        bbox_to_anchor=(0.5, 1.5),
                        use_figure_transform=False,
                    )

        # --- Total chi-squared annotation at the bottom ----------------
        if self.show_total_chi2 and self.error_data is not None:
            full_fit_mask = self._get_fit_mask(self.wave_data)
            g_raw, g_adj, g_n = self._compute_chi2(
                self.flux_data, self.model_flux, self._model_flux_adj,
                self.error_data, self._error_adj, full_fit_mask,
            )
            global_dof = max(g_n - self.n_free_params, 1)

            if _has_nuisance:
                chi2_text = (
                    f"$\\chi^2 = {g_raw:.1f}$"
                    f"    $\\mathrm{{dof}} = {global_dof}$"
                    f"    $\\chi^2_\\nu = {g_raw / global_dof:.2f}$"
                    f"\n"
                    f"$\\chi^2_{{\\mathrm{{adj}}}} = {g_adj:.1f}$"
                    f"    $\\chi^2_{{\\nu,\\mathrm{{adj}}}} = {g_adj / global_dof:.2f}$"
                )
                extras = []
                if self._has_noise_floor:
                    extras.append(f"$s = {self.noise_floor:.4f}$ Jy")
                if self._has_continuum:
                    extras.append(
                        f"$c_0 = {self.continuum_c0:.4f}$ Jy"
                        f"  $c_1 = {self.continuum_c1:.5f}$ Jy/μm"
                    )
                if extras:
                    chi2_text += "\n" + "    ".join(extras)
            else:
                chi2_text = (
                    f"$\\chi^2 = {g_raw:.1f}$"
                    f"    $\\mathrm{{dof}} = {global_dof}$"
                    f"    $\\chi^2_\\nu = {g_raw / global_dof:.2f}$"
                )

            chi2_box = dict(
                boxstyle="round,pad=0.4",
                facecolor="white",
                edgecolor="gray",
                alpha=0.85,
                linewidth=0.8,
            )
            self.fig.text(
                0.52, 0.01, chi2_text,
                ha="center", va="bottom",
                fontsize=10, color=fg,
                fontweight="bold", bbox=chi2_box,
            )

        # Apply theme
        self.apply_theme_to_figure()

    # ------------------------------------------------------------------
    # Convenience setters
    # ------------------------------------------------------------------
    def set_model_flux(self, model_flux: np.ndarray) -> None:
        """Replace the model flux array."""
        self.model_flux = np.asarray(model_flux)

    def set_model_components(
        self, components: List[Dict[str, Any]]
    ) -> None:
        """Replace the model component list."""
        self.model_components = components