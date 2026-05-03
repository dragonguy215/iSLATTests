"""
MainPlotGrid — Three-panel composite plot for iSLAT spectral analysis.

Layout (GridSpec 2x2):
    Row 0, spanning both columns : full spectrum with molecule overlays
    Row 1, left column           : line inspection (zoomed region)
    Row 1, right column          : population (rotation) diagram

The grid can be used entirely standalone in a notebook / script, or
its axes can be shared with the GUI embedding layer.
"""

from typing import Optional, Tuple, List, Dict, Any, Union, TYPE_CHECKING
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from .BasePlot import BasePlot
from .LineInspectionPlot import LineInspectionPlot
from .PopulationDiagramPlot import PopulationDiagramPlot

try:
    import iSLAT.Constants as c
except ImportError:  # pragma: no cover
    c = None  # constants only needed for population diagram scatter

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine

class MainPlotGrid(BasePlot):
    """
    Three-panel composite plot mirroring the default iSLAT GUI layout.

    Parameters
    ----------
    wave_data : np.ndarray
        Full observed wavelength array.
    flux_data : np.ndarray
        Full observed flux array.
    molecules : MoleculeDict, optional
        Molecules to overlay.
    active_molecule : Molecule, optional
        The *active* molecule used for the population diagram and highlighted
        in the inspection panel.
    error_data : np.ndarray, optional
        Flux uncertainties.
    spectrum_range : tuple[float, float], optional
        ``(xmin, xmax)`` wavelength limits for the top spectrum panel.
        When *None* the full data range is used.  Can be updated later
        via :meth:`set_spectrum_range`.
    inspection_range : tuple[float, float], optional
        ``(xmin, xmax)`` for the line inspection sub-plot. If *None* the
        inspection and population panels are left empty until
        :meth:`set_inspection_range` is called.
    inspection_molecules : list[str] or bool, optional
        Controls which molecules appear in the line-inspection panel:

        - *None* or ``False`` (default) -- only the **active molecule** is
          shown.
        - ``True`` -- all visible molecules from the *MoleculeDict* are shown
          (same behaviour as standalone ``LineInspectionPlot``).
        - A list of molecule name strings (e.g. ``["H2O", "CO"]``) -- only
          those molecules are shown.
    line_data : list, optional
        List of ``(MoleculeLine, intensity, tau)`` tuples for line markers
        in the inspection panel.
    line_list : pd.DataFrame, optional
        Saved-line list for annotations on the main spectrum.
    atomic_lines : pd.DataFrame, optional
        Atomic-line annotations.
    figsize : tuple, optional
        Figure size.  Defaults to ``(15, 8.5)``.
    fig : Figure, optional
        Existing figure for GUI embedding.
    ax_spectrum : Axes, optional
        Pre-existing axes for the top spectrum panel.  When provided
        together with *ax_inspection* and *ax_popdiagram* the grid
        operates in **borrowed-axes mode**: :meth:`generate_plot` will
        *not* call ``fig.clf()`` or create new axes — it renders
        directly onto the supplied axes.  This is the mode used by
        ``ThreePanelView`` for GUI embedding.
    ax_inspection : Axes, optional
        Pre-existing axes for the line-inspection panel.
    ax_popdiagram : Axes, optional
        Pre-existing axes for the population-diagram panel.
    """

    def __init__(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        molecules: Optional["MoleculeDict"] = None,
        active_molecule: Optional["Molecule"] = None,
        error_data: Optional[np.ndarray] = None,
        spectrum_range: Optional[Tuple[float, float]] = None,
        inspection_range: Optional[Tuple[float, float]] = None,
        inspection_molecules: Optional[Union[bool, List[str]]] = None,
        line_data: Optional[List[Tuple["MoleculeLine", float, Any]]] = None,
        line_list: Optional[pd.DataFrame] = None,
        atomic_lines: Optional[pd.DataFrame] = None,
        figsize: Optional[Tuple[float, float]] = None,
        wave_data_obs: Optional[np.ndarray] = None,
        ax_spectrum: Optional[Axes] = None,
        ax_inspection: Optional[Axes] = None,
        ax_popdiagram: Optional[Axes] = None,
        **kwargs,
    ):
        # When all three axes are supplied we are in "borrowed-axes" mode:
        # the caller owns the figure and axes — we just render onto them.
        # Pass the figure extracted from the first supplied axes so
        # BasePlot.fig is set correctly even in borrowed mode.
        borrowed = (
            ax_spectrum is not None
            and ax_inspection is not None
            and ax_popdiagram is not None
        )
        if borrowed:
            kwargs.setdefault("fig", ax_spectrum.figure)

        super().__init__(figsize=figsize or (15, 8.5), **kwargs)
        self.wave_data = np.asarray(wave_data)
        self.flux_data = np.asarray(flux_data)
        self.error_data = np.asarray(error_data) if error_data is not None else None
        self.molecules = molecules
        self.active_molecule = active_molecule
        self.spectrum_range = spectrum_range
        self.inspection_range = inspection_range
        self.inspection_molecules = inspection_molecules
        self.line_data = line_data
        self.line_list = line_list
        self.atomic_lines = atomic_lines
        # Observer-frame wavelengths for MoleculeDict calls that apply
        # stellar RV internally.  Falls back to wave_data when not given.
        self.wave_data_obs = (np.asarray(wave_data_obs)
                              if wave_data_obs is not None
                              else self.wave_data)

        # Borrowed-axes mode flag
        self._borrowed_axes: bool = borrowed

        # Panel axes — either injected (borrowed) or created in generate_plot
        self.ax_spectrum: Optional[Axes] = ax_spectrum
        self.ax_inspection: Optional[Axes] = ax_inspection
        self.ax_popdiagram: Optional[Axes] = ax_popdiagram

    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Build the three-panel layout.

        In **standalone mode** (no pre-existing axes were injected) the
        figure is cleared, three new axes are created via
        :meth:`create_three_panel_axes`, and the full theme is applied.

        In **borrowed-axes mode** (axes were supplied to ``__init__``)
        the existing axes are simply cleared and re-rendered — the
        figure layout and axes identities are *not* touched, so cached
        references held by ``InteractionHandler`` or ``PlotRenderer``
        remain valid.
        """
        if self._borrowed_axes:
            # Borrowed mode — render onto the pre-existing axes.
            # Do NOT call fig.clf() or create_three_panel_axes().
            self._render_spectrum_panel()
            self._render_inspection_panel()
            self._render_population_panel()
            # Apply theme only to the axes (figure chrome is owned by caller)
            self.apply_theme_to_figure()
            return

        # Standalone mode — full figure rebuild.
        self._ensure_figure()
        self.fig.clf()

        (self.ax_spectrum,
         self.ax_inspection,
         self.ax_popdiagram) = self.create_three_panel_axes(self.fig)

        self._render_spectrum_panel()
        self._render_inspection_panel()
        self._render_population_panel()
        self.apply_theme_to_figure()

    # ------------------------------------------------------------------
    # Panel renderers
    # ------------------------------------------------------------------
    def _render_spectrum_panel(self) -> None:
        ax = self.ax_spectrum

        # In borrowed-axes mode, preserve user zoom / pan limits across
        # re-renders so the display doesn't jump.  The sentinel (0, 1)
        # means "axes have never been rendered into yet" — skip restore.
        _restore_lims = False
        if self._borrowed_axes:
            _prev_xlim = ax.get_xlim()
            _prev_ylim = ax.get_ylim()
            if _prev_xlim != (0.0, 1.0) and _prev_ylim != (0.0, 1.0):
                _restore_lims = True

        ax.clear()

        if self.wave_data is None or len(self.wave_data) == 0:
            ax.set_title("No spectrum data loaded")
            return

        # Observed spectrum
        self._plot_observed_spectrum(
            ax, self.wave_data, self.flux_data, self.error_data
        )

        # Molecule models + summed
        if self.molecules is not None:
            self._plot_visible_molecules(ax, self.molecules, wave_data=self.wave_data,
                                         wave_data_obs=self.wave_data_obs)
            try:
                s_wave, s_flux = self.molecules.get_summed_flux(
                    self.wave_data_obs, visible_only=True
                )
                self._plot_summed_spectrum(ax, s_wave, s_flux)
            except Exception:
                pass

        # Apply spectrum_range if set, otherwise use full data range
        if self.spectrum_range is not None:
            ax.set_xlim(*self.spectrum_range)
            xr = self.spectrum_range
        else:
            xr = (float(np.nanmin(self.wave_data)), float(np.nanmax(self.wave_data)))

        # Line annotations
        ymin = float(ax.get_ylim()[0]) if ax.get_ylim()[0] != 0 else -0.005
        ymax = float(ax.get_ylim()[1])

        if self.line_list is not None:
            self._plot_line_annotations(ax, self.line_list, xr, ymin, ymax)
        if self.atomic_lines is not None:
            self._plot_atomic_lines(ax, self.atomic_lines, xr=xr)

        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Flux density (Jy)")
        ax.set_title("Full Spectrum with Line Inspection")
        self._update_legend(ax)

        # Restore user zoom / pan limits in borrowed-axes mode
        if _restore_lims:
            ax.set_xlim(_prev_xlim)
            ax.set_ylim(_prev_ylim)

    # ------------------------------------------------------------------
    def _render_inspection_panel(self) -> None:
        ax = self.ax_inspection
        ax.clear()

        if self.inspection_range is None:
            ax.set_title("Line Inspection -- select a range")
            return

        xmin, xmax = self.inspection_range

        # Resolve which molecule(s) to show in the inspection panel.
        # Default: only the active molecule.
        lip_molecule = None
        lip_molecules = None

        if self.inspection_molecules is True:
            # Show all visible molecules
            lip_molecules = self.molecules
        elif isinstance(self.inspection_molecules, (list, tuple)):
            # Show a specific subset by name — build a lightweight copy
            if self.molecules is not None:
                from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict as _MD
                lip_molecules = _MD.__new__(_MD)
                dict.__init__(lip_molecules)
                # Initialise the internal caches that MoleculeDict.__init__
                # normally creates — without these, get_visible_molecules()
                # and other helpers will raise AttributeError.
                lip_molecules._visible_molecules = set()
                lip_molecules._summed_flux_cache = {}
                lip_molecules._global_parms = getattr(
                    self.molecules, "_global_parms", {}
                )
                for name in self.inspection_molecules:
                    if name in self.molecules:
                        lip_molecules[name] = self.molecules[name]
            if not lip_molecules:
                lip_molecule = self.active_molecule
                lip_molecules = None
        else:
            # Default (None / False): only the active molecule
            lip_molecule = self.active_molecule

        # Use a temporary LineInspectionPlot (renders onto our axes)
        lip = LineInspectionPlot(
            wave_data=self.wave_data,
            flux_data=self.flux_data,
            xmin=xmin,
            xmax=xmax,
            error_data=self.error_data,
            molecule=lip_molecule,
            molecules=lip_molecules,
            line_data=self.line_data,
            wave_data_obs=self.wave_data_obs,
            ax=ax,
            theme=self.theme,
        )
        lip.generate_plot()

    # ------------------------------------------------------------------
    def _render_population_panel(self) -> None:
        ax = self.ax_popdiagram
        ax.clear()

        if self.active_molecule is None and self.molecules is None:
            ax.set_title("Population Diagram -- no molecule selected")
            return

        # Default: always show only the active molecule.
        # Fall back to all visible molecules only when there is no active molecule to display.
        if self.active_molecule is not None:
            pdp = PopulationDiagramPlot(
                molecule=self.active_molecule,
                highlight_lines=self.line_data,
                ax=ax,
                theme=self.theme,
            )
        elif self.molecules is not None and len(self.molecules) > 0:
            pdp = PopulationDiagramPlot(
                molecules=self.molecules,
                highlight_lines=self.line_data,
                ax=ax,
                theme=self.theme,
            )
        else:
            ax.set_title("Population Diagram -- no molecule selected")
            return
        pdp.generate_plot()

    # ------------------------------------------------------------------
    # Public update helpers
    # ------------------------------------------------------------------
    def set_inspection_range(self, xmin: float, xmax: float) -> None:
        """Update the inspection range and refresh the bottom panels."""
        self.inspection_range = (xmin, xmax)
        if self.ax_inspection is not None:
            self._render_inspection_panel()
        if self.ax_popdiagram is not None:
            self._render_population_panel()

    def set_spectrum_range(
        self,
        xmin: Optional[float] = None,
        xmax: Optional[float] = None,
    ) -> None:
        """Set the wavelength range for the top spectrum panel.

        Pass *None* for both to reset to the full data range.
        """
        if xmin is None and xmax is None:
            self.spectrum_range = None
        else:
            lo = xmin if xmin is not None else float(np.nanmin(self.wave_data))
            hi = xmax if xmax is not None else float(np.nanmax(self.wave_data))
            self.spectrum_range = (lo, hi)
        if self.ax_spectrum is not None:
            self._render_spectrum_panel()

    def set_inspection_molecules(
        self, molecules: Optional[Union[bool, List[str]]] = None,
    ) -> None:
        """Control which molecules appear in the inspection panel.

        Parameters
        ----------
        molecules : None | False | True | list[str]
            *None* / *False* -- only the active molecule (default).
            *True* -- all visible molecules.
            A list of name strings -- only those molecules.
        """
        self.inspection_molecules = molecules
        if self.ax_inspection is not None:
            self._render_inspection_panel()

    def refresh(self) -> None:
        """Full refresh of all three panels."""
        if self.ax_spectrum is not None:
            self._render_spectrum_panel()
        if self.ax_inspection is not None:
            self._render_inspection_panel()
        if self.ax_popdiagram is not None:
            self._render_population_panel()

    # ==================================================================
    # Data update helpers (for GUI embedding)
    # ==================================================================
    def update_data(
        self,
        wave_data: Optional[np.ndarray] = None,
        flux_data: Optional[np.ndarray] = None,
        molecules: Optional["MoleculeDict"] = None,
        active_molecule: Optional["Molecule"] = None,
        error_data: Optional[np.ndarray] = None,
        wave_data_obs: Optional[np.ndarray] = None,
    ) -> None:
        """Replace data arrays in-place and refresh all panels.

        Only non-*None* arguments are updated; the rest keep their
        current values.  After updating, all three panels are
        re-rendered via :meth:`refresh`.

        Parameters
        ----------
        wave_data, flux_data, molecules, active_molecule, error_data,
        wave_data_obs :
            Same meaning as the constructor arguments.
        """
        if wave_data is not None:
            self.wave_data = np.asarray(wave_data)
        if flux_data is not None:
            self.flux_data = np.asarray(flux_data)
        if molecules is not None:
            self.molecules = molecules
        if active_molecule is not None:
            self.active_molecule = active_molecule
        if error_data is not None:
            self.error_data = np.asarray(error_data)
        if wave_data_obs is not None:
            self.wave_data_obs = np.asarray(wave_data_obs)
        elif wave_data is not None:
            # Keep observer-frame wavelengths in sync when only
            # wave_data was updated (caller did not pass wave_data_obs).
            self.wave_data_obs = self.wave_data

        self.refresh()

    # ------------------------------------------------------------------
    # Fast molecule visibility helpers (no full rebuild)
    # ------------------------------------------------------------------
    def set_molecule_visibility(
        self,
        molecule_name: str,
        visible: bool,
        ax: Optional[Axes] = None,
    ) -> bool:
        """Toggle visibility of an already-plotted molecule on *ax*.

        Searches for ``Line2D`` artists tagged with ``_molecule_name``
        and sets their visibility.  Returns *True* if any artists were
        found (i.e. the fast-path succeeded).

        Parameters
        ----------
        molecule_name : str
            Internal molecule name (``molecule.name``).
        visible : bool
            Target visibility state.
        ax : Axes, optional
            Axes to search.  Defaults to ``ax_spectrum``.
        """
        target = ax if ax is not None else self.ax_spectrum
        if target is None:
            return False
        found = False
        for line in target.lines:
            if getattr(line, "_molecule_name", None) == molecule_name:
                line.set_visible(visible)
                found = True
        return found

    def remove_molecule_lines(
        self,
        molecule_name: str,
        ax: Optional[Axes] = None,
    ) -> None:
        """Remove all ``Line2D`` artists for *molecule_name* from *ax*.

        Parameters
        ----------
        molecule_name : str
            Internal molecule name.
        ax : Axes, optional
            Defaults to ``ax_spectrum``.
        """
        target = ax if ax is not None else self.ax_spectrum
        if target is None:
            return
        for line in target.lines[:]:
            if getattr(line, "_molecule_name", None) == molecule_name:
                line.remove()

    def handle_molecule_visibility_change(
        self,
        molecule_name: str,
        is_visible: bool,
        *,
        force_rerender: bool = False,
    ) -> None:
        """Handle a single molecule's visibility toggle on the spectrum panel.

        1. Fast-toggle existing artists (or create from scratch).
        2. Recompute the summed spectrum fill.
        3. Rebuild the legend.

        This is the standalone equivalent of
        :meth:`PlotRenderer.handle_molecule_visibility_change`.

        Parameters
        ----------
        molecule_name : str
            The molecule whose visibility changed.
        is_visible : bool
            New visibility state.
        force_rerender : bool
            When *True*, stale artists are destroyed and recreated from
            current parameters (e.g. after a parameter change while the
            molecule was hidden).
        """
        if self.ax_spectrum is None or self.molecules is None:
            return
        if molecule_name not in self.molecules:
            return

        molecule = self.molecules[molecule_name]
        ax = self.ax_spectrum

        # Destroy stale artists when forcing a re-render
        if force_rerender and is_visible:
            self.remove_molecule_lines(molecule_name, ax)

        # Try fast visibility toggle
        toggled = self.set_molecule_visibility(molecule_name, is_visible, ax)

        # If turning ON but no artists exist, plot from scratch
        if is_visible and not toggled:
            self._plot_molecule_spectrum(
                ax, molecule,
                wave_data=self.wave_data,
                interpolate_to_input=False,
            )

        # Recompute summed spectrum
        self._update_summed_spectrum()

        # Remove lines for molecules that are no longer visible
        if self.molecules is not None:
            visible_names = {
                getattr(m, "name", None)
                for m in self.molecules.get_visible_molecules(return_objects=True)
            }
            for line in ax.lines[:]:
                name = getattr(line, "_molecule_name", None)
                if name is not None and name not in visible_names:
                    line.remove()

        # Rebuild legend
        self._update_legend(ax)

    def _update_summed_spectrum(self) -> None:
        """Recompute and redraw the summed-spectrum fill on ``ax_spectrum``.

        Removes existing ``_islat_summed`` collections, recomputes the
        visible-molecule sum, and re-plots if any visible molecules exist.
        """
        ax = self.ax_spectrum
        if ax is None:
            return
        # Remove old summed fill
        self._clear_tagged_artists(ax, "_islat_summed", lines=False)

        if self.molecules is None:
            return
        try:
            s_wave, s_flux = self.molecules.get_summed_flux(
                self.wave_data_obs, visible_only=True,
            )
            if s_wave is not None and len(s_flux) > 0 and np.any(s_flux > 0):
                self._plot_summed_spectrum(ax, s_wave, s_flux, deduplicate=True)
        except Exception:
            pass

    # ==================================================================
    # Active-line marker rendering (for GUI interaction layer)
    # ==================================================================
    def render_active_line_markers(
        self,
        line_data: List[Tuple["MoleculeLine", float, Any]],
        active_lines: List[Any],
        max_y: float,
        *,
        threshold: float = 0.0,
        color: str = "green",
        molecule_name: str = "",
        molecule_color: str = "",
    ) -> None:
        """Plot vertical line markers in the inspection panel.

        Creates ``vlines`` + text labels on ``ax_inspection`` for each
        molecular line in *line_data* that exceeds the intensity
        *threshold* (as a fraction of the strongest line).

        Each entry is appended to *active_lines* as
        ``[vline_artist, text_artist, None, info_dict]``.

        Parameters
        ----------
        line_data : list[tuple]
            ``(MoleculeLine, intensity, tau)`` triples.
        active_lines : list
            Mutable list that receives ``[line, text, scatter, info]``
            entries.
        max_y : float
            Scaling height for the strongest line.
        threshold : float
            Fraction (0-1) of max intensity below which lines are hidden.
        color : str
            Colour for the markers and labels.
        molecule_name : str
            Name of the molecule owning these lines (stored in info dict).
        molecule_color : str
            Color of the molecule (stored in info dict for pick handling).
        """
        if not line_data or self.ax_inspection is None:
            return

        intensities = [i for _, i, _ in line_data]
        max_intensity = max(intensities) if intensities else 1.0
        if max_intensity <= 0:
            return

        ax = self.ax_inspection
        for line, intensity, tau_val in line_data:
            frac = intensity / max_intensity
            if frac < threshold:
                continue
            lineheight = frac * max_y
            if lineheight <= 0:
                continue

            vline = ax.vlines(
                line.lam, 0, lineheight,
                color=color, linestyle="dashed", linewidth=1, picker=True,
            )
            text = ax.text(
                line.lam, lineheight,
                f"{line.e_up:.0f},{line.a_stein:.3f}",
                fontsize="x-small", color=color, rotation=45,
            )
            info = LineInspectionPlot.get_line_info(line, intensity, tau_val)
            info["lineheight"] = lineheight
            info["intensity_percent"] = frac * 100
            info["molecule_name"] = molecule_name
            info["molecule_color"] = molecule_color or color

            active_lines.append([vline, text, None, info])

    def render_active_line_scatter(
        self,
        line_data: List[Tuple["MoleculeLine", float, Any]],
        active_lines: List[Any],
        molecule: "Molecule",
        *,
        threshold: float = 0.0,
        color: str = "green",
    ) -> Any:
        """Plot scatter points for active lines on the population diagram.

        Returns the :class:`PathCollection` (scatter artist) so the GUI
        can enable pick-events on it.

        Parameters
        ----------
        line_data : list[tuple]
            ``(MoleculeLine, intensity, tau)`` triples.
        active_lines : list
            Mutable list — existing entries are updated in-place, new
            entries appended.
        molecule : Molecule
            Active molecule (used for ``radius``, ``distance``).
        threshold : float
            Fraction (0-1) of max intensity below which lines are hidden.
        color : str
            Scatter point colour.

        Returns
        -------
        PathCollection or None
            The scatter artist, or *None* if nothing was plotted.
        """
        if (
            not line_data
            or self.ax_popdiagram is None
            or molecule is None
            or c is None
        ):
            return None

        intensities = [i for _, i, _ in line_data]
        max_intensity = max(intensities) if intensities else 1.0
        if max_intensity <= 0:
            return None

        radius = getattr(molecule, "radius", 1.0)
        distance = getattr(molecule, "distance", 140.0)
        area = np.pi * (radius * c.ASTRONOMICAL_UNIT_M * 1e2) ** 2
        dist = distance * c.PARSEC_CM
        beam_s = area / dist ** 2

        e_ups: List[float] = []
        rd_yaxs: List[float] = []
        value_data_list: List[Dict[str, Any]] = []

        for line, intensity, tau_val in line_data:
            frac = intensity / max_intensity
            if frac < threshold:
                continue
            if any(
                x is None
                for x in [intensity, line.a_stein, line.g_up, line.lam]
            ):
                continue

            F = intensity * beam_s
            freq = c.SPEED_OF_LIGHT_MICRONS / line.lam
            rd_yax = np.log(
                4 * np.pi * F / (line.a_stein * c.PLANCK_CONSTANT * freq * line.g_up)
            )
            e_ups.append(line.e_up)
            rd_yaxs.append(rd_yax)

            info = LineInspectionPlot.get_line_info(line, intensity, tau_val)
            info["rd_yax"] = rd_yax
            info["intensity_percent"] = frac * 100
            value_data_list.append(info)

        if not e_ups:
            return None

        ax = self.ax_popdiagram
        sc = ax.scatter(
            e_ups, rd_yaxs, s=30,
            color=color, edgecolors="black", picker=True,
        )

        for idx, info in enumerate(value_data_list):
            info["_scatter_point_index"] = idx
            if idx < len(active_lines):
                active_lines[idx][2] = sc
                active_lines[idx][3].update(info)
            else:
                active_lines.append([None, None, sc, info])

        return sc

    # ==================================================================
    # GUI line-inspection & population-diagram rendering
    # ==================================================================
    # These methods were previously on PlotRenderer but are pure
    # rendering with no controller state, so they belong on the grid.
    # ------------------------------------------------------------------

    def render_line_inspection_plot(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        xmin: float,
        xmax: float,
        active_molecule: Optional["Molecule"] = None,
        fit_result: Optional[Any] = None,
        molecules: Optional["MoleculeDict"] = None,
        wave_data_obs: Optional[np.ndarray] = None,
        legend_visible: bool = True,
        additional_molecules: Optional[List["Molecule"]] = None,
    ) -> None:
        """Render the line-inspection panel (ax_inspection) with observed
        data and the active molecule model.

        Reuses an internal :class:`LineInspectionPlot` instance and
        optionally overlays fit results.  When *additional_molecules* are
        provided their model spectra are overlaid after the primary
        molecule (each in its own color).

        Parameters
        ----------
        wave_data, flux_data : np.ndarray
            Full observed spectrum arrays.
        xmin, xmax : float
            Wavelength bounds for the zoomed inspection region.
        active_molecule : Molecule, optional
            Active molecule to show in the panel.
        fit_result : tuple, optional
            ``(gauss_fit, fitted_wave, fitted_flux)`` to overlay.
        molecules : MoleculeDict, optional
            For matched-spectral-sampling look-ups.
        wave_data_obs : np.ndarray, optional
            Observer-frame wavelengths.
        legend_visible : bool
            Whether the legend should be shown.
        additional_molecules : list[Molecule], optional
            Comparison molecules whose model spectra are also overlaid.
        """
        ax = self.ax_inspection
        if ax is None:
            return
        ax.clear()

        if xmin is None or xmax is None or (xmax - xmin) < 0.0001:
            return

        if molecules is None:
            molecules = self.molecules
        if wave_data_obs is None:
            wave_data_obs = self.wave_data_obs

        # Reuse / create a LineInspectionPlot delegate
        if not hasattr(self, '_lip') or self._lip is None:
            self._lip = LineInspectionPlot(
                wave_data=wave_data,
                flux_data=flux_data,
                xmin=xmin,
                xmax=xmax,
                molecule=active_molecule,
                molecules=molecules,
                wave_data_obs=wave_data_obs,
                ax=ax,
                fig=self.fig,
                theme=self.theme,
                render_all_visible=False,
            )
        else:
            self._lip.wave_data = wave_data
            self._lip.flux_data = flux_data
            self._lip.xmin = xmin
            self._lip.xmax = xmax
            self._lip.molecule = active_molecule
            self._lip.molecules = molecules
            self._lip.wave_data_obs = (
                np.asarray(wave_data_obs) if wave_data_obs is not None
                else wave_data
            )
            self._lip.render_all_visible = False
            self._lip.theme = self.theme
        self._lip.generate_plot()

        # Overlay additional (comparison) molecules in their own colors
        if additional_molecules:
            use_interp = False
            target_wave = None
            if molecules is not None and hasattr(molecules, 'get_matched_sampling_wavelengths'):
                use_interp, target_wave = molecules.get_matched_sampling_wavelengths(
                    wave_data_obs if wave_data_obs is not None else wave_data
                )
                if not use_interp:
                    target_wave = None
            for mol in additional_molecules:
                if mol is not None:
                    self._lip._overlay_molecule(ax, mol, use_interp, target_wave)
            # Rebuild legend now that all molecules have been plotted
            self._lip._update_legend(ax)

        # Overlay fit results if present
        if fit_result is not None:
            current_ylim = ax.get_ylim()
            max_y = current_ylim[1] / 1.1 if current_ylim[1] > 0 else 0.15
            self.render_fit_results(fit_result, xmin, xmax, max_y,
                                    legend_visible=legend_visible)

    def render_population_diagram_for_molecule(
        self,
        molecule: "Molecule",
        force_redraw: bool = False,
    ) -> None:
        """Render the base population diagram on ``ax_popdiagram``.

        Includes a simple cache so repeated calls with the same molecule
        and unchanged parameters are skipped.

        Parameters
        ----------
        molecule : Molecule
            The molecule to visualise.
        force_redraw : bool
            When *True* the cache is bypassed.
        """
        ax = self.ax_popdiagram
        if ax is None:
            return

        # Cache check
        current_hash = None
        if molecule is not None and hasattr(molecule, '_compute_intensity_hash'):
            current_hash = (molecule.name, molecule._compute_intensity_hash())

        if not hasattr(self, '_pdp'):
            self._pdp = None
            self._pdp_molecule = None
            self._pdp_cache_key = None

        if (not force_redraw
                and self._pdp_molecule is molecule
                and current_hash is not None
                and self._pdp_cache_key == current_hash):
            return

        self._pdp_molecule = molecule
        self._pdp_cache_key = current_hash

        try:
            if self._pdp is None:
                self._pdp = PopulationDiagramPlot(
                    molecule=molecule,
                    ax=ax,
                    fig=self.fig,
                    theme=self.theme,
                )
            else:
                self._pdp.molecule = molecule
                self._pdp.theme = self.theme
            self._pdp.generate_plot()
        except Exception as e:
            ax.clear()
            mol_label = (
                BasePlot.get_molecule_display_name(molecule) if molecule
                else "Unknown"
            )
            fg = self._get_theme_value("foreground", "black")
            ax.set_title(f"{mol_label} - Error in calculation", color=fg)

    def render_fit_results(
        self,
        fit_result: Any,
        xmin: float,
        xmax: float,
        max_y: float,
        *,
        user_settings: Optional[Dict[str, Any]] = None,
        legend_visible: bool = True,
    ) -> None:
        """Overlay Gaussian-fit results on the line-inspection panel.

        Parameters
        ----------
        fit_result : tuple
            ``(gauss_fit, fitted_wave, fitted_flux)``.
        xmin, xmax : float
            Wavelength range of the current selection.
        max_y : float
            Maximum y value for scaling (unused but kept for API compat).
        user_settings : dict, optional
            iSLAT user settings — used for ``clear_old_fits`` and
            ``fit_line_uncertainty``.  Falls back to safe defaults.
        legend_visible : bool
            Whether the legend should be visible after rendering.
        """
        ax = self.ax_inspection
        if ax is None:
            return

        if user_settings is None:
            user_settings = {}

        # Clear old fit results if requested
        if user_settings.get('clear_old_fits', True):
            self._clear_old_fit_results(ax, xmin, xmax)

        try:
            gauss_fit, fitted_wave, fitted_flux = fit_result
            if gauss_fit is None or fitted_wave is None or fitted_flux is None:
                return

            fit_mask = (fitted_wave >= xmin) & (fitted_wave <= xmax)
            if not np.any(fit_mask):
                return

            fit_line = ax.plot(
                fitted_wave[fit_mask], fitted_flux[fit_mask],
                color='red', linewidth=1, label='Total Fit', linestyle='--',
            )[0]
            fit_line._islat_fit_result = True

            # Multi-component handling
            if hasattr(gauss_fit, 'params') and gauss_fit.params:
                prefixes = set()
                for pname in gauss_fit.params:
                    if '_' in pname:
                        pfx = pname.split('_')[0] + '_'
                        if pfx.startswith('g') and pfx[1:-1].isdigit():
                            prefixes.add(pfx)

                if len(prefixes) > 1:
                    try:
                        components = gauss_fit.eval_components(x=fitted_wave[fit_mask])
                        for i, pfx in enumerate(sorted(prefixes)):
                            if pfx in components:
                                comp_line = ax.plot(
                                    fitted_wave[fit_mask], components[pfx],
                                    linestyle='--', linewidth=1,
                                    label=f"Component {i+1}",
                                )[0]
                                comp_line._islat_fit_result = True
                    except Exception:
                        pass
                else:
                    sigma = user_settings.get('fit_line_uncertainty', 1.0)
                    dely = gauss_fit.eval_uncertainty(sigma=sigma)
                    fill = ax.fill_between(
                        fitted_wave, fitted_flux - dely, fitted_flux + dely,
                        color='gray', alpha=0.3,
                        label=r'3-$\sigma$ uncertainty band',
                    )
                    fill._islat_fit_result = True

        except Exception:
            pass

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
            leg = ax.get_legend()
            if leg is not None:
                leg.set_visible(legend_visible)

    # ------------------------------------------------------------------
    @staticmethod
    def _clear_old_fit_results(ax: Axes, xmin: float, xmax: float) -> None:
        """Remove existing fit-result artists that overlap ``[xmin, xmax]``."""
        for line in ax.lines[:]:
            if hasattr(line, '_islat_fit_result'):
                xdata = line.get_xdata()
                if len(xdata) > 0:
                    if np.min(xdata) <= xmax and np.max(xdata) >= xmin:
                        line.remove()

        for coll in ax.collections[:]:
            if hasattr(coll, '_islat_fit_result'):
                try:
                    paths = coll.get_paths()
                    if paths:
                        bounds = paths[0].get_extents()
                        if bounds.xmin <= xmax and bounds.xmax >= xmin:
                            coll.remove()
                except Exception:
                    coll.remove()

    def clear_active_lines(self, active_lines_list: List[Any]) -> None:
        """Remove all active-line artists (vlines, text, scatter) and clear the list.

        Parameters
        ----------
        active_lines_list : list
            Mutable list of ``[line_artist, text_artist, scatter_artist, info_dict]``
            entries.  The list is cleared in-place.
        """
        for entry in active_lines_list:
            for artist in entry[:3]:  # line, text, scatter
                if artist is not None and getattr(artist, 'axes', None) is not None:
                    try:
                        artist.remove()
                    except (ValueError, AttributeError):
                        pass
        active_lines_list.clear()