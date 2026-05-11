"""
PopulationDiagramPlot — Boltzmann / rotation diagram.

Plots ``ln(4πF / (hv A_u g_u))`` vs upper-state energy *E_u* using the
computed intensity data from one or more :class:`Molecule` /
:class:`Intensity` objects.

Supports **multi-molecule / multi-component** overlays with automatic
colour coding, as well as **property-based colour mapping** (e.g.
colour by upper-level energy, transition label, Einstein A coefficient,
etc.).

Can be used standalone (notebook / script) or embedded in a GUI layout.
"""

from __future__ import annotations

from typing import (
    Optional, Tuple, List, Dict, Any, Sequence, Union, TYPE_CHECKING,
)
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import Normalize, LogNorm
import matplotlib

from .BasePlot import BasePlot
import iSLAT.Constants as c
from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine as _MoleculeLine

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.Intensity import Intensity
    from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine

# Default colour cycle used when the caller does not supply colours.
_DEFAULT_COLORS: List[str] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

class PopulationDiagramPlot(BasePlot):
    """
    Boltzmann / rotation diagram for one or more molecules.

    The plot can be created from:

    * A single :class:`Molecule` (original API, fully backward-compatible).
    * A single bare :class:`Intensity` together with physical metadata.
    * A *list* of :class:`Molecule` objects, each rendered in its own
      colour.
    * A :class:`MoleculeDict`, where every visible molecule is plotted.

    After the diagram has been generated, the points can optionally be
    re-coloured based on a stored molecular property (e.g. upper-level
    energy, transition label, Einstein-A coefficient) via
    :meth:`color_by`.

    Parameters
    ----------
    molecule : Molecule, optional
        Single molecule whose intensity data is used.  Mutually
        exclusive with *intensity* and *molecules*.
    molecules : list[Molecule] | MoleculeDict, optional
        Multiple molecules to overlay.  Each is colour-coded
        automatically (or by its ``.color`` attribute).  Mutually
        exclusive with *molecule* and *intensity*.
    intensity : Intensity, optional
        Pre-built :class:`Intensity` object.  When provided, *radius*
        and *distance* must also be supplied.  Mutually exclusive with
        *molecule* and *molecules*.
    name : str, optional
        Display name used in the title when *intensity* is given
        directly.  Ignored when *molecule* / *molecules* is provided.
    color : str, optional
        Marker colour when *intensity* is given directly.
    radius : float, optional
        Emitting radius in AU — required when *intensity* is given
        directly.  Defaults to ``1.0``.
    distance : float, optional
        Distance to source in pc — required when *intensity* is given
        directly.  Defaults to ``160.0``.
    highlight_lines : list, optional
        List of ``(MoleculeLine, intensity, tau)`` tuples rendered as
        larger scatter points.
    figsize : tuple, optional
        Defaults to ``(6, 5)``.
    ax : Axes, optional
        Pre-existing axes for embedding.
    """

    def __init__(
        self,
        molecule: Optional["Molecule"] = None,
        *,
        molecules: Optional[Union[Sequence["Molecule"], "MoleculeDict"]] = None,
        intensity: Optional["Intensity"] = None,
        name: Optional[str] = None,
        color: Optional[str] = None,
        radius: Optional[float] = None,
        distance: Optional[float] = None,
        highlight_lines: Optional[
            List[Tuple["MoleculeLine", float, Optional[float]]]
        ] = None,
        figsize: Optional[Tuple[float, float]] = None,
        ax: Optional[Axes] = None,
        **kwargs,
    ):
        super().__init__(figsize=figsize or (6, 5), **kwargs)

        # ---- Validate mutual-exclusion --------------------------------
        n_sources = sum(
            x is not None for x in (molecule, molecules, intensity)
        )
        if n_sources > 1:
            raise ValueError(
                "Provide only one of 'molecule', 'molecules', or "
                "'intensity', not more than one."
            )

        self.molecule = molecule

        # Multi-molecule storage (normalised to a list later)
        self._molecules_input = molecules

        # Bare-intensity path
        self._intensity_obj: Optional["Intensity"] = intensity
        self._intensity_name: str = name or "Unknown"
        self._intensity_color: Optional[str] = color
        self._intensity_radius: float = radius if radius is not None else 1.0
        self._intensity_distance: float = (
            distance if distance is not None else 160.0
        )

        self.highlight_lines = highlight_lines
        self._external_ax = ax

        # Cached per-component data (populated by generate_plot)
        self._component_data: List[Dict[str, Any]] = []

        # Active colour-map state (set by color_by / clear_color_mapping)
        self._color_mapping: Optional[Dict[str, Any]] = None
        # Tracked colorbar so it can be removed on the next regeneration
        self._colorbar = None
        # Axis configuration
        self._x_prop: str = "eu"
        self._y_prop: str = "rd_yax"
        self._x_log: bool = False
        self._y_log: bool = False
        # When True, the plot is locked to all active/visible molecules in
        # _molecules_dict_ref and regenerates automatically when the set changes.
        self._all_molecules_mode: bool = False
        self._molecules_dict_ref = None   # weak reference target (MoleculeDict)
        self._molecules_change_cb = None  # registered callback handle
        # Cache for the last render_active_lines call so active scatter can
        # be restored after generate_plot clears the axes.
        self._active_lines_cache: Optional[Dict[str, Any]] = None

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------
    @property
    def ax(self) -> Axes:
        return self._ax

    @property
    def component_data(self) -> List[Dict[str, Any]]:
        """Per-component data computed during the last ``generate_plot`` call.

        Each element is a dict with keys: ``'name'``, ``'color'``,
        ``'eu'``, ``'rd_yax'``, ``'wavelength'``, ``'intens'``,
        ``'a_stein'``, ``'g_up'``, ``'lev_up'``, ``'lev_low'``,
        ``'e_low'``, ``'valid_mask'``.
        """
        return self._component_data

    # ------------------------------------------------------------------
    # Active-lines API implementation
    # ------------------------------------------------------------------
    def render_active_lines(
        self,
        line_data: List[Tuple["MoleculeLine", float, Any]],
        active_lines: List[Any],
        *,
        molecule: Optional["Molecule"] = None,
        threshold: float = 0.0,
        color: str = "green",
        **kwargs,
    ) -> Any:
        """Render scatter points for *line_data* on the population diagram.

        Appends / updates entries in *active_lines* and returns the
        :class:`~matplotlib.collections.PathCollection` scatter artist.

        Parameters
        ----------
        line_data : list[tuple]
            ``(MoleculeLine, intensity, tau)`` triples.
        active_lines : list
            Mutable list updated in-place with scatter artist references.
        molecule : Molecule
            Active molecule (required for ``radius`` and ``distance``).
        threshold : float
            Fraction (0-1) of max intensity below which lines are excluded.
        color : str
            Scatter point colour.

        Returns
        -------
        PathCollection or None
            The scatter artist, or *None* if nothing was plotted.
        """
        if (
            not line_data
            or self._ax is None
            or molecule is None
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

        x_vals: List[float] = []
        y_vals: List[float] = []
        # Each entry: (info_dict, al_idx) where al_idx is the index into
        # active_lines for the corresponding vline (= position among the
        # threshold-passing lines in line_data order).
        value_data_list: List[Tuple[Dict[str, Any], int]] = []

        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity

        mol_name = getattr(molecule, "name", "") if molecule else ""

        # al_idx tracks position among threshold-passing lines so it stays
        # aligned with what LineInspectionPlot.render_active_lines added to
        # active_lines (it increments for every threshold-passing line,
        # whether or not the pop-diagram NaN check also passes).
        al_idx = 0

        for line, intensity, tau_val in line_data:
            frac = intensity / max_intensity
            if frac < threshold:
                continue  # skipped by both panels — al_idx unchanged

            # This line passed threshold → a vline may exist at active_lines[al_idx].
            current_al_idx = al_idx
            al_idx += 1

            if any(
                x is None
                for x in [intensity, line.a_stein, line.g_up, line.lam]
            ):
                continue  # no scatter point for this line

            F = intensity * beam_s
            freq = c.SPEED_OF_LIGHT_MICRONS / line.lam
            rd_yax = (
                np.log(4 * np.pi * F / (line.a_stein * c.PLANCK_CONSTANT * freq * line.g_up))
                if (F > 0 and line.a_stein > 0) else np.nan
            )

            # Delegate per-line metadata (including FWHM breakdown) to Intensity
            # so all physics computation lives on the data object.
            info = _Intensity.get_line_info(line, intensity, tau_val, molecule=molecule)
            info["rd_yax"] = rd_yax  # computed here using beam_s
            info["intensity_percent"] = frac * 100
            info["molecule_name"] = mol_name  # needed for pick-highlight routing

            # Build the axis-property lookup map; FWHM values come from info.
            _LINE_MAP = {
                "eu":                    line.e_up,
                "e_low":                 line.e_low,
                "rd_yax":                rd_yax,
                "wavelength":            line.lam,
                "intens":                intensity,
                "a_stein":               line.a_stein,
                "g_up":                  line.g_up,
                "g_low":                 line.g_low,
                "tau":                   tau_val,
                "fwhm_instrumental_kms": info.get("fwhm_instrumental_kms"),
                "fwhm_convolved_kms":    info.get("fwhm_convolved_kms"),
            }
            xv = _LINE_MAP.get(self._x_prop)
            yv = _LINE_MAP.get(self._y_prop)
            try:
                xv = float(xv) if xv is not None else None
                yv = float(yv) if yv is not None else None
            except (TypeError, ValueError):
                xv = yv = None
            if xv is None or yv is None or np.isnan(xv) or np.isnan(yv):
                continue  # no scatter point for this line

            x_vals.append(xv)
            y_vals.append(yv)
            value_data_list.append((info, current_al_idx))

        if not x_vals:
            return None

        ax = self._ax
        sc = ax.scatter(
            x_vals, y_vals, s=30,
            color=color, edgecolors="black", picker=True,
        )

        for scatter_idx, (info, al_idx) in enumerate(value_data_list):
            info["_scatter_point_index"] = scatter_idx
            if al_idx < len(active_lines):
                # Update the existing vline entry with scatter artist + pop-diagram info
                active_lines[al_idx][2] = sc
                active_lines[al_idx][3].update(info)
            else:
                # No corresponding vline (pop-diagram rendered standalone)
                active_lines.append([None, None, sc, info])

        # Cache params so generate_plot can restore the scatter after a
        # full axes clear (e.g. when set_axes or color_by is called).
        self._active_lines_cache = {
            "line_data": line_data,
            "active_lines": active_lines,
            "molecule": molecule,
            "threshold": threshold,
            "color": color,
        }

        return sc

    def clear_active_lines(self, active_lines: List[Any]) -> None:
        """Remove scatter artists from *active_lines* and clear the list.

        Parameters
        ----------
        active_lines : list
            Mutable list of ``[vline, text, scatter, info]`` entries.
        """
        seen_scatter = set()
        for entry in active_lines:
            sc = entry[2] if len(entry) > 2 else None
            if sc is not None and id(sc) not in seen_scatter:
                seen_scatter.add(id(sc))
                if getattr(sc, 'axes', None) is not None:
                    try:
                        sc.remove()
                    except (ValueError, AttributeError):
                        pass
        active_lines.clear()
        # Clear the cache so generate_plot won't re-render stale scatter.
        self._active_lines_cache = None

    # ------------------------------------------------------------------
    # Plot generation
    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Generate the population diagram."""
        # Remove any existing colorbar before (re-)drawing so stale
        # colorbar axes don't accumulate on the figure.
        if self._colorbar is not None:
            try:
                self._colorbar.remove()
            except Exception:
                pass
            self._colorbar = None

        if self._external_ax is not None:
            self._ax = self._external_ax
        else:
            self._ensure_figure()
            self.fig.clf()
            self._ax = self.fig.add_subplot(111)

        ax = self._ax
        ax.clear()

        fg = self._get_theme_value("foreground", "black")

        # ---- Resolve data sources ------------------------------------
        components = self._resolve_components()

        if not components:
            ax.set_title("No molecule selected")
            return

        # ---- Compute population-diagram values for each component ----
        self._component_data = []
        all_valid_x: List[np.ndarray] = []
        all_valid_y: List[np.ndarray] = []

        for comp in components:
            cdata = self._compute_component(comp)
            if cdata is not None:
                self._component_data.append(cdata)
                mask = cdata["valid_mask"]
                x_arr = self._get_axis_array(cdata, self._x_prop)
                y_arr = self._get_axis_array(cdata, self._y_prop)
                if np.any(mask) and x_arr is not None and y_arr is not None:
                    all_valid_x.append(x_arr[mask])
                    all_valid_y.append(y_arr[mask])

        if not self._component_data:
            title = components[0]["name"] if components else "Unknown"
            ax.set_title(f"{title} - No intensity data", color=fg)
            return

        if not all_valid_x or not all_valid_y:
            names = ", ".join(cd["name"] for cd in self._component_data)
            ax.set_title(
                f"{names} - No valid data for population diagram",
                color=fg,
            )
            return

        # ---- Render scatter points -----------------------------------
        if self._color_mapping is not None:
            self._render_colormapped(ax)
        else:
            self._render_by_component(ax)

        # ---- Highlighted lines (e.g. from an inspection selection) ---
        if self.highlight_lines and self._component_data:
            first = self._component_data[0]
            self._render_highlights(ax, first["beam_s"])

        # ---- Axis limits & labels ------------------------------------
        cat_valid_x = np.concatenate(all_valid_x)
        cat_valid_y = np.concatenate(all_valid_y)
        _all_x_arrs = [self._get_axis_array(cd, self._x_prop) for cd in self._component_data]
        _all_y_arrs = [self._get_axis_array(cd, self._y_prop) for cd in self._component_data]
        all_x_data = np.concatenate([a for a in _all_x_arrs if a is not None]) if any(a is not None for a in _all_x_arrs) else cat_valid_x
        all_y_data = np.concatenate([a for a in _all_y_arrs if a is not None]) if any(a is not None for a in _all_y_arrs) else cat_valid_y

        # Set scale BEFORE computing limits so matplotlib never receives
        # non-positive limits on a log axis (which causes the squish bug).
        ax.set_xscale("log" if self._x_log else "linear")
        ax.set_yscale("log" if self._y_log else "linear")

        # For log axes restrict the data used for limit computation to
        # strictly positive values; fall back to auto-limits if none exist.
        if self._x_log:
            pos_x = cat_valid_x[cat_valid_x > 0]
            pos_all_x = all_x_data[all_x_data > 0]
            if len(pos_x) == 0 or len(pos_all_x) == 0:
                ax.autoscale(axis="x")
            else:
                _xlo = np.nanmin(pos_all_x)
                _xhi = np.nanmax(pos_x)
                _factor = (_xhi / _xlo) ** 0.05 if _xlo > 0 and _xhi > _xlo else 1.05
                ax.set_xlim(_xlo / _factor, _xhi * _factor)
        else:
            # X limits (linear)
            if self._x_prop == "eu":
                ax.set_xlim(np.nanmin(all_x_data) - 50, np.nanmax(cat_valid_x))
            else:
                _xr = np.nanmax(cat_valid_x) - np.nanmin(all_x_data)
                _xpad = _xr * 0.05 if _xr > 0 else max(abs(np.nanmin(all_x_data)) * 0.05, 1)
                ax.set_xlim(np.nanmin(all_x_data) - _xpad, np.nanmax(cat_valid_x) + _xpad)

        if self._y_log:
            pos_y = cat_valid_y[cat_valid_y > 0]
            pos_all_y = all_y_data[all_y_data > 0]
            if len(pos_y) == 0 or len(pos_all_y) == 0:
                ax.autoscale(axis="y")
            else:
                _ylo = np.nanmin(pos_y)
                _yhi = np.nanmax(pos_all_y)
                _factor = (_yhi / _ylo) ** 0.05 if _ylo > 0 and _yhi > _ylo else 1.05
                ax.set_ylim(_ylo / _factor, _yhi * _factor)
        else:
            # Y limits (linear)
            if self._y_prop == "rd_yax":
                ax.set_ylim(np.nanmin(cat_valid_y), np.nanmax(all_y_data) + 0.5)
            else:
                _yr = np.nanmax(all_y_data) - np.nanmin(cat_valid_y)
                _ypad = _yr * 0.05 if _yr > 0 else max(abs(np.nanmin(cat_valid_y)) * 0.05, 1)
                ax.set_ylim(np.nanmin(cat_valid_y) - _ypad, np.nanmax(all_y_data) + _ypad)

        ax.set_ylabel(self._get_axis_label(self._y_prop), color=fg, labelpad=-1)
        ax.set_xlabel(self._get_axis_label(self._x_prop), color=fg)

        # Title: single component → use its name; multi → generic
        if len(self._component_data) == 1:
            title = f"{self._component_data[0]['name']} Population diagram"
        else:
            title = "Population diagram"
        ax.set_title(title, fontsize="medium", color=fg)

        # Legend for multi-component mode
        if len(self._component_data) > 1 and self._color_mapping is None:
            self._render_component_legend(ax)

        # Re-render any previously plotted active-line scatter so it
        # persists across axes changes and color_by calls.
        if self._active_lines_cache is not None:
            cache = self._active_lines_cache
            al = cache["active_lines"]
            # Do NOT call al.clear() here — that would destroy the vline
            # entries added by LineInspectionPlot.render_active_lines,
            # making the vlines in the inspection panel un-pickable after
            # an axis change.  Instead, just null out the stale scatter
            # references so render_active_lines updates them in-place.
            for entry in al:
                if len(entry) > 2:
                    entry[2] = None
            self.render_active_lines(
                cache["line_data"],
                al,
                molecule=cache["molecule"],
                threshold=cache["threshold"],
                color=cache["color"],
            )

        if self.fig is not None:
            self.apply_theme_to_figure()

    # ------------------------------------------------------------------
    # Component resolution
    # ------------------------------------------------------------------
    def _resolve_components(self) -> List[Dict[str, Any]]:
        """Normalise the various input modes into a list of component dicts.

        Each dict has keys: ``'molecule'``, ``'intensity'``, ``'name'``,
        ``'color'``, ``'radius'``, ``'distance'``.
        """
        components: List[Dict[str, Any]] = []

        # --- Single Molecule ---
        if self.molecule is not None:
            components.append(self._mol_to_component(self.molecule, idx=0))
            return components

        # --- Multiple Molecules (list or MoleculeDict) ---
        if self._molecules_input is not None:
            mol_seq = self._molecules_input
            # MoleculeDict → use active set when in all-molecules mode,
            # otherwise fall back to visible molecules.
            if hasattr(mol_seq, "get_active_set") and self._all_molecules_mode:
                mol_seq = mol_seq.get_active_set()
            elif hasattr(mol_seq, "get_visible_molecules"):
                mol_seq = list(
                    mol_seq.get_visible_molecules(return_objects=True)
                )
            elif hasattr(mol_seq, "values") and not isinstance(mol_seq, list):
                mol_seq = list(mol_seq.values())
            for idx, mol in enumerate(mol_seq):
                components.append(self._mol_to_component(mol, idx=idx))
            # Ensure distinct colours when multiple components share a
            # colour — assign from the default cycle in that case.
            if len(components) > 1:
                self._deduplicate_colors(components)
            return components

        # --- Bare Intensity ---
        if self._intensity_obj is not None:
            components.append(
                {
                    "molecule": None,
                    "intensity": self._intensity_obj,
                    "name": self._intensity_name,
                    "color": self._intensity_color
                    or _DEFAULT_COLORS[0],
                    "radius": self._intensity_radius,
                    "distance": self._intensity_distance,
                }
            )
            return components

        return components

    @staticmethod
    def _deduplicate_colors(components: List[Dict[str, Any]]) -> None:
        """Assign distinct colours from the default cycle when components
        would otherwise share the same colour."""
        seen: Dict[str, int] = {}
        has_duplicates = False
        for comp in components:
            c = comp["color"]
            seen[c] = seen.get(c, 0) + 1
            if seen[c] > 1:
                has_duplicates = True
        if has_duplicates:
            for idx, comp in enumerate(components):
                comp["color"] = _DEFAULT_COLORS[
                    comp.get("color_idx", idx) % len(_DEFAULT_COLORS)
                ]

    @staticmethod
    def _mol_to_component(mol: "Molecule", idx: int = 0) -> Dict[str, Any]:
        """Convert a single Molecule into a component dict."""
        name = BasePlot.get_molecule_display_name(mol)
        color = BasePlot.get_molecule_color(mol)
        return {
            "molecule": mol,
            "intensity": None,
            "name": name,
            "color": color,
            "color_idx": idx,
            "radius": getattr(mol, "radius", 1.0),
            "distance": getattr(mol, "distance", 160.0),
        }

    # ------------------------------------------------------------------
    # Per-component computation
    # ------------------------------------------------------------------
    def _compute_component(
        self, comp: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """Compute population-diagram values for a single component.

        Delegates all physics to
        :meth:`~iSLAT.Modules.DataTypes.Intensity.Intensity.get_population_diagram_data`
        and annotates the result with the rendering metadata (name, color)
        that only the plotting layer knows.

        Returns a data dict or ``None`` if no valid data is available.
        """
        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity

        mol          = comp.get("molecule")
        intensity_obj = comp.get("intensity")
        radius       = comp["radius"]
        distance     = comp["distance"]

        if mol is not None:
            # Ensure intensity has been calculated
            if not hasattr(mol, "intensity") or mol.intensity is None:
                if hasattr(mol, "calculate_intensity"):
                    mol.calculate_intensity()
            intensity_from_mol = getattr(mol, "intensity", None)
            if intensity_from_mol is None:
                return None

            data = intensity_from_mol.get_population_diagram_data(
                radius, distance, molecule=mol, full_range=True
            )

            # Fall back to active-range when full-range returns all-NaN intensities
            if data is not None:
                intens_arr = data.get("intens")
                if (
                    intens_arr is not None
                    and not (np.isfinite(intens_arr) & (intens_arr > 0)).any()
                ):
                    data_active = intensity_from_mol.get_population_diagram_data(
                        radius, distance, molecule=mol, full_range=False
                    )
                    if data_active is not None:
                        active_intens = data_active.get("intens")
                        if (
                            active_intens is not None
                            and (np.isfinite(active_intens) & (active_intens > 0)).any()
                        ):
                            data = data_active

        elif intensity_obj is not None:
            data = intensity_obj.get_population_diagram_data(radius, distance)
        else:
            return None

        if data is None:
            return None

        return {
            "name":  comp["name"],
            "color": comp["color"],
            **data,
        }

    # ------------------------------------------------------------------
    # Rendering helpers
    # ------------------------------------------------------------------
    def _render_by_component(self, ax: Axes) -> None:
        """Render each component as a single-colour scatter series.

        When only one component is present, the default scatter colour
        comes from the theme (``scatter_main_color``, typically a
        neutral grey) to match the original single-molecule appearance.
        In multi-component mode each molecule uses its own colour so
        that the user can visually distinguish them.
        """
        single = len(self._component_data) == 1
        for cdata in self._component_data:
            color = (
                self._get_theme_value("scatter_main_color", "#838B8B")
                if single
                else cdata["color"]
            )
            x_arr = self._get_axis_array(cdata, self._x_prop)
            y_arr = self._get_axis_array(cdata, self._y_prop)
            if x_arr is None or y_arr is None:
                continue
            sc = ax.scatter(
                x_arr,
                y_arr,
                s=0.5 if single else 5,
                color=color,
                label=cdata["name"],
                alpha=0.8,
                picker=True,
            )
            # Tag with per-point wavelengths so shift+click can look up
            # the line wavelength regardless of the current x/y axis choice.
            wav = cdata.get("wavelength")
            if wav is not None:
                sc._islat_scatter_wavelengths = np.asarray(wav)

    def _render_component_legend(self, ax: Axes) -> None:
        """Add a legend showing component names and colours."""
        from matplotlib.lines import Line2D

        handles = []
        for cdata in self._component_data:
            handles.append(
                Line2D(
                    [0], [0],
                    marker="o",
                    color="w",
                    markerfacecolor=cdata["color"],
                    markersize=6,
                    label=cdata["name"],
                    linewidth=0,
                )
            )
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize="small",
            framealpha=0.7,
        )

    # ------------------------------------------------------------------
    # Colour-mapping API
    # ------------------------------------------------------------------
    def color_by(
        self,
        prop: str,
        *,
        cmap: str = "viridis",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        pmin: Optional[float] = None,
        pmax: Optional[float] = None,
        log_scale: bool = False,
        regenerate: bool = True,
    ) -> None:
        """Re-colour the scatter points by a molecular property.

        Parameters
        ----------
        prop : str
            The property name to colour by.  Supported values:

            * ``'e_up'``   — upper-level energy (K)
            * ``'e_low'``  — lower-level energy (K)
            * ``'a_stein'``— Einstein A coefficient
            * ``'g_up'``   — upper-state degeneracy
            * ``'g_low'``  — lower-state degeneracy
            * ``'wavelength'`` / ``'lam'``  — line wavelength
            * ``'intens'`` — model intensity
            * ``'lev_up'`` — upper quantum-state label (categorical)
            * ``'lev_low'``— lower quantum-state label (categorical)
            * ``'tau'``     — line-center opacity
            * ``'component'`` — which molecule / component the point
              belongs to (categorical)
            * ``'molecule'`` — alias for ``'component'``; colours each
              point by which molecule it belongs to (categorical)

        cmap : str
            Matplotlib colourmap name (default ``'viridis'``).
        vmin, vmax : float, optional
            Explicit colour-scale limits for continuous properties.
            Ignored when *pmin* / *pmax* are set.  Ignored for
            categorical properties.
        pmin, pmax : float, optional
            Percentile cutoffs (0-100) for the colour-scale limits.
            When set, the colour scale minimum / maximum is computed
            as the *pmin*-th / *pmax*-th percentile of the plotted
            values, overriding *vmin* / *vmax*.  For example,
            ``pmax=75`` caps the top colour at the 75th-percentile
            value so that the upper 25 % of the distribution all
            receive the maximum colour.
        log_scale : bool
            Use a logarithmic colour-norm (default ``False``).
        regenerate : bool
            If ``True`` (default) the plot is regenerated immediately.
        """
        self._color_mapping = {
            "prop": prop,
            "cmap": cmap,
            "vmin": vmin,
            "vmax": vmax,
            "pmin": pmin,
            "pmax": pmax,
            "log_scale": log_scale,
        }
        if regenerate:
            self.generate_plot()

    def clear_color_mapping(self, *, regenerate: bool = True) -> None:
        """Remove any property-based colour mapping and revert to
        per-component colouring.

        Parameters
        ----------
        regenerate : bool
            If ``True`` (default) the plot is regenerated immediately.
        """
        self._color_mapping = None
        if regenerate:
            self.generate_plot()

    # ------------------------------------------------------------------
    def _render_colormapped(self, ax: Axes) -> None:
        """Render all components with a single property-based colourmap."""
        mapping = self._color_mapping
        if mapping is None:
            return

        prop = mapping["prop"]
        cmap_name = mapping["cmap"]
        vmin = mapping.get("vmin")
        vmax = mapping.get("vmax")
        pmin = mapping.get("pmin")
        pmax = mapping.get("pmax")
        log_scale = mapping.get("log_scale", False)

        # Resolve aliases between user-facing names and internal keys
        _ALIASES = {
            "lam": "wavelength",
            "e_up": "eu",
        }
        internal_key = _ALIASES.get(prop, prop)
        display_prop = prop  # keep the user-facing name for labels

        # ---- Categorical properties ----------------------------------
        _CATEGORICAL = {"lev_up", "lev_low", "component", "molecule"}
        if prop in _CATEGORICAL:
            # 'molecule' is a convenience alias for 'component'
            render_prop = "component" if prop == "molecule" else prop
            self._render_categorical_colormap(ax, render_prop, cmap_name)
            return

        # ---- Continuous properties -----------------------------------
        all_x: List[np.ndarray] = []
        all_y: List[np.ndarray] = []
        all_vals: List[np.ndarray] = []
        all_wav: List[np.ndarray] = []

        for cdata in self._component_data:
            vals = cdata.get(internal_key)
            x_arr = self._get_axis_array(cdata, self._x_prop)
            y_arr = self._get_axis_array(cdata, self._y_prop)
            if vals is None or x_arr is None or y_arr is None:
                continue
            vals_arr = np.asarray(vals, dtype=float)

            # Apply the valid_mask so that only lines above the 1%-of-max
            # flux threshold are plotted and used for scale computation.
            # Without this, near-zero dim lines dominate the vmin/vmax for
            # log-distributed properties (intens, tau, a_stein) making all
            # the visible points appear as the same colour.
            mask = cdata.get("valid_mask")
            if mask is not None:
                mask = np.asarray(mask, dtype=bool)
                x_arr = x_arr[mask]
                y_arr = y_arr[mask]
                vals_arr = vals_arr[mask]
                wav_arr = cdata.get("wavelength")
                if wav_arr is not None:
                    all_wav.append(np.asarray(wav_arr)[mask])
            else:
                wav_arr = cdata.get("wavelength")
                if wav_arr is not None:
                    all_wav.append(np.asarray(wav_arr))

            if len(x_arr) == 0:
                continue
            all_x.append(x_arr)
            all_y.append(y_arr)
            all_vals.append(vals_arr)

        if not all_vals:
            self._render_by_component(ax)
            return

        eu_cat = np.concatenate(all_x)
        rd_cat = np.concatenate(all_y)
        val_cat = np.concatenate(all_vals)

        # Percentile cutoffs override explicit vmin/vmax when set.
        # Use only finite values so that -inf / nan (from log of 0) don't
        # corrupt the scale.
        finite_vals = val_cat[np.isfinite(val_cat)]
        if len(finite_vals) == 0:
            self._render_by_component(ax)
            return

        if pmin is not None:
            vmin = float(np.nanpercentile(finite_vals, float(pmin)))
        elif vmin is None:
            vmin = float(np.nanmin(finite_vals))
        if pmax is not None:
            vmax = float(np.nanpercentile(finite_vals, float(pmax)))
        elif vmax is None:
            vmax = float(np.nanmax(finite_vals))

        if log_scale:
            # LogNorm requires strictly positive limits; derive from the
            # positive-only subset if vmin ended up ≤ 0.
            pos_vals = finite_vals[finite_vals > 0]
            if len(pos_vals) == 0:
                pos_vals = np.array([1e-30])
            safe_vmin = max(vmin, float(np.nanmin(pos_vals))) if vmin <= 0 else vmin
            safe_vmin = max(safe_vmin, 1e-300)
            safe_vmax = max(vmax, safe_vmin * 10)
            norm = LogNorm(vmin=safe_vmin, vmax=safe_vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)
        cmap_obj = matplotlib.colormaps.get_cmap(cmap_name)

        sc = ax.scatter(
            eu_cat,
            rd_cat,
            c=val_cat,
            s=5,
            cmap=cmap_obj,
            norm=norm,
            alpha=0.8,
            picker=True,
        )
        # Wavelengths are already masked and concatenated above in all_wav.
        _wav_parts = all_wav
        if _wav_parts:
            sc._islat_scatter_wavelengths = np.concatenate(_wav_parts)

        # Add a colourbar and track it for later removal
        label = self._property_label(display_prop)
        fig = ax.get_figure()
        if fig is not None:
            self._colorbar = fig.colorbar(sc, ax=ax, label=label, pad=0.02)

    # ------------------------------------------------------------------
    def _render_categorical_colormap(
        self, ax: Axes, prop: str, cmap_name: str
    ) -> None:
        """Render scatter with categorical colour mapping."""
        all_x: List[np.ndarray] = []
        all_y: List[np.ndarray] = []
        all_labels: List[np.ndarray] = []

        for cdata in self._component_data:
            x_arr = self._get_axis_array(cdata, self._x_prop)
            y_arr = self._get_axis_array(cdata, self._y_prop)
            if x_arr is None or y_arr is None:
                continue
            mask = cdata.get("valid_mask")
            if mask is not None:
                mask = np.asarray(mask, dtype=bool)
                x_arr = x_arr[mask]
                y_arr = y_arr[mask]
            n = len(x_arr)
            if n == 0:
                continue
            if prop == "component":
                labels = np.full(n, cdata["name"], dtype=object)
            else:
                vals = cdata.get(prop)
                if vals is None:
                    labels = np.full(n, "unknown", dtype=object)
                else:
                    vals_arr = np.asarray(vals, dtype=str)
                    if mask is not None:
                        vals_arr = vals_arr[mask]
                    labels = vals_arr
            all_x.append(x_arr)
            all_y.append(y_arr)
            all_labels.append(labels)

        if not all_x:
            return
        eu_cat = np.concatenate(all_x)
        rd_cat = np.concatenate(all_y)
        label_cat = np.concatenate(all_labels)

        unique_labels = np.unique(label_cat)

        # In all-molecules mode, use each molecule's own .color attribute
        # so the population diagram matches the control-panel colours.
        if self._all_molecules_mode and prop == "component":
            label_to_color = {
                cdata["name"]: cdata["color"]
                for cdata in self._component_data
            }
            # Fall back to colormap for any label that somehow lacks a colour
            cmap_obj = matplotlib.colormaps.get_cmap(cmap_name).resampled(
                max(len(unique_labels), 1)
            )
            for i, lbl in enumerate(unique_labels):
                if lbl not in label_to_color:
                    label_to_color[lbl] = cmap_obj(i / max(len(unique_labels) - 1, 1))
        else:
            cmap_obj = matplotlib.colormaps.get_cmap(cmap_name).resampled(
                max(len(unique_labels), 1)
            )
            label_to_color = {
                lbl: cmap_obj(i / max(len(unique_labels) - 1, 1))
                for i, lbl in enumerate(unique_labels)
            }

        colors = np.array([label_to_color[lbl] for lbl in label_cat])
        ax.scatter(eu_cat, rd_cat, c=colors, s=5, alpha=0.8)

        # Legend for categories (limit to 20 entries max for readability)
        from matplotlib.lines import Line2D

        entries = list(unique_labels)
        if len(entries) > 20:
            entries = entries[:20]

        handles = [
            Line2D(
                [0], [0],
                marker="o", color="w",
                markerfacecolor=label_to_color[lbl],
                markersize=6, label=str(lbl), linewidth=0,
            )
            for lbl in entries
        ]
        if len(unique_labels) > 20:
            handles.append(
                Line2D(
                    [0], [0], marker="", color="w",
                    label=f"… +{len(unique_labels) - 20} more",
                    linewidth=0,
                )
            )
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize="x-small",
            framealpha=0.7,
            ncol=max(1, len(entries) // 10),
        )

    # ------------------------------------------------------------------
    @staticmethod
    def _property_label(prop: str) -> str:
        """LaTeX label for a property name (used for colourbar labels).

        Delegates to :meth:`~iSLAT.Modules.DataTypes.Intensity.Intensity.get_axis_label`
        so the label registry lives on the data object.
        """
        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity
        return _Intensity.get_axis_label(prop)

    # Axis-property helpers — the label dict and lookup now live on Intensity;
    # _AXIS_LABELS is kept here as a convenience reference for callers that
    # read it directly (e.g. GUI dropdown population), but it is sourced from
    # Intensity.AXIS_LABELS so there is a single source of truth.
    @classmethod
    def _get_all_axis_labels(cls) -> Dict[str, str]:
        """Return the full axis-label registry (sourced from Intensity)."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity
        return _Intensity.AXIS_LABELS

    @classmethod
    def _get_axis_label(cls, prop: str) -> str:
        """Human-readable axis label for a given property name.

        Delegates to :attr:`~iSLAT.Modules.DataTypes.Intensity.Intensity.AXIS_LABELS`
        so the label registry is stored on the data object.
        """
        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity
        return _Intensity.get_axis_label(prop)

    @staticmethod
    def _get_axis_array(
        cdata: Dict[str, Any], prop: str
    ) -> Optional[np.ndarray]:
        """Extract a property array from a component data dict.

        Returns ``None`` if the property is absent or cannot be cast to float.
        """
        arr = cdata.get(prop)
        if arr is None:
            return None
        try:
            return np.asarray(arr, dtype=float)
        except (ValueError, TypeError):
            return None

    def set_axes(
        self,
        x_prop: str = "eu",
        y_prop: str = "rd_yax",
        x_log: bool = False,
        y_log: bool = False,
        *,
        regenerate: bool = True,
    ) -> None:
        """Configure the plot axes.

        Parameters
        ----------
        x_prop : str
            Property key to use for the X axis (e.g. ``'eu'``, ``'wavelength'``).
        y_prop : str
            Property key to use for the Y axis (e.g. ``'rd_yax'``, ``'intens'``).
        x_log : bool
            Whether to use a logarithmic scale for the X axis.
        y_log : bool
            Whether to use a logarithmic scale for the Y axis.
        regenerate : bool
            If ``True`` (default) the plot is regenerated immediately.
        """
        self._x_prop = x_prop
        self._y_prop = y_prop
        self._x_log = x_log
        self._y_log = y_log
        if regenerate:
            self.generate_plot()

    # ------------------------------------------------------------------
    # Highlighted-lines overlay
    # ------------------------------------------------------------------
    def _render_highlights(self, ax: Axes, beam_s: float) -> None:
        """Overlay highlighted lines as larger scatter points."""
        if not self.highlight_lines:
            return
        color = self._get_theme_value("active_scatter_line_color", "green")

        from iSLAT.Modules.DataTypes.Intensity import Intensity as _Intensity

        x_vals: List[float] = []
        y_vals: List[float] = []

        for line_obj, intensity, tau in self.highlight_lines:
            if any(
                v is None
                for v in [intensity, line_obj.a_stein, line_obj.g_up, line_obj.lam]
            ):
                continue

            # Pre-compute derived values used by multiple axis options
            F = intensity * beam_s
            freq = c.SPEED_OF_LIGHT_MICRONS / line_obj.lam
            rd = np.log(
                4 * np.pi * F
                / (line_obj.a_stein * c.PLANCK_CONSTANT * freq * line_obj.g_up)
            ) if (F > 0 and line_obj.a_stein > 0) else np.nan

            # Resolve the active molecule for FWHM lookup.
            _hl_mol = self.molecule or (
                list(self._molecules_input.get_visible_molecules(return_objects=True))[0]
                if self._molecules_input is not None
                and hasattr(self._molecules_input, "get_visible_molecules")
                else None
            )

            # Delegate per-line metadata (including FWHM) to Intensity.get_line_info
            # so all physics computation lives on the data object.
            _hl_info = _Intensity.get_line_info(line_obj, intensity, tau, molecule=_hl_mol)

            def _resolve(prop: str) -> Optional[float]:
                """Map a property key to a scalar value for this line."""
                _LINE_MAP = {
                    "eu":                    line_obj.e_up,
                    "e_low":                 line_obj.e_low,
                    "rd_yax":                rd,
                    "wavelength":            line_obj.lam,
                    "intens":                intensity,
                    "a_stein":               line_obj.a_stein,
                    "g_up":                  line_obj.g_up,
                    "g_low":                 line_obj.g_low,
                    "tau":                   tau,
                    "fwhm_instrumental_kms": _hl_info.get("fwhm_instrumental_kms"),
                    "fwhm_convolved_kms":    _hl_info.get("fwhm_convolved_kms"),
                }
                val = _LINE_MAP.get(prop)
                try:
                    return float(val) if val is not None else None
                except (TypeError, ValueError):
                    return None

            xv = _resolve(self._x_prop)
            yv = _resolve(self._y_prop)
            if xv is None or yv is None or np.isnan(xv) or np.isnan(yv):
                continue
            x_vals.append(xv)
            y_vals.append(yv)

        if x_vals:
            ax.scatter(
                x_vals, y_vals, s=30, color=color,
                edgecolors="black", zorder=5,
            )

    # ------------------------------------------------------------------
    # Bare-intensity table builder
    # ------------------------------------------------------------------
    def _build_table_from_intensity(self) -> Optional[pd.DataFrame]:
        """Build an intensity table from the stored bare Intensity object.

        Backward-compatible wrapper around
        :meth:`_build_table_from_intensity_obj`.
        """
        return self._build_table_from_intensity_obj(self._intensity_obj)

    @staticmethod
    def _build_table_from_intensity_obj(
        intensity_obj: Optional["Intensity"],
    ) -> Optional[pd.DataFrame]:
        """Build an intensity table from an arbitrary Intensity object."""
        if intensity_obj is None:
            return None
        if hasattr(intensity_obj, "build_table"):
            df = intensity_obj.build_table(full_range=True)
        else:
            return None
        if df is not None and hasattr(df, "index") and not df.empty:
            return df
        return None

    # ------------------------------------------------------------------
    # Public mutators (backward-compatible)
    # ------------------------------------------------------------------
    def set_molecule(self, molecule: "Molecule") -> None:
        """Switch to a single molecule and regenerate."""
        self._exit_all_molecules_mode()
        self.molecule = molecule
        self._molecules_input = None
        self._intensity_obj = None
        self._color_mapping = None
        self.generate_plot()

    def set_molecules(
        self, molecules: Union[Sequence["Molecule"], "MoleculeDict"],
    ) -> None:
        """Switch to multiple molecules and regenerate.

        When *molecules* is a :class:`MoleculeDict` the plot enters
        *all-molecules mode*: it automatically colours each component by
        its own molecule colour and re-renders whenever the active/visible
        set changes.

        Parameters
        ----------
        molecules : list[Molecule] | MoleculeDict
            The molecules to overlay on the population diagram.
        """
        self.molecule = None
        self._molecules_input = molecules
        self._intensity_obj = None
        self.generate_plot()

        # If we were handed a MoleculeDict, enter persistent all-molecules mode
        if hasattr(molecules, 'add_active_molecule_change_callback'):
            self._enter_all_molecules_mode(molecules)
        else:
            # Plain list — not persistent, clear any prior mode
            self._exit_all_molecules_mode()

    # ------------------------------------------------------------------
    # All-molecules mode helpers
    # ------------------------------------------------------------------
    def _enter_all_molecules_mode(self, molecules_dict) -> None:
        """Lock the plot to a MoleculeDict and register a change callback."""
        # Unregister any previous callback first
        self._exit_all_molecules_mode()
        self._all_molecules_mode = True
        self._molecules_dict_ref = molecules_dict
        # Automatically colour by molecule (uses each mol's own colour)
        self._color_mapping = {"prop": "molecule", "cmap": "tab10"}
        # Register for active-molecule and comparison-molecule changes
        cb = self._on_molecules_changed
        self._molecules_change_cb = cb
        try:
            molecules_dict.add_active_molecule_change_callback(cb)
        except Exception:
            pass
        try:
            molecules_dict.add_comparison_molecule_change_callback(
                self._on_comparison_changed
            )
        except Exception:
            pass

    def _exit_all_molecules_mode(self) -> None:
        """Unregister change callbacks and clear all-molecules mode."""
        if self._all_molecules_mode and self._molecules_dict_ref is not None:
            try:
                self._molecules_dict_ref.remove_active_molecule_change_callback(
                    self._molecules_change_cb
                )
            except Exception:
                pass
            try:
                self._molecules_dict_ref.remove_comparison_molecule_change_callback(
                    self._on_comparison_changed
                )
            except Exception:
                pass
        self._all_molecules_mode = False
        self._molecules_dict_ref = None
        self._molecules_change_cb = None

    def _on_molecules_changed(self, old_molecule=None, new_molecule=None) -> None:
        """Callback fired when the active molecule changes in all-molecules mode."""
        if not self._all_molecules_mode or self._molecules_dict_ref is None:
            return
        self._molecules_input = self._molecules_dict_ref
        self.generate_plot()

    def _on_comparison_changed(self, comparison_molecules_list=None) -> None:
        """Callback fired when comparison molecules change in all-molecules mode."""
        if not self._all_molecules_mode or self._molecules_dict_ref is None:
            return
        self._molecules_input = self._molecules_dict_ref
        self.generate_plot()

    def set_intensity(
        self,
        intensity: "Intensity",
        *,
        name: Optional[str] = None,
        color: Optional[str] = None,
        radius: Optional[float] = None,
        distance: Optional[float] = None,
    ) -> None:
        """Switch to a bare :class:`Intensity` object and regenerate.

        Parameters
        ----------
        intensity : Intensity
            The Intensity object to use.
        name, color, radius, distance
            Optional metadata overrides.  Values that are ``None`` keep
            their previously stored value.
        """
        self.molecule = None
        self._molecules_input = None
        self._intensity_obj = intensity
        if name is not None:
            self._intensity_name = name
        if color is not None:
            self._intensity_color = color
        if radius is not None:
            self._intensity_radius = radius
        if distance is not None:
            self._intensity_distance = distance
        self._color_mapping = None
        self.generate_plot()

    def add_molecule(self, molecule: "Molecule") -> None:
        """Add a molecule to an existing multi-molecule diagram.

        If the diagram currently shows a single molecule, it is
        converted to multi-molecule mode automatically.
        """
        existing: List["Molecule"] = []
        if self.molecule is not None:
            existing.append(self.molecule)
            self.molecule = None
        if self._molecules_input is not None:
            seq = self._molecules_input
            if hasattr(seq, "values"):
                existing.extend(seq.values())
            else:
                existing.extend(seq)
        existing.append(molecule)
        self._molecules_input = existing
        self._intensity_obj = None
        self.generate_plot()

    def remove_molecule(self, name: str) -> None:
        """Remove a molecule by name from the multi-molecule diagram."""
        if self._molecules_input is None:
            if self.molecule is not None and getattr(self.molecule, "name", None) == name:
                self.molecule = None
                self.generate_plot()
            return
        seq = self._molecules_input
        if hasattr(seq, "values"):
            seq = list(seq.values())
        else:
            seq = list(seq)
        filtered = [m for m in seq if getattr(m, "name", None) != name]
        if len(filtered) == 0:
            self._molecules_input = None
        elif len(filtered) == 1:
            self.molecule = filtered[0]
            self._molecules_input = None
        else:
            self._molecules_input = filtered
        self.generate_plot()