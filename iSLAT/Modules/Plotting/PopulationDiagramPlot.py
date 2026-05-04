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
        # When True, the plot is locked to all active/visible molecules in
        # _molecules_dict_ref and regenerates automatically when the set changes.
        self._all_molecules_mode: bool = False
        self._molecules_dict_ref = None   # weak reference target (MoleculeDict)
        self._molecules_change_cb = None  # registered callback handle

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
        all_valid_rd: List[np.ndarray] = []
        all_valid_eu: List[np.ndarray] = []

        for comp in components:
            cdata = self._compute_component(comp)
            if cdata is not None:
                self._component_data.append(cdata)
                mask = cdata["valid_mask"]
                if np.any(mask):
                    all_valid_rd.append(cdata["rd_yax"][mask])
                    all_valid_eu.append(cdata["eu"][mask])

        if not self._component_data:
            title = components[0]["name"] if components else "Unknown"
            ax.set_title(f"{title} - No intensity data", color=fg)
            return

        if not all_valid_rd:
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
        cat_valid_rd = np.concatenate(all_valid_rd)
        cat_valid_eu = np.concatenate(all_valid_eu)
        all_rd = np.concatenate([cd["rd_yax"] for cd in self._component_data])
        all_eu = np.concatenate([cd["eu"] for cd in self._component_data])

        ax.set_ylim(np.nanmin(cat_valid_rd), np.nanmax(all_rd) + 0.5)
        ax.set_xlim(np.nanmin(all_eu) - 50, np.nanmax(cat_valid_eu))
        ax.set_ylabel(
            r"ln(4πF/(hν$A_{u}$$g_{u}$))", color=fg, labelpad=-1
        )
        ax.set_xlabel(r"$E_{u}$ (K)", color=fg)

        # Title: single component → use its name; multi → generic
        if len(self._component_data) == 1:
            title = f"{self._component_data[0]['name']} Population diagram"
        else:
            title = "Population diagram"
        ax.set_title(title, fontsize="medium", color=fg)

        # Legend for multi-component mode
        if len(self._component_data) > 1 and self._color_mapping is None:
            self._render_component_legend(ax)

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

        Returns a data dict or ``None`` if no valid data is available.
        """
        mol = comp.get("molecule")
        intensity_obj = comp.get("intensity")

        if mol is not None:
            int_pars = self.get_intensity_data(mol, full_range=True)
            # Fall back to active-range if full_range returns all-NaN
            if int_pars is not None and not int_pars.empty:
                intens_vals = int_pars["intens"]
                if intens_vals.isna().all():
                    int_pars_active = self.get_intensity_data(
                        mol, full_range=False
                    )
                    if (
                        int_pars_active is not None
                        and not int_pars_active.empty
                        and not int_pars_active["intens"].isna().all()
                    ):
                        int_pars = int_pars_active
        elif intensity_obj is not None:
            int_pars = self._build_table_from_intensity_obj(intensity_obj)
        else:
            return None

        if int_pars is None or int_pars.empty:
            return None

        wavelength = np.asarray(int_pars["lam"])
        intens_mod = np.asarray(int_pars["intens"])
        Astein_mod = np.asarray(int_pars["a_stein"])
        gu = np.asarray(int_pars["g_up"])
        eu = np.asarray(int_pars["e_up"])

        # Extra columns for colour-mapping
        lev_up = np.asarray(int_pars["lev_up"]) if "lev_up" in int_pars.columns else None
        lev_low = np.asarray(int_pars["lev_low"]) if "lev_low" in int_pars.columns else None
        e_low = np.asarray(int_pars["e_low"]) if "e_low" in int_pars.columns else None
        g_low = np.asarray(int_pars["g_low"]) if "g_low" in int_pars.columns else None
        tau = np.asarray(int_pars["tau"], dtype=float) if "tau" in int_pars.columns else None

        radius = comp["radius"]
        distance = comp["distance"]

        area = np.pi * (radius * c.ASTRONOMICAL_UNIT_M * 1e2) ** 2
        dist = distance * c.PARSEC_CM
        beam_s = area / dist**2

        F = intens_mod * beam_s
        frequency = c.SPEED_OF_LIGHT_MICRONS / wavelength

        with np.errstate(divide="ignore", invalid="ignore"):
            rd_yax = np.log(
                4 * np.pi * F
                / (Astein_mod * c.PLANCK_CONSTANT * frequency * gu)
            )

        threshold = np.nanmax(F) / 100 if np.any(F > 0) else 0
        valid_mask = F > threshold

        return {
            "name": comp["name"],
            "color": comp["color"],
            "eu": eu,
            "rd_yax": rd_yax,
            "wavelength": wavelength,
            "intens": intens_mod,
            "a_stein": Astein_mod,
            "g_up": gu,
            "g_low": g_low,
            "lev_up": lev_up,
            "lev_low": lev_low,
            "e_low": e_low,
            "tau": tau,
            "valid_mask": valid_mask,
            "beam_s": beam_s,
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
            ax.scatter(
                cdata["eu"],
                cdata["rd_yax"],
                s=0.5 if single else 5,
                color=color,
                label=cdata["name"],
                alpha=0.8,
            )

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
        all_eu: List[np.ndarray] = []
        all_rd: List[np.ndarray] = []
        all_vals: List[np.ndarray] = []

        for cdata in self._component_data:
            vals = cdata.get(internal_key)
            if vals is None:
                continue
            all_eu.append(cdata["eu"])
            all_rd.append(cdata["rd_yax"])
            all_vals.append(np.asarray(vals, dtype=float))

        if not all_vals:
            self._render_by_component(ax)
            return

        eu_cat = np.concatenate(all_eu)
        rd_cat = np.concatenate(all_rd)
        val_cat = np.concatenate(all_vals)

        # Percentile cutoffs override explicit vmin/vmax when set
        if pmin is not None:
            vmin = float(np.nanpercentile(val_cat, float(pmin)))
        elif vmin is None:
            vmin = float(np.nanmin(val_cat))
        if pmax is not None:
            vmax = float(np.nanpercentile(val_cat, float(pmax)))
        elif vmax is None:
            vmax = float(np.nanmax(val_cat))

        if log_scale:
            # LogNorm requires strictly positive limits
            safe_vmin = max(vmin, 1e-30) if vmin is not None else max(float(np.nanmin(val_cat[val_cat > 0])), 1e-30)
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
        )

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
        all_eu: List[np.ndarray] = []
        all_rd: List[np.ndarray] = []
        all_labels: List[np.ndarray] = []

        for cdata in self._component_data:
            if prop == "component":
                labels = np.full(len(cdata["eu"]), cdata["name"], dtype=object)
            else:
                vals = cdata.get(prop)
                if vals is None:
                    labels = np.full(len(cdata["eu"]), "unknown", dtype=object)
                else:
                    labels = np.asarray(vals, dtype=str)
            all_eu.append(cdata["eu"])
            all_rd.append(cdata["rd_yax"])
            all_labels.append(labels)

        eu_cat = np.concatenate(all_eu)
        rd_cat = np.concatenate(all_rd)
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
        """Human-readable label for a property name."""
        _LABELS = {
            "e_up": r"$E_{u}$ (K)",
            "e_low": r"$E_{low}$ (K)",
            "a_stein": r"$A_{u}$ (s$^{-1}$)",
            "g_up": r"$g_{u}$",
            "g_low": r"$g_{low}$",
            "wavelength": r"$\lambda$ (μm)",
            "intens": "Intensity",
            "tau": r"Opacity ($\tau$)",
        }
        return _LABELS.get(prop, prop)

    # ------------------------------------------------------------------
    # Highlighted-lines overlay
    # ------------------------------------------------------------------
    def _render_highlights(self, ax: Axes, beam_s: float) -> None:
        """Overlay highlighted lines as larger scatter points."""
        if not self.highlight_lines:
            return
        color = self._get_theme_value("active_scatter_line_color", "green")
        e_ups: List[float] = []
        rd_vals: List[float] = []
        for line_obj, intensity, _ in self.highlight_lines:
            if any(
                v is None
                for v in [intensity, line_obj.a_stein, line_obj.g_up, line_obj.lam]
            ):
                continue
            F = intensity * beam_s
            freq = c.SPEED_OF_LIGHT_MICRONS / line_obj.lam
            rd = np.log(
                4 * np.pi * F
                / (line_obj.a_stein * c.PLANCK_CONSTANT * freq * line_obj.g_up)
            )
            e_ups.append(line_obj.e_up)
            rd_vals.append(rd)
        if e_ups:
            ax.scatter(
                e_ups, rd_vals, s=30, color=color,
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