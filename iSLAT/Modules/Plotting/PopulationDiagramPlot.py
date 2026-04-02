"""
PopulationDiagramPlot — Boltzmann / rotation diagram for a single molecule.

Plots ``ln(4πF / (hv A_u g_u))`` vs upper-state energy *E_u* using the
computed intensity data from a :class:`Molecule` or :class:`Intensity`
object.

Can be used standalone (notebook / script) or embedded in a GUI layout.
"""

from typing import Optional, Tuple, List, Dict, Any, Union, TYPE_CHECKING
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .BasePlot import BasePlot
import iSLAT.Constants as c

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.Intensity import Intensity
    from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine

class PopulationDiagramPlot(BasePlot):
    """
    Boltzmann / rotation diagram for a molecule.

    The plot can be created from either a full :class:`Molecule` object **or**
    a bare :class:`Intensity` object together with the physical metadata
    (``name``, ``color``, ``radius``, ``distance``) that a ``Molecule``
    would normally carry.

    Parameters
    ----------
    molecule : Molecule, optional
        Molecule whose intensity data is used.  Intensity will be
        automatically calculated if not already done.  Mutually
        exclusive with *intensity*.
    intensity : Intensity, optional
        Pre-built :class:`Intensity` object. When provided, *radius* and
        *distance* must also be supplied.
    name : str, optional
        Display name used in the title when *intensity* is given
        directly.  Ignored when *molecule* is provided (the molecule's
        own display label is used instead).
    color : str, optional
        Marker colour when *intensity* is given directly.
    radius : float, optional
        Emitting radius in AU — required when *intensity* is given
        directly.  Defaults to ``1.0``.
    distance : float, optional
        Distance to the source in pc — required when *intensity* is
        given directly.  Defaults to ``160.0``.
    highlight_lines : list, optional
        List of ``(MoleculeLine, intensity, tau)`` tuples.  These are
        rendered as larger coloured scatter points on top of the base
        diagram.
    figsize : tuple, optional
        Defaults to ``(6, 5)``.
    ax : Axes, optional
        Pre-existing axes for embedding.
    """

    def __init__(
        self,
        molecule: Optional["Molecule"] = None,
        *,
        intensity: Optional["Intensity"] = None,
        name: Optional[str] = None,
        color: Optional[str] = None,
        radius: Optional[float] = None,
        distance: Optional[float] = None,
        highlight_lines: Optional[List[Tuple["MoleculeLine", float, Optional[float]]]] = None,
        figsize: Optional[Tuple[float, float]] = None,
        ax: Optional[Axes] = None,
        **kwargs,
    ):
        super().__init__(figsize=figsize or (6, 5), **kwargs)

        if molecule is not None and intensity is not None:
            raise ValueError(
                "Provide either 'molecule' or 'intensity', not both."
            )

        self.molecule = molecule

        # Store an Intensity object for the "bare intensity" path
        self._intensity_obj: Optional["Intensity"] = intensity
        self._intensity_name: str = name or "Unknown"
        self._intensity_color: Optional[str] = color
        self._intensity_radius: float = radius if radius is not None else 1.0
        self._intensity_distance: float = distance if distance is not None else 160.0

        self.highlight_lines = highlight_lines
        self._external_ax = ax

    @property
    def ax(self) -> Axes:
        return self._ax

    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Generate the population diagram."""
        if self._external_ax is not None:
            self._ax = self._external_ax
        else:
            self._ensure_figure()
            # Clear previous axes so regeneration doesn't stack on top
            self.fig.clf()
            self._ax = self.fig.add_subplot(111)

        ax = self._ax
        ax.clear()

        fg = self._get_theme_value("foreground", "black")

        # ---- Resolve data source (Molecule *or* bare Intensity) ----
        mol = self.molecule
        use_bare_intensity = mol is None and self._intensity_obj is not None

        if mol is None and not use_bare_intensity:
            ax.set_title("No molecule selected")
            return

        if use_bare_intensity:
            # Build the table directly from the Intensity object
            int_pars = self._build_table_from_intensity()
            display_name = self._intensity_name
            radius = self._intensity_radius
            distance = self._intensity_distance
        else:
            int_pars = self.get_intensity_data(mol, full_range=True)
            display_name = self.get_molecule_display_name(mol)
            radius = getattr(mol, "radius", 1.0)
            distance = getattr(mol, "distance", 160.0)

        if int_pars is None:
            ax.set_title(
                f"{display_name} - No intensity data",
                color=fg,
            )
            return

        # Extract arrays from the intensity table
        wavelength = np.asarray(int_pars["lam"])
        intens_mod = np.asarray(int_pars["intens"])
        Astein_mod = np.asarray(int_pars["a_stein"])
        gu = np.asarray(int_pars["g_up"])
        eu = np.asarray(int_pars["e_up"])

        area = np.pi * (radius * c.ASTRONOMICAL_UNIT_M * 1e2) ** 2
        dist = distance * c.PARSEC_CM
        beam_s = area / dist ** 2

        F = intens_mod * beam_s
        frequency = c.SPEED_OF_LIGHT_MICRONS / wavelength

        # Suppress divide-by-zero warnings for log
        with np.errstate(divide="ignore", invalid="ignore"):
            rd_yax = np.log(4 * np.pi * F / (Astein_mod * c.PLANCK_CONSTANT * frequency * gu))

        threshold = np.nanmax(F) / 100

        valid_rd = rd_yax[F > threshold]
        valid_eu = eu[F > threshold]

        if len(valid_rd) == 0 or len(valid_eu) == 0:
            ax.set_title(
                f"{display_name} - No valid data for population diagram",
                color=fg,
            )
            return

        # Base scatter
        ax.scatter(
            eu,
            rd_yax,
            s=0.5,
            color=self._get_theme_value("scatter_main_color", "#838B8B"),
        )

        # Highlighted lines (e.g. from an inspection selection)
        if self.highlight_lines:
            self._render_highlights(ax, beam_s)

        ax.set_ylim(np.nanmin(valid_rd), np.nanmax(rd_yax) + 0.5)
        ax.set_xlim(np.nanmin(eu) - 50, np.nanmax(valid_eu))
        ax.set_ylabel(r"ln(4πF/(hν$A_{u}$$g_{u}$))", color=fg, labelpad=-1)
        ax.set_xlabel(r"$E_{u}$ (K)", color=fg)
        ax.set_title(
            f"{display_name} Population diagram",
            fontsize="medium",
            color=fg,
        )

        # Apply full theme (backgrounds, spines, etc.) to the figure
        if self.fig is not None:
            self.apply_theme_to_figure()

    # ------------------------------------------------------------------
    def _build_table_from_intensity(self) -> Optional[pd.DataFrame]:
        """Build an intensity table directly from a bare Intensity object.

        Returns the same DataFrame that ``BasePlot.get_intensity_data``
        would produce, but without requiring a full :class:`Molecule`
        wrapper.
        """
        intensity_obj = self._intensity_obj
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
    def _render_highlights(self, ax: Axes, beam_s: float) -> None:
        """Overlay highlighted lines as larger scatter points."""
        if not self.highlight_lines:
            return
        color = self._get_theme_value("active_scatter_line_color", "green")
        e_ups, rd_vals = [], []
        for line_obj, intensity, _ in self.highlight_lines:
            if any(v is None for v in [intensity, line_obj.a_stein, line_obj.g_up, line_obj.lam]):
                continue
            F = intensity * beam_s
            freq = c.SPEED_OF_LIGHT_MICRONS / line_obj.lam
            rd = np.log(4 * np.pi * F / (line_obj.a_stein * c.PLANCK_CONSTANT * freq * line_obj.g_up))
            e_ups.append(line_obj.e_up)
            rd_vals.append(rd)
        if e_ups:
            ax.scatter(e_ups, rd_vals, s=30, color=color, edgecolors="black", zorder=5)

    # ------------------------------------------------------------------
    def set_molecule(self, molecule: "Molecule") -> None:
        """Switch to a different molecule and regenerate."""
        self.molecule = molecule
        self._intensity_obj = None  # clear bare intensity path
        self.generate_plot()

    # ------------------------------------------------------------------
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
        self.molecule = None  # clear molecule path
        self._intensity_obj = intensity
        if name is not None:
            self._intensity_name = name
        if color is not None:
            self._intensity_color = color
        if radius is not None:
            self._intensity_radius = radius
        if distance is not None:
            self._intensity_distance = distance
        self.generate_plot()