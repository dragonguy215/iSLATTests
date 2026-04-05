"""
StackedSpectralPanel -- Abstract composer for vertically stacked spectral panels.

Manages a collection of :class:`SpectralPanel` instances laid out in a
vertical stack using matplotlib's :class:`~matplotlib.gridspec.GridSpec`.
Each row (cell) may contain an **arbitrary number** of sub-panels
arranged vertically with configurable height ratios -- for example a
single spectrum axes (as in :class:`FullSpectrumPlot`) or a spectrum +
residual pair (as in :class:`ResidualSpectrumPlot`).

Concrete subclasses must implement :meth:`_create_cell` to produce the
panels and axes for each wavelength row, and may override hooks like
:meth:`_post_render_cell` and :meth:`_cell_height_ratios`.

Usage sketch::

    class MyStackedPlot(StackedSpectralPanel):
        def _create_cell(self, idx, xmin, xmax, gs_slot, **kw):
            ax = self.fig.add_subplot(gs_slot)
            panel = MyPanel(wave, flux, xmin, xmax, ax=ax)
            return [panel]

    plot = MyStackedPlot(wave, flux, n_panels=6)
    plot.generate_plot()
    plot.show()
"""

from abc import abstractmethod
from typing import (
    Callable,
    Optional,
    Tuple,
    List,
    Dict,
    Any,
    TYPE_CHECKING,
)

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure as MplFigure
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator

from .BasePlot import BasePlot
from .SpectralPanel import SpectralPanel

if TYPE_CHECKING:
    from matplotlib.gridspec import SubplotSpec
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict


class StackedSpectralPanel(BasePlot):
    """
    Abstract base for vertically stacked :class:`SpectralPanel` figures.

    Each row in the outer GridSpec is called a **cell**.  A cell may
    contain one or more sub-panels (axes) -- the concrete subclass
    decides via :meth:`_create_cell` and :meth:`_cell_height_ratios`.

    Parameters
    ----------
    wave_data : np.ndarray
        Full observed wavelength array.
    flux_data : np.ndarray
        Full observed flux array.
    error_data : np.ndarray, optional
        Flux uncertainties.
    molecules : MoleculeDict, optional
        Molecule collection forwarded to each child panel.
    n_panels : int
        Target number of wavelength rows.  Default 5.
    step : float, optional
        Fixed wavelength width per row.  Overrides *n_panels*.
    xlim_range : tuple[float, float], optional
        ``(start, end)`` wavelength sub-range to display.  Defaults to
        the full extent of *wave_data*.
    ymax_factor : float
        Fractional headroom above the peak flux per panel (0.2 = 20 %).
    uniform_ylim : bool
        When *True* all panels share the same vertical scale.
        Default *False*.
    hspace : float
        Vertical spacing between cells in the outer GridSpec.
        Default 0.15.
    row_height : float
        Base height (inches) per cell used for auto-figsize.
        Default 2.0.
    figsize : tuple, optional
        Explicit figure size.  When *None* it is computed from
        ``(14, row_height * n_cells)``.
    **kwargs
        Forwarded to :class:`BasePlot`.
    """

    def __init__(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        error_data: Optional[np.ndarray] = None,
        molecules: Optional["MoleculeDict"] = None,
        n_panels: int = 5,
        step: Optional[float] = None,
        xlim_range: Optional[Tuple[float, float]] = None,
        ymax_factor: float = 0.2,
        uniform_ylim: bool = False,
        hspace: float = 0.15,
        row_height: float = 2.0,
        figsize: Optional[Tuple[float, float]] = None,
        **kwargs,
    ):
        super().__init__(figsize=figsize, **kwargs)
        self.wave_data = np.asarray(wave_data)
        self.flux_data = np.asarray(flux_data)
        self.error_data = (
            np.asarray(error_data) if error_data is not None else None
        )
        self.molecules = molecules
        self.n_panels = n_panels
        self.ymax_factor = ymax_factor
        self.uniform_ylim = uniform_ylim
        self.hspace = hspace
        self.row_height = row_height

        # Wavelength range
        if xlim_range is not None:
            self._xlim_start, self._xlim_end = xlim_range
        else:
            self._xlim_start = float(np.nanmin(self.wave_data))
            self._xlim_end = float(np.nanmax(self.wave_data))

        # Panel step
        if step is not None:
            self._step = step
        else:
            self._step = (self._xlim_end - self._xlim_start) / max(
                self.n_panels, 1
            )

        # Pre-compute panel edges
        self._panel_edges: np.ndarray = np.arange(
            self._xlim_start, self._xlim_end, self._step
        )

        # Auto figsize
        if self._figsize is None:
            self._figsize = (14, self.row_height * len(self._panel_edges))

        # Storage: {row_index: list[SpectralPanel]}
        self.panels: Dict[int, List[SpectralPanel]] = {}

    # ------------------------------------------------------------------
    # Abstract factory -- subclasses produce the concrete cell contents
    # ------------------------------------------------------------------
    @abstractmethod
    def _create_cell(
        self,
        idx: int,
        xmin: float,
        xmax: float,
        gs_slot: "SubplotSpec",
        **kwargs,
    ) -> List[SpectralPanel]:
        """Create and return one or more :class:`SpectralPanel` instances
        for a single cell (row) in the stacked figure.

        The subclass is responsible for creating the axes (via
        ``self.fig.add_subplot`` or a nested ``GridSpecFromSubplotSpec``)
        and attaching them to the panels.

        Parameters
        ----------
        idx : int
            Row index (0-based).
        xmin, xmax : float
            Wavelength bounds for this row.
        gs_slot : SubplotSpec
            The outer-GridSpec slot allocated to this row.  Use it
            directly with ``fig.add_subplot(gs_slot)`` for a single-axes
            cell, or split it further with ``GridSpecFromSubplotSpec``
            for multi-axes cells.
        **kwargs
            Extra keyword arguments forwarded from :meth:`generate_plot`.

        Returns
        -------
        list[SpectralPanel]
            One or more panels that will be rendered in this row.
        """
        ...

    # ------------------------------------------------------------------
    # Optional per-row hook
    # ------------------------------------------------------------------
    def _post_render_cell(
        self,
        idx: int,
        cell_panels: List[SpectralPanel],
        is_last: bool,
    ) -> None:
        """Called after all panels in a cell have been generated.

        Subclasses can override this to add per-row tick formatting,
        labels, annotations, chi-squared boxes, etc.  The default
        implementation applies common tick locators and hides x-axis
        labels on all but the last row.
        """
        fg = self._get_theme_value("foreground", "black")
        for p_idx, panel in enumerate(cell_panels):
            ax = panel.ax
            if ax is None:
                continue
            ax.tick_params(axis="x", labelsize=7)
            ax.tick_params(axis="y", labelsize=7)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=8, prune="both"))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))
            is_bottom_ax = (p_idx == len(cell_panels) - 1)
            if is_last and is_bottom_ax:
                ax.set_xlabel(
                    "Wavelength (\u03bcm)", fontsize=8, color=fg,
                )
            else:
                ax.tick_params(axis="x", labelbottom=is_bottom_ax)

    # ------------------------------------------------------------------
    # y-limit computation
    # ------------------------------------------------------------------
    @staticmethod
    def _default_ylim_fn(
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        mask: np.ndarray,
        ymax_factor: float = 0.2,
    ) -> Tuple[float, float]:
        """Compute default ``(ymin, ymax)`` from observed flux."""
        if np.any(mask):
            peak = float(np.nanmax(flux_data[mask]))
            return (-0.005, peak + peak * ymax_factor)
        return (-0.005, 0.1)

    def _compute_panel_ylims(
        self,
        uniform: Optional[bool] = None,
        ylim_fn: Optional[
            Callable[[np.ndarray], Tuple[float, float]]
        ] = None,
    ) -> List[Tuple[float, float]]:
        """Compute per-cell y-limits.

        Parameters
        ----------
        uniform : bool, optional
            When *True* all entries share a global min/max.  Defaults to
            :attr:`uniform_ylim`.
        ylim_fn : callable, optional
            ``fn(mask) -> (ymin, ymax)``.  Defaults to a closure around
            :meth:`_default_ylim_fn`.
        """
        if uniform is None:
            uniform = self.uniform_ylim

        if ylim_fn is None:
            def ylim_fn(mask: np.ndarray) -> Tuple[float, float]:
                return self._default_ylim_fn(
                    self.wave_data, self.flux_data, mask, self.ymax_factor,
                )

        ylims: List[Tuple[float, float]] = []
        for xlim_start in self._panel_edges:
            panel_end = xlim_start + self._step
            mask = (self.wave_data >= xlim_start) & (
                self.wave_data <= panel_end
            )
            ylims.append(ylim_fn(mask))

        if uniform and ylims:
            g_ymin = min(yl[0] for yl in ylims)
            g_ymax = max(yl[1] for yl in ylims)
            ylims = [(g_ymin, g_ymax)] * len(ylims)

        return ylims

    # ------------------------------------------------------------------
    # Post-render annotation helpers (delegate to each panel)
    # ------------------------------------------------------------------
    def plot_atomic_lines(
        self,
        atomic_df: Any,
        tag: str = "_islat_atomic_line",
    ) -> None:
        """Draw atomic-line markers on every panel that supports it."""
        if atomic_df is None or (hasattr(atomic_df, "__len__") and len(atomic_df) == 0):
            return
        for cell_panels in self.panels.values():
            for panel in cell_panels:
                if hasattr(panel, "plot_atomic_lines"):
                    panel.plot_atomic_lines(atomic_df, tag=tag)

    def remove_atomic_lines(
        self,
        tag: str = "_islat_atomic_line",
    ) -> None:
        """Remove atomic-line artists from every panel."""
        for cell_panels in self.panels.values():
            for panel in cell_panels:
                if hasattr(panel, "remove_atomic_lines"):
                    panel.remove_atomic_lines(tag=tag)

    def plot_saved_lines(
        self,
        line_data: Any,
        tag: str = "_islat_saved_line",
    ) -> None:
        """Draw saved-line annotations on every panel that supports it."""
        if line_data is None or (hasattr(line_data, "__len__") and len(line_data) == 0):
            return
        for cell_panels in self.panels.values():
            for panel in cell_panels:
                if hasattr(panel, "plot_saved_lines"):
                    panel.plot_saved_lines(line_data, tag=tag)

    def remove_saved_lines(
        self,
        tag: str = "_islat_saved_line",
    ) -> None:
        """Remove saved-line artists from every panel."""
        for cell_panels in self.panels.values():
            for panel in cell_panels:
                if hasattr(panel, "remove_saved_lines"):
                    panel.remove_saved_lines(tag=tag)

    # ------------------------------------------------------------------
    # Top-level plot generation
    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Build the stacked multi-panel figure.

        For each panel edge a cell is created (via :meth:`_create_cell`),
        each child panel's :meth:`~SpectralPanel.generate_plot` is
        invoked, and :meth:`_post_render_cell` is called for per-row
        decoration.
        """
        n = len(self._panel_edges)
        self._ensure_figure()
        self.fig.clf()
        self.panels.clear()

        # Disable constrained_layout for uniform cell heights.
        self.fig.set_layout_engine(None)

        gs = GridSpec(
            nrows=n,
            ncols=1,
            figure=self.fig,
            hspace=self.hspace,
        )

        self.fig.subplots_adjust(
            left=0.06, right=0.94, top=0.93, bottom=0.06,
        )

        fg = self._get_theme_value("foreground", "black")

        for idx, xlim_start in enumerate(self._panel_edges):
            is_last = idx == n - 1
            panel_end = xlim_start + self._step
            xmin, xmax = xlim_start, panel_end

            # Delegate cell creation to the concrete subclass.
            cell_panels = self._create_cell(
                idx, xmin, xmax, gs[idx, 0], **kwargs,
            )
            self.panels[idx] = cell_panels

            # Render each sub-panel in the cell.
            for panel in cell_panels:
                panel.generate_plot(**kwargs)

            # Per-row decoration hook.
            self._post_render_cell(idx, cell_panels, is_last)

        # --- Global figure labels --------------------------------------
        self.fig.supylabel(
            "Flux Density (Jy)", fontsize=10, color=fg, x=0.01,
        )

        # Apply theme
        self.apply_theme_to_figure()
