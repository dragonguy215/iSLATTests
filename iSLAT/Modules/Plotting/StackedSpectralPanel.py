"""
StackedSpectralPanel -- Abstract composer for vertically stacked spectral panels.

Manages a collection of :class:`SpectralPanel` instances laid out in a
vertical stack using matplotlib's :class:`~matplotlib.gridspec.GridSpec`.
This provides the same panel-stacking logic currently used by
:class:`FullSpectrumPlot` and :class:`ResidualSpectrumPlot` but
generalised so that **any** :class:`SpectralPanel` subclass can be
composed into a multi-panel figure.

Concrete subclasses must implement :meth:`_create_panel` to produce
the appropriate :class:`SpectralPanel` for each wavelength range, and
may override :meth:`_post_render_panel` to add per-row decorations.

Usage sketch::

    class MyStackedPlot(StackedSpectralPanel):
        def _create_panel(self, xmin, xmax, ax, **kw):
            return MyConcretePanel(wave, flux, xmin, xmax, ax=ax)

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
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict


class StackedSpectralPanel(BasePlot):
    """
    Abstract base for vertically stacked :class:`SpectralPanel` figures.

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
        Vertical spacing between rows in the outer GridSpec.
        Default 0.15.
    row_height : float
        Base height (inches) per row used for auto-figsize.  Default 2.0.
    figsize : tuple, optional
        Explicit figure size.  When *None* it is computed from
        ``(14, row_height * n_panels)``.
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

        # Child panel storage: {row_index: SpectralPanel}
        self.panels: Dict[int, SpectralPanel] = {}
        # Axes storage (mirrors child panels for external access)
        self.subplots: Dict[int, Axes] = {}

    # ------------------------------------------------------------------
    # Abstract factory -- subclasses produce the concrete panel type
    # ------------------------------------------------------------------
    @abstractmethod
    def _create_panel(
        self,
        xmin: float,
        xmax: float,
        ax: Axes,
        **kwargs,
    ) -> SpectralPanel:
        """Create and return a :class:`SpectralPanel` for one row.

        This is the factory method that concrete subclasses must
        implement.  The returned panel will have its
        :meth:`~SpectralPanel.generate_plot` called immediately
        afterwards.

        Parameters
        ----------
        xmin, xmax : float
            Wavelength bounds for the row.
        ax : Axes
            Pre-created axes that the panel should draw into.
        **kwargs
            Extra keyword arguments forwarded from :meth:`generate_plot`.

        Returns
        -------
        SpectralPanel
        """
        ...

    # ------------------------------------------------------------------
    # Optional per-row hook
    # ------------------------------------------------------------------
    def _post_render_panel(
        self,
        idx: int,
        panel: SpectralPanel,
        ax: Axes,
        is_last: bool,
    ) -> None:
        """Called after each panel is generated.

        Subclasses can override this to add per-row tick formatting,
        labels, annotations, chi-squared boxes, etc.  The default
        implementation applies common tick locators and hides x-axis
        labels on all but the last row.

        Parameters
        ----------
        idx : int
            Row index (0-based).
        panel : SpectralPanel
            The child panel that was just rendered.
        ax : Axes
            The axes the panel was drawn on.
        is_last : bool
            *True* when this is the bottom row.
        """
        fg = self._get_theme_value("foreground", "black")
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=8, prune="both"))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune="both"))
        if is_last:
            ax.set_xlabel("Wavelength (\u03bcm)", fontsize=8, color=fg)
        else:
            ax.tick_params(axis="x", labelbottom=False)

    # ------------------------------------------------------------------
    # y-limit computation (same algorithm as FSP._compute_panel_ylims)
    # ------------------------------------------------------------------
    @staticmethod
    def _default_ylim_fn(
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        mask: np.ndarray,
        ymax_factor: float = 0.2,
    ) -> Tuple[float, float]:
        """Compute default ``(ymin, ymax)`` from observed flux.

        Parameters
        ----------
        wave_data, flux_data : np.ndarray
            Full data arrays.
        mask : np.ndarray
            Boolean mask selecting points in the current panel.
        ymax_factor : float
            Fractional headroom.
        """
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
        """Compute per-panel y-limits (mirrors FSP._compute_panel_ylims).

        Parameters
        ----------
        uniform : bool, optional
            When *True* all entries share a global min/max.  Defaults to
            :attr:`uniform_ylim`.
        ylim_fn : callable, optional
            ``fn(mask) -> (ymin, ymax)``.  Defaults to a closure around
            :meth:`_default_ylim_fn` using this instance's data and
            *ymax_factor*.

        Returns
        -------
        list[tuple[float, float]]
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
    # Top-level plot generation
    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Build the stacked multi-panel figure.

        For each panel edge a new :class:`SpectralPanel` is created
        (via :meth:`_create_panel`), its :meth:`generate_plot` is
        invoked, and :meth:`_post_render_panel` is called for per-row
        decoration.
        """
        n = len(self._panel_edges)
        self._ensure_figure()
        self.fig.clf()
        self.panels.clear()
        self.subplots.clear()

        # Disable constrained_layout for uniform panel heights.
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

        # --- Pre-compute per-panel y-limits ----------------------------
        panel_ylims = self._compute_panel_ylims()

        fg = self._get_theme_value("foreground", "black")

        for idx, xlim_start in enumerate(self._panel_edges):
            is_last = idx == n - 1
            panel_end = xlim_start + self._step
            xmin, xmax = xlim_start, panel_end

            ax = self.fig.add_subplot(gs[idx, 0])
            self.subplots[idx] = ax

            # --- Delegate to the concrete factory ----------------------
            panel = self._create_panel(xmin, xmax, ax, **kwargs)
            self.panels[idx] = panel

            # Let the child panel render onto the supplied axes.
            panel.generate_plot(**kwargs)

            # Apply pre-computed y-limits
            ax.set_xlim(xmin, xmax)
            ymin, ymax = panel_ylims[idx]
            ax.set_ylim(ymin, ymax)

            # Per-row decoration hook
            self._post_render_panel(idx, panel, ax, is_last)

        # --- Global figure labels --------------------------------------
        self.fig.supylabel("Flux Density (Jy)", fontsize=10, color=fg, x=0.01)

        # Apply theme
        self.apply_theme_to_figure()
