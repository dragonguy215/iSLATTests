"""
CompositeStackedPanel -- Concrete :class:`StackedSpectralPanel` produced by
:meth:`~StackedSpectralPanel.stack_with` or the ``+`` operator.

Composites rows from two (or more) source stacked-spectral plots into a
single figure.  Each row is backed by a specific ``(source_plot, cell_index)``
pair; cell creation, rendering, and post-render decoration are all delegated
to the source plot that owns the row.

A composite of a :class:`FullSpectrumPlot` and a
:class:`ResidualSpectrumPlot` correctly produces the expected single-axes
and spectrum+residual rows, respectively.
"""
from __future__ import annotations

from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Any,
    TYPE_CHECKING,
)

import numpy as np
from matplotlib.gridspec import GridSpec

from .StackedSpectralPanel import StackedSpectralPanel, _extract_cell_kwargs
from .SpectralPanel import SpectralPanel, GapMode

if TYPE_CHECKING:
    from matplotlib.gridspec import SubplotSpec

class CompositeStackedPanel(StackedSpectralPanel):
    """Concrete :class:`StackedSpectralPanel` that composites rows from two (or more) source stacked-spectral plots into a single figure.

    Each row in the composite figure is backed by a specific ``(source_plot, cell_index)`` pair.
    Cell creation, rendering, and post-render decoration are all delegated to the source plot that
    owns the row, so a composite of a :class:`FullSpectrumPlot` and a :class:`ResidualSpectrumPlot` correctly produces the expected
    single-axes and spectrum+residual rows, respectively.

    Use :meth:`from_pair` (or the more convenient :meth:`StackedSpectralPanel.stack_with` / ``+`` operator) to build instances.

    Attributes
    ----------
    row_plan : list[tuple[StackedSpectralPanel, int]]
        Ordered sequence of ``(source_plot, cell_index)`` describing
        each row in the composite figure.
    labels : tuple[str, str] or None
        Optional descriptive labels for each source, rendered in the
        left margin of each row.
    sources : tuple[StackedSpectralPanel, StackedSpectralPanel]
        The two source plots that were combined.
    """
    def __init__(
        self,
        row_plan: List[Tuple[StackedSpectralPanel, int]],
        sources: Tuple[StackedSpectralPanel, StackedSpectralPanel],
        *,
        hspace: float = 0.25,
        row_height: float = 2.0,
        figsize: Optional[Tuple[float, float]] = None,
        labels: Optional[Tuple[str, str]] = None,
        **kwargs,
    ) -> None:
        # We need a dummy wave/flux to satisfy the parent __init__;
        # the composite never uses them for ylim / mask computation
        # because every cell is delegated to the real source plot.
        first_src = sources[0]
        super().__init__(
            wave_data=first_src.wave_data,
            flux_data=first_src.flux_data,
            error_data=first_src.error_data,
            molecules=getattr(first_src, "molecules", None),
            n_panels=len(row_plan),
            hspace=hspace,
            row_height=row_height,
            figsize=figsize,
            **kwargs,
        )

        self.row_plan = row_plan
        self.sources = sources
        self.labels = labels

        # Override the panel_edges computed by the parent - they are
        # meaningless for a composite; we only need the *count*.
        self._panel_edges = np.arange(len(row_plan), dtype=float)
        self._panel_ends = np.arange(1, len(row_plan) + 1, dtype=float)
        self._step = 1.0  # dummy

        # Pre-compute kwargs for each source once.
        self._source_kwargs: Dict[int, Dict[str, Any]] = {
            id(sources[0]): _extract_cell_kwargs(sources[0]),
            id(sources[1]): _extract_cell_kwargs(sources[1]),
        }

        # Figsize auto-compute
        if figsize is None:
            self._figsize = (14, row_height * len(row_plan))

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------
    @classmethod
    def from_pair(
        cls,
        plot_a: StackedSpectralPanel,
        plot_b: StackedSpectralPanel,
        *,
        hspace: float = 0.25,
        row_height: Optional[float] = None,
        figsize: Optional[Tuple[float, float]] = None,
        labels: Optional[Tuple[str, str]] = None,
    ) -> CompositeStackedPanel:
        """Build a composite from two stacked-spectral plots.

        Cells are matched by **closest wavelength-range midpoint**.
        Matched pairs are interleaved (A-row, B-row); unmatched cells
        from either source are appended at the end.
        """
        # Panel edges / ends are computed in __init__, so we do NOT
        # need to call generate_plot() here.  Doing so would create
        # pyplot-managed figures that leak as unwanted widget outputs
        # in Jupyter notebooks.

        # Midpoints for matching.
        mid_a = (plot_a._panel_edges + plot_a._panel_ends) / 2.0
        mid_b = (plot_b._panel_edges + plot_b._panel_ends) / 2.0

        # Greedy closest-midpoint matching (one-to-one).
        matched: List[Tuple[int, int]] = []
        used_b: set[int] = set()

        for a_idx in range(len(mid_a)):
            dists = np.abs(mid_b - mid_a[a_idx])
            # Iterate in order of increasing distance.
            for b_idx in np.argsort(dists):
                b_idx_int = int(b_idx)
                if b_idx_int not in used_b:
                    matched.append((a_idx, b_idx_int))
                    used_b.add(b_idx_int)
                    break

        unmatched_a = sorted(
            set(range(len(mid_a))) - {p[0] for p in matched}
        )
        unmatched_b = sorted(
            set(range(len(mid_b))) - used_b
        )

        # Build the row plan.
        row_plan: List[Tuple[StackedSpectralPanel, int]] = []
        for a_idx, b_idx in sorted(matched, key=lambda p: p[0]):
            row_plan.append((plot_a, a_idx))
            row_plan.append((plot_b, b_idx))
        for a_idx in unmatched_a:
            row_plan.append((plot_a, a_idx))
        for b_idx in unmatched_b:
            row_plan.append((plot_b, b_idx))

        rh = row_height or max(plot_a.row_height, plot_b.row_height)

        # Inherit gap_mode from sources if either uses SKIP.
        src_gap_mode = GapMode.CONNECT
        src_gap_threshold = None
        if plot_a.gap_mode is GapMode.SKIP or plot_b.gap_mode is GapMode.SKIP:
            src_gap_mode = GapMode.SKIP
            src_gap_threshold = (
                plot_a.gap_threshold
                if plot_a.gap_threshold is not None
                else plot_b.gap_threshold
            )

        return cls(
            row_plan=row_plan,
            sources=(plot_a, plot_b),
            hspace=hspace,
            row_height=rh,
            figsize=figsize,
            labels=labels,
            gap_mode=src_gap_mode,
            gap_threshold=src_gap_threshold,
        )

    # ------------------------------------------------------------------
    # Gap-skip override - delegate to the source that owns the row
    # ------------------------------------------------------------------
    def _cell_has_data(self, xmin: float, xmax: float) -> bool:
        """Delegate to the real source plot for the row at this position.

        The composite uses dummy panel_edges, so the base-class implementation (which checks ``self.wave_data``) would give wrong results.
        Instead we look up the source that owns each row and ask *it* whether the cell has data.
        """
        # xmin in the composite is a dummy integer (0, 1, 2, ...).
        # Map it back to the row_plan entry.
        row_idx = int(round(xmin))
        if row_idx < 0 or row_idx >= len(self.row_plan):
            return True  # Fallback: keep the row
        owner, cell_idx = self.row_plan[row_idx]
        edge = owner._panel_edges[cell_idx]
        return owner._cell_has_data(edge, owner._panel_ends[cell_idx])

    # ------------------------------------------------------------------
    # _create_cell - delegate to the source plot
    # ------------------------------------------------------------------
    def _create_cell(
        self,
        idx: int,
        xmin: float,
        xmax: float,
        gs_slot: "SubplotSpec",
        **kwargs,
    ) -> List[SpectralPanel]:
        """Delegate cell creation to the source plot that owns this row."""
        owner, cell_idx = self.row_plan[idx]
        edge = owner._panel_edges[cell_idx]
        real_xmin = edge
        real_xmax = owner._panel_ends[cell_idx]
        kw = self._source_kwargs.get(id(owner), {})
        return owner._create_cell(cell_idx, real_xmin, real_xmax, gs_slot, **kw)

    # ------------------------------------------------------------------
    # _post_render_cell - delegate to the source plot
    # ------------------------------------------------------------------
    def _post_render_cell(
        self,
        idx: int,
        cell_panels: List[SpectralPanel],
        is_last: bool,
    ) -> None:
        """Delegate post-render decoration to the source plot."""
        owner, cell_idx = self.row_plan[idx]
        owner._post_render_cell(cell_idx, cell_panels, is_last)

    # ------------------------------------------------------------------
    # generate_plot - override to handle fig redirection + auto-labels
    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:
        """Build the composite figure.

        Temporarily redirects each source plot's ``fig`` reference to
        the composite figure so that ``_create_cell`` (which calls
        ``self.fig.add_subplot``) renders into the correct figure.

        When *labels* are set, each panel's top-left corner is annotated with the corresponding source label.
        """
        n = len(self.row_plan)
        self._ensure_figure()
        self.fig.clf()
        self.panels.clear()

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

        # Temporarily redirect source figs.
        src_a, src_b = self.sources
        orig_fig_a = src_a.fig
        orig_fig_b = src_b.fig
        orig_panels_a = src_a.panels
        orig_panels_b = src_b.panels

        try:
            src_a.fig = self.fig
            src_b.fig = self.fig
            src_a.panels = {}
            src_b.panels = {}

            for row_idx in range(n):
                is_last = row_idx == n - 1
                cell_panels = self._create_cell(
                    row_idx, 0.0, 1.0, gs[row_idx, 0],
                )
                self.panels[row_idx] = cell_panels

                owner, cell_idx = self.row_plan[row_idx]
                kw = self._source_kwargs.get(id(owner), {})
                for panel in cell_panels:
                    panel.generate_plot(**kw)

                self._post_render_cell(row_idx, cell_panels, is_last)

        finally:
            src_a.fig = orig_fig_a
            src_b.fig = orig_fig_b
            src_a.panels = orig_panels_a
            src_b.panels = orig_panels_b

        # --- Auto-label panels with source names ----------------------
        # Build cell_labels dict from the labels tuple and row_plan.
        if self.labels is not None:
            label_a, label_b = self.labels
            cell_labels: Dict[int, str] = {}
            for row_idx, (owner, _) in enumerate(self.row_plan):
                cell_labels[row_idx] = (
                    label_a if owner is src_a else label_b
                )
            self._apply_cell_labels(cell_labels)

        # Global labels and theming.
        fg = self._get_theme_value("foreground", "black")
        self.fig.supylabel(
            "Flux Density (Jy)", fontsize=10, color=fg, x=0.01,
        )
        self.apply_theme_to_figure()