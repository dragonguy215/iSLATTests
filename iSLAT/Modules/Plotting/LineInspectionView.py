"""
LineInspectionView — :class:`SpectrumPanelView` for a standalone full-size
line inspection plot.

Shows a :class:`LineInspectionPlot` for the most recent span-selector
selection in a dedicated figure that fills the main plotting area.
When no selection has been made yet the view displays a placeholder.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, List, Optional, Tuple

import matplotlib.pyplot as plt

from .SpectrumPanelView import SpectrumPanelView
from .LineInspectionPlot import LineInspectionPlot
from .LineInspectionContextMixin import LineInspectionContextMixin

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine

try:
    from iSLAT.Modules.Debug import debug_config
except ImportError:
    class _FallbackDebug:
        def verbose(self, *a, **kw): pass
        def info(self, *a, **kw): pass
        def warning(self, *a, **kw): print(f"WARNING: {kw or a}")
        def error(self, *a, **kw): print(f"ERROR: {kw or a}")
        def trace(self, *a, **kw): pass
    debug_config = _FallbackDebug()


class LineInspectionView(SpectrumPanelView, LineInspectionContextMixin):
    """
    Standalone line-inspection view.

    Renders a :class:`LineInspectionPlot` for the currently selected
    wavelength range (taken from the main controller's
    ``toggle_state["current_selection"]``).  The view owns its own
    :class:`~matplotlib.figure.Figure` and
    :class:`~matplotlib.backends.backend_tkagg.FigureCanvasTkAgg` so it
    can be swapped in / out of the main window independently of the
    three-panel and full-spectrum views.

    When no selection has been made a placeholder message is shown.

    Parameters
    ----------
    plot_manager : iSLATPlot
        The main controller (``iSLAT.Modules.Plotting.MainPlot.iSLATPlot``).
    """

    def __init__(self, plot_manager: Any) -> None:
        super().__init__(plot_manager)

        # Inspection-specific state
        self._plot: Optional[LineInspectionPlot] = None
        self._current_selection: Optional[Tuple[float, float]] = None

    # ==================================================================
    # Internal helpers
    # ==================================================================

    def _get_active_molecule(self) -> Optional["Molecule"]:
        return getattr(self._islat, "active_molecule", None)

    def _get_molecules(self) -> Optional["MoleculeDict"]:
        return getattr(self._islat, "molecules_dict", None)

    def _build_plot(self) -> None:
        """Create (or recreate) the owned figure and :class:`LineInspectionPlot`."""
        # Close the previous figure to avoid memory leaks
        if self._fig is not None:
            plt.close(self._fig)
            self._fig = None
            self._plot = None

        self._fig = plt.Figure(figsize=(10, 5), constrained_layout=True)
        ax = self._fig.add_subplot(111)

        sel = self._current_selection
        if sel is None:
            # No selection yet — draw a placeholder
            ax.set_title("Line Inspection")
            ax.set_xlabel("Wavelength (µm)")
            ax.set_ylabel("Flux density (Jy)")
            ax.text(
                0.5, 0.5,
                "Drag to select a wavelength range\nin the main spectrum",
                ha="center", va="center",
                transform=ax.transAxes,
                color="gray", fontsize=12,
            )
            return

        xmin, xmax = sel
        wave = getattr(self._islat, "wave_data", None)
        flux = getattr(self._islat, "flux_data", None)
        if wave is None or flux is None or len(wave) == 0:
            ax.set_title("Line Inspection")
            ax.set_xlabel("Wavelength (µm)")
            ax.set_ylabel("Flux density (Jy)")
            ax.text(
                0.5, 0.5, "No spectral data available",
                ha="center", va="center",
                transform=ax.transAxes,
                color="gray", fontsize=12,
            )
            return

        error = getattr(self._islat, "error_data", None)
        mol = self._get_active_molecule()
        mols = self._get_molecules()

        # Gather line data for the active molecule (same helper used by ThreePanelView)
        line_data: Optional[List] = None
        if mol is not None and hasattr(self._pm, "get_molecule_line_data"):
            try:
                line_data = self._pm.get_molecule_line_data(mol, xmin, xmax)
            except Exception as exc:
                debug_config.warning(
                    "line_inspection_view",
                    f"get_molecule_line_data failed: {exc}",
                )

        try:
            self._plot = LineInspectionPlot(
                wave_data=wave,
                flux_data=flux,
                xmin=xmin,
                xmax=xmax,
                error_data=error,
                molecule=mol,
                molecules=mols,
                line_data=line_data,
                ax=ax,
                render_all_visible=True,
                theme=self._pm.theme,
            )
            self._plot.generate_plot()
        except Exception as exc:
            debug_config.warning(
                "line_inspection_view",
                f"LineInspectionPlot.generate_plot failed: {exc}",
            )

    # ==================================================================
    # PlotView — lifecycle
    # ==================================================================

    def activate(self, parent_frame: Any) -> None:
        """Pack the line-inspection canvas into *parent_frame*."""
        self._parent_frame = parent_frame

        # Sync selection from the controller's shared state
        ctrl_sel = self._pm.toggle_state.get("current_selection")
        if ctrl_sel != self._current_selection:
            self._current_selection = ctrl_sel
            self._needs_refresh = True

        if not self._initialised or self._needs_refresh:
            self._build_plot()
            self._initialised = True
            self._needs_refresh = False

        self._ensure_canvas()
        self.apply_theme(self._pm.theme)
        self._register_control_fields()

        if self._canvas is not None:
            self._canvas.get_tk_widget().pack(
                fill="both", expand=True, padx=0, pady=0
            )
            self._canvas.draw_idle()

    # ==================================================================
    # PlotView — interaction context
    # ==================================================================

    def build_context_menu(self, event: Any, canvas_widget: Any) -> Any:
        """Return the line-inspection right-click menu."""
        return self._build_line_inspection_menu(canvas_widget)

    # ==================================================================
    # PlotView — core rendering
    # ==================================================================

    def update_model_plot(
        self,
        wave_data: Any = None,
        flux_data: Any = None,
        molecules_dict: "MoleculeDict" = None,
        error_data: Optional[Any] = None,
        **kwargs: Any,
    ) -> None:
        """Rebuild the inspection plot from current iSLAT data."""
        self._needs_refresh = True
        if self._initialised:
            self._build_plot()
            self._needs_refresh = False
            self._ensure_canvas()
            self.apply_theme(self._pm.theme)
            if self._canvas is not None:
                self._canvas.draw_idle()

    # ==================================================================
    # PlotView — selection (override base no-ops)
    # ==================================================================

    def on_selection(self, xmin: float, xmax: float) -> None:
        """Update the displayed wavelength range when the user makes a new selection."""
        self._current_selection = (xmin, xmax)
        self._needs_refresh = True
        if self._initialised:
            self._build_plot()
            self._needs_refresh = False
            self._ensure_canvas()
            self.apply_theme(self._pm.theme)
            if self._canvas is not None:
                self._canvas.draw_idle()

    def clear_selection(self) -> None:
        """Clear the selection and show the placeholder."""
        self._current_selection = None
        self._needs_refresh = True
        if self._initialised:
            self._build_plot()
            self._needs_refresh = False
            self._ensure_canvas()
            if self._canvas is not None:
                self._canvas.draw_idle()

    # ==================================================================
    # PlotView — molecule lifecycle (override to preserve selection)
    # ==================================================================

    def on_active_molecule_changed(
        self,
        new_molecule: Optional["Molecule"] = None,
        current_selection: Optional[Tuple[float, float]] = None,
    ) -> None:
        """Rebuild for the newly selected active molecule."""
        if current_selection is not None:
            self._current_selection = current_selection
        self.update_model_plot()

    def get_selected_line(self) -> Optional[Any]:
        """No persistent selected-line state in this view."""
        return None