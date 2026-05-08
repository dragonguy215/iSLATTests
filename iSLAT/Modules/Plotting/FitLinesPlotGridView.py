"""
FitLinesPlotGridView — :class:`PlotView` for displaying a
:class:`FitLinesPlotGrid` within the main iSLAT plotting area.

Shows a grid of individual line-fit results for the current fitting
session.  When multiple grids are available (e.g. after a batch fit)
they are presented in a :class:`ttk.Notebook` with one tab per grid,
mirroring the layout of
:class:`~iSLAT.Modules.GUI.PlotGridWindow.PlotGridWindow`.  Each tab is
scrollable to accommodate grids that exceed the frame height.

When no plot grids have been set a placeholder message is shown.

Typical usage::

    # inside TopBar / fitting callback
    plot_grid = FitLinesPlotGrid(fit_data=fit_data, ...)
    plot_grid.generate_plot()
    main_plot.fit_lines_grid_view.set_plot_grids([plot_grid])
    main_plot.switch_view("Fit Lines Grid")
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt

from .PlotView import PlotView
from .FitLinesPlotGrid import FitLinesPlotGrid

if TYPE_CHECKING:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.Molecule import Molecule

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

try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    import tkinter as tk
    from tkinter import ttk
except ImportError:
    FigureCanvasTkAgg = None  # type: ignore[assignment,misc]
    tk = None  # type: ignore[assignment]
    ttk = None  # type: ignore[assignment]

# Size constants mirrored from PlotGridWindow
SUBPLOT_WIDTH_INCHES = 1.8
SUBPLOT_HEIGHT_INCHES = 1.5


class FitLinesPlotGridView(PlotView):
    """
    Embedded fit-lines plot grid view.

    Renders one or more :class:`FitLinesPlotGrid` objects inside the main
    iSLAT plotting area.  When multiple grids are available they are
    presented in a :class:`ttk.Notebook` with one tab per grid,
    mirroring the layout of
    :class:`~iSLAT.Modules.GUI.PlotGridWindow.PlotGridWindow`.  Each tab
    is scrollable to accommodate grids that exceed the frame height.

    When no plot grids have been set a placeholder message is displayed.

    Call :meth:`set_plot_grids` to push a new list of
    :class:`FitLinesPlotGrid` objects into the view; the view will
    automatically refresh if it is currently active.

    Parameters
    ----------
    plot_manager : iSLATPlot
        The main controller (``iSLAT.Modules.Plotting.MainPlot.iSLATPlot``).
    """

    def __init__(self, plot_manager: Any) -> None:
        self._pm = plot_manager
        self._islat = plot_manager.islat

        self._parent_frame: Any = None

        # Primary canvas / figure references (point at first tab or placeholder)
        self._canvas: Optional["FigureCanvasTkAgg"] = None
        self._fig: Optional["Figure"] = None

        # Placeholder figure (shown when _plot_grids is empty)
        self._placeholder_fig: Optional["Figure"] = None
        self._placeholder_canvas: Optional["FigureCanvasTkAgg"] = None

        # Notebook and per-tab canvases for multi-grid display
        self._notebook: Optional[Any] = None            # ttk.Notebook
        self._tab_canvases: List["FigureCanvasTkAgg"] = []
        self._container_frame: Optional[Any] = None    # ttk.Frame packed into parent

        # The list of FitLinesPlotGrid objects to display
        self._plot_grids: List[FitLinesPlotGrid] = []

        self._initialised: bool = False
        self._needs_refresh: bool = True

    # ==================================================================
    # Public API
    # ==================================================================

    def set_plot_grids(self, plot_grids: List[FitLinesPlotGrid]) -> None:
        """Set the list of :class:`FitLinesPlotGrid` objects to display.

        Calling this method marks the view as needing a refresh.  If the
        view is already active the display is updated immediately.

        Parameters
        ----------
        plot_grids : list of FitLinesPlotGrid
            The new grids to display.  Pass an empty list (or ``None``) to
            show the placeholder.
        """
        self._plot_grids = list(plot_grids) if plot_grids else []
        self._needs_refresh = True
        if self._initialised and self._parent_frame is not None:
            self._rebuild_display()

    # ==================================================================
    # Internal helpers
    # ==================================================================

    def _build_line_list_preview(self) -> Optional[FitLinesPlotGrid]:
        """Build a :class:`FitLinesPlotGrid` from the loaded line list with no fit data.

        Returns a grid whose panels show the observed spectrum around every
        line in ``islat.input_line_list`` but with *no* Gaussian overlay —
        identical to the post-fit grid except fit results are absent.

        Returns ``None`` if the required data (line list file, wave/flux
        arrays) is unavailable.
        """
        try:
            import pandas as pd
        except ImportError:
            return None

        line_list_path = getattr(self._islat, "input_line_list", None)
        if not line_list_path:
            return None

        wave = getattr(self._islat, "wave_data", None)
        flux = getattr(self._islat, "flux_data", None)
        if wave is None or flux is None or len(wave) == 0:
            return None

        try:
            df = pd.read_csv(line_list_path)
        except Exception as exc:
            debug_config.warning(
                "fit_lines_grid_view",
                f"Could not read line list for preview: {exc}",
            )
            return None

        if df.empty:
            return None

        # Build fit_csv_dict from the line list rows.  Only the keys consumed
        # by FitLinesPlotGrid._create_panels and _post_render are required:
        # xmin, xmax, species, lam, Fit_det.
        fit_csv_dict: dict = {}
        for idx, row in enumerate(df.itertuples(index=False)):
            xmin = getattr(row, "xmin", None)
            xmax = getattr(row, "xmax", None)
            if xmin is None or xmax is None:
                # Fall back to a ±0.015 µm window around lam if bounds missing
                lam = getattr(row, "lam", None)
                if lam is None:
                    continue
                xmin = float(lam) - 0.015
                xmax = float(lam) + 0.015
            fit_csv_dict[idx] = {
                "xmin": float(xmin),
                "xmax": float(xmax),
                "species": str(getattr(row, "species", "")),
                "lam": float(getattr(row, "lam", 0.0)),
                "Fit_det": False,
            }

        if not fit_csv_dict:
            return None

        n = len(fit_csv_dict)
        # Null fit tuples — _post_render handles None gracefully
        null_fits = [None] * n
        fit_data_tuple_list = (null_fits, null_fits, null_fits)

        spectrum_name = getattr(
            self._islat, "loaded_spectrum_name",
            "Line List Preview",
        )

        try:
            err = getattr(self._islat, "err_data", None)
            grid = FitLinesPlotGrid(
                fit_data=(fit_csv_dict, fit_data_tuple_list),
                wave_data=wave,
                flux_data=flux,
                err_data=err,
                spectrum_name=spectrum_name,
            )
            grid.generate_plot()
            return grid
        except Exception as exc:
            debug_config.warning(
                "fit_lines_grid_view",
                f"Failed to build line list preview grid: {exc}",
            )
            return None

    def _build_placeholder(self) -> None:
        """Create a simple placeholder figure when no grid data is available."""
        if self._placeholder_fig is not None:
            plt.close(self._placeholder_fig)

        self._placeholder_fig = plt.Figure(figsize=(8, 5), constrained_layout=True)
        ax = self._placeholder_fig.add_subplot(111)
        ax.set_title("Fit Lines Plot Grid", fontsize=14)
        ax.text(
            0.5, 0.5,
            "No fit data available.\n\nRun 'Fit Lines' to populate this view.",
            ha="center", va="center",
            transform=ax.transAxes,
            color="gray", fontsize=12,
        )
        ax.axis("off")

    def _apply_theme_to_fig(self, fig: "Figure") -> None:
        """Apply the controller's current theme to *fig*."""
        theme = self._pm.theme
        bg = theme.get("background", "#181A1B")
        fg = theme.get("foreground", "#F0F0F0")

        fig.patch.set_facecolor(bg)
        for ax in fig.axes:
            ax.set_facecolor(theme.get("axes_background", bg))
            ax.tick_params(colors=fg)
            ax.xaxis.label.set_color(fg)
            ax.yaxis.label.set_color(fg)
            ax.title.set_color(fg)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg)

    def _teardown_display(self) -> None:
        """Destroy all Tk widgets owned by this view."""
        for canvas in self._tab_canvases:
            try:
                canvas.get_tk_widget().destroy()
            except Exception:
                pass
        self._tab_canvases.clear()

        if self._placeholder_canvas is not None:
            try:
                self._placeholder_canvas.get_tk_widget().destroy()
            except Exception:
                pass
            self._placeholder_canvas = None

        # Clear primary canvas reference (widget already destroyed above)
        self._canvas = None
        self._fig = None

        if self._notebook is not None:
            try:
                self._notebook.destroy()
            except Exception:
                pass
            self._notebook = None

        if self._container_frame is not None:
            try:
                self._container_frame.destroy()
            except Exception:
                pass
            self._container_frame = None

    def _rebuild_display(self) -> None:
        """Tear down and recreate all display widgets for the current plot grids."""
        self._teardown_display()

        if self._parent_frame is None or tk is None:
            return

        if not self._plot_grids:
            self._show_placeholder()
            return

        # Outer container frame packed into the parent
        self._container_frame = ttk.Frame(self._parent_frame)
        self._container_frame.pack(fill="both", expand=True)

        if len(self._plot_grids) == 1:
            # Single grid — no notebook overhead, just one scrollable canvas
            self._create_scrollable_tab(self._container_frame, self._plot_grids[0])
        else:
            # Multiple grids — tabbed notebook, one tab per grid
            self._notebook = ttk.Notebook(self._container_frame)
            self._notebook.pack(fill="both", expand=True)
            for plot_grid in self._plot_grids:
                tab_frame = ttk.Frame(self._notebook)
                self._notebook.add(tab_frame, text=plot_grid.spectrum_name)
                self._create_scrollable_tab(tab_frame, plot_grid)

    def _show_placeholder(self) -> None:
        """Display the placeholder figure inside the parent frame.

        If a line list is loaded but no fit has been run, automatically
        renders a preview grid instead of a blank placeholder.
        """
        # Try to build a line-list preview before falling back to blank
        preview = self._build_line_list_preview()
        if preview is not None:
            # Temporarily set as the only grid and render it
            self._plot_grids = [preview]
            self._rebuild_display()
            return

        self._build_placeholder()
        self._apply_theme_to_fig(self._placeholder_fig)

        if FigureCanvasTkAgg is not None and self._parent_frame is not None:
            self._placeholder_canvas = FigureCanvasTkAgg(
                self._placeholder_fig, master=self._parent_frame
            )
            self._placeholder_canvas.draw()
            self._placeholder_canvas.get_tk_widget().pack(fill="both", expand=True)

        # Expose as primary canvas / figure for get_canvas / get_figure
        self._canvas = self._placeholder_canvas
        self._fig = self._placeholder_fig

    def _create_scrollable_tab(
        self,
        parent: Any,
        plot_grid: FitLinesPlotGrid,
    ) -> None:
        """Build a scrollable matplotlib canvas for *plot_grid* inside *parent*.

        The layout mirrors
        :meth:`~iSLAT.Modules.GUI.PlotGridWindow.PlotGridWindow._create_tab`.

        Parameters
        ----------
        parent : Tk widget
            The frame (tab or container) to build the scrollable view in.
        plot_grid : FitLinesPlotGrid
            The grid whose figure to embed.
        """
        # --- Save toolbar ---------------------------------------------------
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill="x", side="top", pady=2)

        btn_save = ttk.Button(
            toolbar_frame, text="Save Figure",
            command=lambda pg=plot_grid: self._save_figure(pg),
        )
        btn_save.pack(side="left", padx=5)

        # --- Scrollable area ------------------------------------------------
        canvas_frame = ttk.Frame(parent)
        canvas_frame.pack(fill="both", expand=True)

        scrollbar = ttk.Scrollbar(canvas_frame, orient="vertical")
        scrollbar.pack(side="right", fill="y")

        scroll_canvas = tk.Canvas(
            canvas_frame, highlightthickness=0,
            yscrollcommand=scrollbar.set,
        )
        scroll_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=scroll_canvas.yview)

        inner_frame = ttk.Frame(scroll_canvas)
        canvas_window = scroll_canvas.create_window(
            (0, 0), window=inner_frame, anchor="nw"
        )

        # --- Configure figure size ------------------------------------------
        fig_width = SUBPLOT_WIDTH_INCHES * plot_grid.cols
        fig_height = SUBPLOT_HEIGHT_INCHES * plot_grid.rows
        plot_grid.fig.set_size_inches(fig_width, fig_height)
        plot_grid.fig.set_dpi(150)
        # Disable any active layout engine before calling subplots_adjust,
        # which is incompatible with constrained_layout / tight_layout.
        # Using None (not the string "none") fully removes the engine.
        try:
            plot_grid.fig.set_layout_engine(None)
        except Exception:
            try:
                plot_grid.fig.set_constrained_layout(False)
            except Exception:
                pass
            try:
                plot_grid.fig.set_tight_layout(False)
            except Exception:
                pass
        plot_grid.fig.subplots_adjust(
            left=0.06, right=0.98,
            top=0.96, bottom=0.05,
            wspace=0.35, hspace=0.65,
        )

        # Compact font sizes for dense subplot grids
        for ax in plot_grid.axs.flat:
            ax.tick_params(axis="both", labelsize=6, pad=1)
            ax.title.set_fontsize(7)
            ax.title.set_position((0.5, 1.0))
            if ax.yaxis.label:
                ax.yaxis.label.set_fontsize(6)

        self._apply_theme_to_fig(plot_grid.fig)

        # --- Embed matplotlib figure ----------------------------------------
        fig_canvas = FigureCanvasTkAgg(plot_grid.fig, master=inner_frame)
        fig_canvas.draw()
        fig_canvas.get_tk_widget().pack(fill="both", expand=True)
        self._tab_canvases.append(fig_canvas)

        # First tab's canvas / figure become the primary references
        if self._canvas is None:
            self._canvas = fig_canvas
            self._fig = plot_grid.fig

        # --- Scroll behaviour -----------------------------------------------
        def update_scroll_region(event=None):
            scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

        inner_frame.bind("<Configure>", update_scroll_region)

        def configure_inner_frame(event):
            scroll_canvas.itemconfig(canvas_window, width=event.width)

        scroll_canvas.bind("<Configure>", configure_inner_frame)

        def on_mousewheel(event):
            scroll_canvas.yview_scroll(int(-1 * (event.delta / 60)), "units")

        def bind_mousewheel(event):
            scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)

        def unbind_mousewheel(event):
            scroll_canvas.unbind_all("<MouseWheel>")

        scroll_canvas.bind("<Enter>", bind_mousewheel)
        scroll_canvas.bind("<Leave>", unbind_mousewheel)

    def _save_figure(self, plot_grid: FitLinesPlotGrid) -> None:
        """Open a save-as dialog and write *plot_grid*'s figure to disk."""
        from tkinter import filedialog

        filetypes = [
            ("PNG files", "*.png"),
            ("PDF files", "*.pdf"),
            ("SVG files", "*.svg"),
            ("All files", "*.*"),
        ]
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=filetypes,
            initialfile=f"{plot_grid.spectrum_name}_fit_grid",
        )
        if filepath:
            plot_grid.fig.savefig(filepath, dpi=300, bbox_inches="tight")

    # ==================================================================
    # PlotView — lifecycle
    # ==================================================================

    def activate(self, parent_frame: Any) -> None:
        """Pack the fit-lines grid view into *parent_frame* and refresh if needed.

        If no plot grids have been set but a line list is loaded, a preview
        grid is built automatically from the line list so the view is never
        shown empty when data are available.
        """
        self._parent_frame = parent_frame

        # When no explicit grids are present, always re-attempt a preview
        # (the line list or spectrum data may have changed since last activation).
        if not self._plot_grids:
            self._needs_refresh = True

        if not self._initialised or self._needs_refresh:
            self._rebuild_display()
            self._initialised = True
            self._needs_refresh = False

        self.apply_theme(self._pm.theme)
        self._register_control_fields()

    def deactivate(self) -> None:
        """Tear down the view's widgets and unregister any ControlBus fields."""
        bus = getattr(self._pm, "control_bus", None)
        if bus is not None:
            bus.unregister_owner(self)

        self._teardown_display()

    # ==================================================================
    # PlotView — interaction context
    # ==================================================================

    def build_context_menu(self, event: Any, canvas_widget: Any) -> Any:
        """No custom right-click menu for this view."""
        return None

    # ==================================================================
    # PlotView — theme
    # ==================================================================

    def apply_theme(self, theme: dict) -> None:
        """Apply the current theme to all embedded figures and canvases."""
        self._pm.theme = theme

        for plot_grid in self._plot_grids:
            try:
                self._apply_theme_to_fig(plot_grid.fig)
            except Exception as exc:
                debug_config.warning("fit_lines_grid_view", f"apply_theme error: {exc}")

        if not self._plot_grids and self._placeholder_fig is not None:
            try:
                self._apply_theme_to_fig(self._placeholder_fig)
            except Exception:
                pass

        for canvas in self._tab_canvases:
            try:
                canvas.draw_idle()
            except Exception:
                pass

        if self._canvas is not None and self._canvas not in self._tab_canvases:
            try:
                self._canvas.draw_idle()
            except Exception:
                pass

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
        """No-op — fit-lines grids are self-contained and do not re-render on model changes.

        Use :meth:`set_plot_grids` to push new grids into this view.
        """
        pass

    def on_molecule_visibility_changed(
        self,
        molecule_name: str,
        is_visible: bool,
        molecules_dict: "MoleculeDict",
        wave_data: Any,
        active_molecule: Optional["Molecule"] = None,
        current_selection: Optional[Tuple[float, float]] = None,
        force_rerender: bool = False,
    ) -> None:
        """No-op — fit-lines grids do not react to molecule visibility changes."""
        pass

    # ==================================================================
    # PlotView — selection / line inspection (no-ops)
    # ==================================================================

    def on_selection(self, xmin: float, xmax: float) -> None:
        """No-op — this view does not respond to span-selector events."""
        pass

    def clear_selection(self) -> None:
        """No-op."""
        pass

    def clear_active_lines(self) -> None:
        """No-op — no active-line markers in this view."""
        pass

    # ==================================================================
    # PlotView — molecule lifecycle (no-ops)
    # ==================================================================

    def on_active_molecule_changed(
        self,
        new_molecule: Optional["Molecule"] = None,
        current_selection: Optional[Tuple[float, float]] = None,
    ) -> None:
        """No-op."""
        pass

    def on_molecule_parameter_changed(
        self,
        molecule_name: str,
        parameter_name: str,
        current_selection: Optional[Tuple[float, float]] = None,
    ) -> None:
        """No-op."""
        pass

    def on_molecule_deleted(self, molecule_name: str) -> None:
        """No-op."""
        pass

    # ==================================================================
    # PlotView — toggle helpers (all no-ops for this view)
    # ==================================================================

    def sync_toggle_state(self, toggle_state: dict) -> None:
        """No-op — no toggleable overlays in this view."""
        pass

    def toggle_summed_spectrum(self, visible: bool) -> None:
        """No-op."""
        pass

    def toggle_legend(self, visible: Optional[bool] = None) -> None:
        """No-op."""
        pass

    def toggle_saved_lines(self, show: bool, loaded_lines: Any = None) -> None:
        """No-op."""
        pass

    def toggle_atomic_lines(self, show: bool) -> None:
        """No-op."""
        pass

    # ==================================================================
    # PlotView — canvas / drawing
    # ==================================================================

    def draw(self) -> None:
        """Flush all pending changes to every embedded canvas."""
        for canvas in self._tab_canvases:
            try:
                canvas.draw_idle()
            except Exception:
                pass
        # Also flush the placeholder canvas if it's active
        if self._canvas is not None and self._canvas not in self._tab_canvases:
            try:
                self._canvas.draw_idle()
            except Exception:
                pass

    def get_canvas(self) -> "FigureCanvasTkAgg":
        """Return the primary :class:`FigureCanvasTkAgg` for this view."""
        return self._canvas

    def get_figure(self) -> "Figure":
        """Return the primary :class:`~matplotlib.figure.Figure` for this view."""
        return self._fig
