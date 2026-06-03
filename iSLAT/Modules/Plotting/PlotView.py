"""
PlotView - abstract interface for switchable plot views in the iSLAT GUI.

The iSLATPlot controller owns one *active_view* at a time.

Every user-facing action (toggle molecule, toggle summed spectrum, toggle legend, …) is forwarded to the active view's implementation, eliminating scattered ``if is_full_spectrum`` checks throughout the codebase.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.Molecule import Molecule

class PlotView(ABC):
    """
    Abstract base for a swappable plot view inside the main iSLAT window.

    Each view owns a matplotlib *Figure* and a Tk *FigureCanvasTkAgg*.
    The controller (:class:`iSLATPlot`) calls these methods without knowing which concrete view is active.
    """
    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------
    @abstractmethod
    def activate(self, parent_frame: Any) -> None:
        """
        Make this view the visible one.

        Pack / display the view's canvas inside *parent_frame* and ensure the rendered content is up-to-date.
        """
        ...

    @abstractmethod
    def deactivate(self) -> None:
        """
        Hide this view (pack_forget the canvas).

        The view should **not** destroy its resources - it may be reactivated later.
        """
        ...

    # ------------------------------------------------------------------
    # Core rendering
    # ------------------------------------------------------------------
    @abstractmethod
    def update_model_plot(
        self,
        wave_data: Any,
        flux_data: Any,
        molecules_dict: "MoleculeDict",
        error_data: Optional[Any] = None,
        **kwargs: Any,
    ) -> None:
        """
        Full re-render of the model spectrum (observed + molecules + sum).

        Called when a new spectrum is loaded, molecule parameters change globally, or the user switches into this view.
        """
        ...

    @abstractmethod
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
        """
        Lightweight update after a single molecule's visibility is toggled.

        Implementations should **not** reload data from disk - only update the artists that changed.

        Parameters
        ----------
        force_rerender : bool
            When *True* the molecule's artists must be re-rendered from current parameters (e.g. because parameters changed while the molecule was hidden).
            Implementations should remove the stale artists and create fresh ones instead of just toggling visibility.
        """
        ...

    # ------------------------------------------------------------------
    # Selection & line-inspection
    # ------------------------------------------------------------------
    @abstractmethod
    def on_selection(self, xmin: float, xmax: float) -> None:
        """Handle a wavelength-range selection (span selector drag).

        The view should render the line inspection panel, populate the population diagram scatter points, and highlight the strongest line.
        Views that do not support interactive selection (e.g. :class:`FullSpectrumView`) should implement this as a no-op or as a deferred switch back to the three-panel view.
        """
        ...

    @abstractmethod
    def clear_selection(self) -> None:
        """Clear the current line-inspection selection and reset panels."""
        ...

    @abstractmethod
    def clear_active_lines(self) -> None:
        """Remove all active-line artists (vlines, text, scatter) from the view."""
        ...

    # ------------------------------------------------------------------
    # Molecule lifecycle callbacks
    # ------------------------------------------------------------------
    @abstractmethod
    def on_active_molecule_changed(
        self,
        new_molecule: Optional["Molecule"] = None,
        current_selection: Optional[Tuple[float, float]] = None,
    ) -> None:
        """The user selected a different active molecule.

        If a selection is active the view should re-run the line inspection / population diagram for the new molecule.
        Otherwise only the population-diagram title and base diagram need updating.
        """
        ...

    @abstractmethod
    def on_molecule_parameter_changed(
        self,
        molecule_name: str,
        parameter_name: str,
        current_selection: Optional[Tuple[float, float]] = None,
    ) -> None:
        """A single molecule's parameter changed (not visibility).

        The view should update the spectrum if the molecule is visible, and refresh the line-inspection / population diagram if the molecule is the active one.
        """
        ...

    @abstractmethod
    def on_molecule_deleted(self, molecule_name: str) -> None:
        """A molecule was removed from the dict - update all panels."""
        ...

    # ------------------------------------------------------------------
    # Theme
    # ------------------------------------------------------------------
    @abstractmethod
    def apply_theme(self, theme: dict) -> None:
        """
        Apply a new theme dictionary to this view.

        Implementations must propagate the theme to every owned figure, axes, canvas widget, and sub-plot delegate so that switching between views always reflects the current theme.

        Called by the controller's :meth:`iSLATPlot.apply_theme` and automatically on :meth:`activate` when the theme has changed since the view was last visible.
        """
        ...

    # ------------------------------------------------------------------
    # Toggle helpers
    # ------------------------------------------------------------------
    @abstractmethod
    def sync_toggle_state(self, toggle_state: dict) -> None:
        """
        Reconcile the view's visual state with the canonical *toggle_state* dict from the controller.

        Called by the controller when this view is **activated** so that
        overlays (atomic lines, saved lines, summed spectrum, legend) match
        the state the user set while a different view was visible.
        """
        ...

    @abstractmethod
    def toggle_summed_spectrum(self, visible: bool) -> None:
        """Show or hide the summed model fill across all relevant axes."""
        ...

    @abstractmethod
    def toggle_legend(self, visible: Optional[bool] = None) -> None:
        """Toggle the legend visibility on all relevant axes."""
        ...

    @abstractmethod
    def toggle_saved_lines(self, show: bool, loaded_lines: Any = None) -> None:
        """Add or remove saved line annotations."""
        ...

    @abstractmethod
    def toggle_atomic_lines(self, show: bool) -> None:
        """Add or remove atomic line annotations."""
        ...

    # ------------------------------------------------------------------
    # File output
    # ------------------------------------------------------------------
    def save_figure(
        self,
        save_path: str | None = None,
        file_format: str = "pdf",
        dpi: int | None = None,
        rasterized: bool = False,
        **kwargs,
    ) -> str | None:
        """
        Save the current view's figure to a file.

        The default implementation saves the figure returned by
        :meth:`get_figure`.
        Subclasses may override this to produce a *different* figure for export (e.g. with toggle state baked in).

        Parameters
        ----------
        save_path : str or None
            Destination path.  If *None* a file dialog is opened.
        file_format : str
            File format extension (``"pdf"``, ``"png"``, …).
        dpi : int or None
            Resolution.  *None* uses matplotlib's default.
        rasterized : bool
            If *True* axes are rasterized before saving (useful for PDFs
            with very dense data).
        **kwargs
            Extra keyword arguments forwarded to ``fig.savefig()``.

        Returns
        -------
        str or None
            The path that was saved to, or *None* if the user cancelled.
        """
        from pathlib import Path
        from tkinter import filedialog

        fig = self.get_figure()
        if fig is None:
            return None

        if save_path is None:
            save_path = filedialog.asksaveasfilename(
                title="Save Figure",
                defaultextension=f".{file_format}",
                filetypes=[(f"{file_format.upper()} files", f"*.{file_format}")],
            )
        if not save_path:
            return None

        if rasterized:
            for ax in fig.axes:
                ax.set_rasterized(True)

        save_kw = {"bbox_inches": "tight", "format": file_format}
        if dpi is not None:
            save_kw["dpi"] = dpi
        elif rasterized:
            save_kw["dpi"] = 300
        else:
            save_kw["dpi"] = 300  # high-quality default for all exports
        save_kw.update(kwargs)

        fig.savefig(save_path, **save_kw)
        return save_path

    # ------------------------------------------------------------------
    # Line-inspection helpers
    # ------------------------------------------------------------------
    def get_selected_line(self) -> Optional[Dict[str, Any]]:
        """Return the currently selected line info dict, or *None*.

        Views that track a selected-line (e.g. :class:`ThreePanelView`)
        should override this.  The default returns *None*.
        """
        return None

    # ------------------------------------------------------------------
    # Interaction context
    # ------------------------------------------------------------------
    def build_context_menu(self, event: Any, canvas_widget: Any) -> Optional[Any]:
        """Build and return a ``tk.Menu`` for a right-click at *event*.

        The default creates a minimal menu containing only "Save Figure".
        Views that provide additional items should create their own menu,
        append items, and call :meth:`_append_save_figure_item` before
        returning.

        Parameters
        ----------
        event :
            The matplotlib ``MouseEvent`` that triggered the right-click.
        canvas_widget :
            The Tk widget backing the active canvas (used to anchor the menu
            and compute screen coordinates).
        """
        try:
            import tkinter as tk
        except ImportError:
            return None
        menu = tk.Menu(canvas_widget, tearoff=0)
        self._append_save_figure_item(menu)
        return menu

    def _append_save_figure_item(self, menu: Any) -> None:
        """Append a ``Save Figure…`` command to *menu*.

        Calls :meth:`save_figure` with a file-dialog so the user can choose
        the destination and format.
        """
        menu.add_separator()
        menu.add_command(label="Save Figure…", command=self.save_figure)

    def _register_control_fields(self) -> None:
        """Register this view's :class:`ControlField` objects on the :class:`ControlBus`.

        Called from :meth:`activate`.  The default is a no-op.  Views that
        expose view-specific controls (e.g. :class:`FullSpectrumView`'s
        n-panels entry) should override this and call
        ``bus.unregister_owner(self)`` in :meth:`deactivate`.
        """

    # ------------------------------------------------------------------
    # Canvas / drawing
    # ------------------------------------------------------------------
    @abstractmethod
    def draw(self) -> None:
        """
        Flush all pending changes to the screen (``canvas.draw_idle()``).
        """
        ...

    @abstractmethod
    def get_canvas(self) -> "FigureCanvasTkAgg":
        """Return the Tk canvas widget for this view."""
        ...

    @abstractmethod
    def get_figure(self) -> "Figure":
        """Return the matplotlib Figure for this view."""
        ...