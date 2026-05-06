"""
LineInspectionContextMixin — shared right-click context menu for the
line inspection panel (ax2).

:class:`ThreePanelView` inherits this mixin so that the line-inspection
context menu lives in one place and can be reused by any future view that
embeds an ax2-style line inspection panel.

Usage::

    class MyView(PlotView, LineInspectionContextMixin):
        def build_context_menu(self, event, canvas_widget):
            if event.inaxes == self.ax2:
                return self._build_line_inspection_menu(canvas_widget)
            return None
"""

from __future__ import annotations

from typing import Any, Optional


class LineInspectionContextMixin:
    """Mixin that supplies the line inspection panel right-click menu.

    Requires that the host class exposes ``self._islat`` (the main iSLAT
    application object with a ``GUI.top_bar`` attribute).
    """

    def _build_line_inspection_menu(self, canvas_widget: Any) -> Optional[Any]:
        """Return a populated ``tk.Menu`` for the line inspection panel.

        Parameters
        ----------
        canvas_widget :
            The Tk widget used as the menu's parent.
        """
        try:
            import tkinter as tk
        except ImportError:
            return None

        islat = getattr(self, '_islat', None)
        menu = tk.Menu(canvas_widget, tearoff=0)

        def _save_current_line():
            if islat is not None and hasattr(islat, 'GUI') and hasattr(islat.GUI, 'top_bar'):
                islat.GUI.top_bar.save_line(save_type="selected")

        def _fit_current_line():
            if islat is not None and hasattr(islat, 'GUI') and hasattr(islat.GUI, 'top_bar'):
                islat.GUI.top_bar.fit_selected_line(deblend=False)

        def _run_deblender():
            if islat is not None and hasattr(islat, 'GUI') and hasattr(islat.GUI, 'top_bar'):
                islat.GUI.top_bar.fit_selected_line(deblend=True)

        def _save_all_lines_in_range():
            if islat is not None and hasattr(islat, 'GUI') and hasattr(islat.GUI, 'top_bar'):
                islat.GUI.top_bar.save_all_lines_in_range()

        menu.add_command(label="Save Current Line",       command=_save_current_line)
        menu.add_command(label="Fit Current Line",        command=_fit_current_line)
        menu.add_command(label="Run Deblender",           command=_run_deblender)
        menu.add_separator()
        menu.add_command(label="Save All Lines in Range", command=_save_all_lines_in_range)
        return menu
