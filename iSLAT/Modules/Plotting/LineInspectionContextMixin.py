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

        # ------------------------------------------------------------------
        # Spectrum-navigation helpers — only available from ThreePanelView
        # (i.e. when self has ax1 / ax2 / _pm).
        # ------------------------------------------------------------------
        ax1 = getattr(self, 'ax1', None)
        ax2 = getattr(self, 'ax2', None)
        pm  = getattr(self, '_pm', None)

        def _get_ax2_range():
            """Return (xmin, xmax) from ax2, or None if unavailable."""
            if ax2 is None:
                return None
            try:
                x0, x1 = ax2.get_xlim()
                if x0 < x1:
                    return x0, x1
            except Exception:
                pass
            return None

        def _refresh_ax1():
            """Apply new ax1 limits and refresh via the plot manager."""
            if pm is not None and hasattr(pm, 'match_display_range'):
                pm.match_display_range(match_y=True)
            elif pm is not None and hasattr(pm, 'canvas'):
                pm.canvas.draw_idle()

        def _match_spectrum_to_inspection():
            """Set ax1 xlim to the current inspection (ax2) range."""
            r = _get_ax2_range()
            if r is None or ax1 is None:
                return
            xmin, xmax = r
            ax1.set_xlim(xmin, xmax)
            if pm is not None and hasattr(pm, 'islat'):
                pm.islat.display_range = (xmin, xmax)
            _refresh_ax1()

        def _center_spectrum_on_inspection():
            """Center the spectrum view (ax1) on the inspection midpoint,
            keeping the current ax1 x-span unchanged."""
            r = _get_ax2_range()
            if r is None or ax1 is None:
                return
            insp_mid = (r[0] + r[1]) / 2.0
            try:
                cur_lo, cur_hi = ax1.get_xlim()
                half_span = (cur_hi - cur_lo) / 2.0
            except Exception:
                return
            new_lo = insp_mid - half_span
            new_hi = insp_mid + half_span
            ax1.set_xlim(new_lo, new_hi)
            if pm is not None and hasattr(pm, 'islat'):
                pm.islat.display_range = (new_lo, new_hi)
            _refresh_ax1()

        menu.add_command(label="Save Current Line",       command=_save_current_line)
        menu.add_command(label="Fit Current Line",        command=_fit_current_line)
        menu.add_command(label="Run Deblender",           command=_run_deblender)
        menu.add_separator()
        menu.add_command(label="Save All Lines in Range", command=_save_all_lines_in_range)

        if ax1 is not None and ax2 is not None:
            menu.add_separator()
            menu.add_command(
                label="Match Spectrum to Inspection Range",
                command=_match_spectrum_to_inspection,
            )
            menu.add_command(
                label="Center Spectrum on Inspection Location",
                command=_center_spectrum_on_inspection,
            )

        return menu