# -*- coding: utf-8 -*-
"""Tests for GUI widget and utility classes.

Tkinter widgets require a Tk root. Every test that needs one uses the
``tk_root`` fixture which creates and destroys a Tk instance; if the
display server is unavailable (headless CI) the test is skipped.
"""

import platform
import pytest

# ---------------------------------------------------------------------------
# Tk root fixture -- skip gracefully on headless machines
# ---------------------------------------------------------------------------

@pytest.fixture
def tk_root():
    """Create a Tk root for widget tests; skip if no display."""
    try:
        import tkinter as tk
        root = tk.Tk()
        root.withdraw()          # keep window hidden
        yield root
        root.destroy()
    except Exception:
        pytest.skip("Tk display not available (headless environment)")


# ========================================================================
# Tooltip module tests
# ========================================================================

class TestTooltipTheme:
    """Tests for the module-level tooltip theme helpers."""

    def test_set_tooltip_theme_stores_theme(self):
        import sys
        from iSLAT.Modules.GUI.Tooltips import set_tooltip_theme
        theme = {"tooltip_background": "#222", "tooltip_foreground": "#eee"}
        set_tooltip_theme(theme)
        # Access updated module variable via sys.modules
        T = sys.modules['iSLAT.Modules.GUI.Tooltips']
        assert T._current_theme is theme
        # Cleanup
        set_tooltip_theme(None)

    def test_resolve_tooltip_colors_uses_theme(self):
        from iSLAT.Modules.GUI.Tooltips import (
            set_tooltip_theme,
            _resolve_tooltip_colors,
        )
        theme = {"tooltip_background": "#111", "tooltip_foreground": "#fff"}
        set_tooltip_theme(theme)
        bg, fg = _resolve_tooltip_colors()
        assert bg == "#111"
        assert fg == "#fff"
        set_tooltip_theme(None)

    def test_resolve_tooltip_colors_explicit_overrides(self):
        from iSLAT.Modules.GUI.Tooltips import (
            set_tooltip_theme,
            _resolve_tooltip_colors,
        )
        theme = {"tooltip_background": "#111", "tooltip_foreground": "#fff"}
        set_tooltip_theme(theme)
        bg, fg = _resolve_tooltip_colors(explicit_bg="red", explicit_fg="blue")
        assert bg == "red"
        assert fg == "blue"
        set_tooltip_theme(None)

    def test_resolve_tooltip_colors_no_theme_uses_defaults(self):
        from iSLAT.Modules.GUI.Tooltips import (
            set_tooltip_theme,
            _resolve_tooltip_colors,
            _LIGHT_DEFAULTS,
        )
        set_tooltip_theme(None)
        bg, fg = _resolve_tooltip_colors()
        assert bg == _LIGHT_DEFAULTS["tooltip_background"]
        assert fg == _LIGHT_DEFAULTS["tooltip_foreground"]


class TestToolTipWidget:
    """Tests for ToolTip / CreateToolTip that need a Tk root."""

    def test_create_tooltip(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.Tooltips import CreateToolTip

        btn = tk.Button(tk_root, text="hover me")
        btn.pack()
        tt = CreateToolTip(btn, "helpful text")
        assert tt.text == "helpful text"
        assert tt.active is True

    def test_tooltip_disable_enable(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.Tooltips import CreateToolTip

        btn = tk.Button(tk_root, text="test")
        btn.pack()
        tt = CreateToolTip(btn, "tip")
        tt.disable()
        assert tt.active is False
        tt.enable()
        assert tt.active is True

    def test_tooltip_showtip_hidetip(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.Tooltips import CreateToolTip

        btn = tk.Button(tk_root, text="tip")
        btn.pack()
        tk_root.update_idletasks()
        tt = CreateToolTip(btn, "some tip")
        # Manually show / hide
        tt.showtip()
        assert tt.tipwindow is not None
        tt.hidetip()
        assert tt.tipwindow is None


# ========================================================================
# GUIFunctions tests
# ========================================================================

class TestGUIFunctions:
    """Tests for standalone GUI function utilities."""

    def test_configure_all_button_styles(self, tk_root):
        from iSLAT.Modules.GUI.GUIFunctions import configure_all_button_styles
        theme = {
            "foreground": "#000",
            "buttons": {"DefaultBotton": {"background": "lightgray", "active_background": "gray"}},
            "delete_button_bg_color": "#b1403b",
            "delete_button_fg_color": "#000000",
        }
        # Should not raise
        configure_all_button_styles(theme)

    def test_create_button(self, tk_root):
        import tkinter as tk
        from tkinter import ttk
        from iSLAT.Modules.GUI.GUIFunctions import create_button

        frame = ttk.Frame(tk_root)
        frame.pack()
        theme = {
            "foreground": "#000",
            "buttons": {"DefaultBotton": {"background": "#ccc", "active_background": "#aaa"}},
            "delete_button_bg_color": "#b1403b",
            "delete_button_fg_color": "#000",
        }
        btn = create_button(frame, theme, "Test", lambda: None, 0, 0, tip_text="A tip")
        assert isinstance(btn, ttk.Button)

    def test_create_menu_btn(self, tk_root):
        import tkinter as tk
        from tkinter import ttk
        from iSLAT.Modules.GUI.GUIFunctions import create_menu_btn

        frame = ttk.Frame(tk_root)
        frame.pack()
        theme = {
            "foreground": "#000",
            "buttons": {"DefaultBotton": {"background": "#ccc", "active_background": "#aaa"}},
            "delete_button_bg_color": "#b1403b",
            "delete_button_fg_color": "#000",
        }
        mb = create_menu_btn(frame, theme, "Menu", 0, 0)
        assert isinstance(mb, ttk.Menubutton)

    def test_create_wrapper_frame(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.GUIFunctions import create_wrapper_frame

        parent = tk.Frame(tk_root)
        parent.pack()
        wrapper = create_wrapper_frame(parent, row=0, col=0, bg="darkgrey")
        assert isinstance(wrapper, tk.Frame)

    def test_create_scrollable_frame(self, tk_root):
        import tkinter as tk
        from tkinter import ttk
        from iSLAT.Modules.GUI.GUIFunctions import create_scrollable_frame

        parent = ttk.Frame(tk_root)
        parent.pack()
        scroll_frame = create_scrollable_frame(
            parent, height=100, width=200, vertical=True, horizontal=True,
        )
        assert isinstance(scroll_frame, ttk.Frame)


# ========================================================================
# ColorButton tests
# ========================================================================

class TestColorButton:
    """Tests for the ColorButton widget."""

    def test_color_button_init(self, tk_root):
        from iSLAT.Modules.GUI.GUIFunctions import ColorButton
        cb = ColorButton(tk_root, color="#4CAF50")
        cb.pack()
        tk_root.update_idletasks()
        assert cb.color == "#4caf50"

    def test_color_button_darken(self, tk_root):
        from iSLAT.Modules.GUI.GUIFunctions import ColorButton
        cb = ColorButton(tk_root, color="#ffffff")
        cb.pack()
        tk_root.update_idletasks()
        dark = cb._darken_color("#ffffff", 0.5)
        # 255 * 0.5 = 127 -> #7f7f7f
        assert dark == "#7f7f7f"

    def test_color_button_change_color(self, tk_root):
        from iSLAT.Modules.GUI.GUIFunctions import ColorButton
        cb = ColorButton(tk_root, color="#ff0000")
        cb.pack()
        tk_root.update_idletasks()
        cb.change_color("#00ff00")
        assert cb.color == "#00ff00"

    def test_color_button_command(self, tk_root):
        from iSLAT.Modules.GUI.GUIFunctions import ColorButton
        called = []
        cb = ColorButton(tk_root, color="#aabbcc", command=lambda: called.append(1))
        cb.pack()
        assert cb.command is not None
        cb.add_command(lambda: called.append(2))
        cb.command()
        assert called == [2]


# ========================================================================
# DataField tests
# ========================================================================

class TestDataField:
    """Tests for the DataField text output widget."""

    def test_datafield_init(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.DataField import DataField
        df = DataField(value="test", master=tk_root)
        df.pack()
        assert df.value == "test"

    def test_datafield_insert_and_clear(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.DataField import DataField
        df = DataField(value="", master=tk_root)
        df.pack()
        tk_root.update_idletasks()
        df.insert_text("hello world", clear_after=False)
        content = df.text.get("1.0", "end").strip()
        assert "hello world" in content

        df.clear()
        content = df.text.get("1.0", "end").strip()
        assert content == ""

    def test_datafield_clear_before_next(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.DataField import DataField
        df = DataField(value="", master=tk_root)
        df.pack()
        tk_root.update_idletasks()
        df.insert_text("first", clear_after=True)
        assert df.clear_before_next is True
        df.insert_text("second", clear_after=False)
        content = df.text.get("1.0", "end").strip()
        # "first" should have been cleared before inserting "second"
        assert "first" not in content
        assert "second" in content

    def test_datafield_repr(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.DataField import DataField
        df = DataField(value=42, master=tk_root)
        df.pack()
        # DataField.__repr__ references self.name which is not set
        # in __init__, so it raises AttributeError (source code bug)
        with pytest.raises(AttributeError):
            repr(df)


# ========================================================================
# trim_label tests
# ========================================================================

class TestTrimLabel:
    """Tests for the trim_label widget in FileInteractionPane."""

    def test_short_text_not_trimmed(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.FileInteractionPane import trim_label
        lbl = trim_label(tk_root, max_len=25, text="short")
        lbl.pack()
        assert lbl.trimmed is False
        assert lbl.cget("text") == "short"

    def test_long_text_trimmed(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.FileInteractionPane import trim_label
        long_text = "A" * 50
        lbl = trim_label(tk_root, max_len=25, text=long_text)
        lbl.pack()
        assert lbl.trimmed is True
        displayed = lbl.cget("text")
        assert displayed.endswith("...")
        assert len(displayed) == 25

    def test_exact_length_not_trimmed(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.FileInteractionPane import trim_label
        text_25 = "A" * 25
        lbl = trim_label(tk_root, max_len=25, text=text_25)
        lbl.pack()
        assert lbl.trimmed is False


# ========================================================================
# ResizableFrame tests
# ========================================================================

class TestResizableFrame:
    """Tests for the ResizableFrame container widget."""

    def test_init_defaults(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.ResizableFrame import ResizableFrame
        rf = ResizableFrame(tk_root)
        rf.pack()
        assert rf.orientation == "vertical"
        assert rf.sash_size == 4
        assert rf.frames == []
        assert rf.sashes == []

    def test_add_frame(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.Widgets.ResizableFrame import ResizableFrame
        rf = ResizableFrame(tk_root)
        rf.pack()
        child1 = tk.Frame(rf)
        child2 = tk.Frame(rf)
        rf.add_frame(child1, weight=1)
        rf.add_frame(child2, weight=2)
        assert len(rf.frames) == 2
        assert len(rf.sashes) == 1  # one sash between two frames
        assert rf.total_weight == 3

    def test_horizontal_orientation(self, tk_root):
        from iSLAT.Modules.GUI.Widgets.ResizableFrame import ResizableFrame
        rf = ResizableFrame(tk_root, orientation="horizontal")
        rf.pack()
        assert rf.orientation == "horizontal"

    def test_theme_applied_to_sashes(self, tk_root):
        import tkinter as tk
        from iSLAT.Modules.GUI.Widgets.ResizableFrame import ResizableFrame
        rf = ResizableFrame(tk_root)
        rf.pack()
        rf.add_frame(tk.Frame(rf))
        rf.add_frame(tk.Frame(rf))
        theme = {"toolbar": "#555555", "selection_color": "#777777"}
        rf._apply_theme_to_sashes(theme)
        # Sash colour should match theme
        assert rf.sashes[0].cget("bg") == "#555555"


# ========================================================================
# DisplayConfig tests (non-Tk, pure Python)
# ========================================================================

class TestDisplayConfigExtended:
    """Extended tests for DisplayConfig beyond what test_display_config.py covers."""

    def test_display_config_singleton_frozen(self):
        from iSLAT.Modules.GUI.DisplayConfig import display_config
        with pytest.raises(AttributeError):
            display_config.scale_factor = 999.0

    def test_display_config_attributes(self):
        from iSLAT.Modules.GUI.DisplayConfig import display_config
        assert display_config.figure_dpi == 100
        assert display_config.savefig_dpi == 300
        assert display_config.text_antialiased is True
        assert display_config.lines_antialiased is True

    def test_apply_matplotlib_defaults_idempotent(self):
        from iSLAT.Modules.GUI.DisplayConfig import apply_matplotlib_defaults
        import matplotlib
        apply_matplotlib_defaults()
        dpi1 = matplotlib.rcParams["figure.dpi"]
        apply_matplotlib_defaults()
        dpi2 = matplotlib.rcParams["figure.dpi"]
        assert dpi1 == dpi2

    def test_detect_scale_factor_positive(self):
        from iSLAT.Modules.GUI.DisplayConfig import detect_scale_factor
        sf = detect_scale_factor()
        assert sf > 0


# ========================================================================
# InteractionHandler tests (with mocks)
# ========================================================================

class TestInteractionHandler:
    """Tests for InteractionHandler using a mocked plot_manager."""

    @pytest.fixture
    def mock_plot_manager(self):
        """Build a minimal mock that satisfies InteractionHandler.__init__."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from unittest.mock import MagicMock

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        canvas = fig.canvas

        pm = MagicMock()
        pm.fig = fig
        pm.ax1 = ax1
        pm.ax2 = ax2
        pm.ax3 = ax3
        pm.canvas = canvas
        pm.islat = MagicMock()
        pm.parent_frame = None  # no Tk frame
        yield pm
        plt.close(fig)

    def test_init_creates_span_selector(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        assert ih.span_selector is not None
        assert ih.selected_range is None

    def test_init_callback_dicts(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        assert isinstance(ih.selection_callbacks, dict)
        assert isinstance(ih.click_callbacks, dict)
        assert isinstance(ih.zoom_callbacks, dict)

    def test_double_click_threshold(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        assert ih.double_click_threshold == 0.5

    def test_on_span_select_zero_range_clears(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        # Simulate selecting a zero-width range
        ih._on_span_select(15.0, 15.0)
        # Should have cleared the selection (via clear_current_selection)
        assert ih.selected_range is None

    def test_on_span_select_stores_range(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        ih._on_span_select(12.0, 18.0)
        assert ih.selected_range == (12.0, 18.0)

    def test_on_span_select_swaps_reversed_range(self, mock_plot_manager):
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        ih = InteractionHandler(mock_plot_manager)
        ih._on_span_select(18.0, 12.0)
        assert ih.selected_range == (12.0, 18.0)
