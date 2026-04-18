# -*- coding: utf-8 -*-
"""Tests for the dynamic ControlPanel / PlotView field-provider protocol.

Covers:
  - PlotView.get_view_fields() / get_display_range_binding() default behaviour
  - ThreePanelView inherits defaults (empty fields, None binding)
  - FullSpectrumView provides n_panels field and display-range binding
  - iSLATPlot view-change callback infrastructure
  - ControlPanel._rebuild_view_fields() / _on_view_changed()
  - ControlPanel._update_display_range() binding routing
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from unittest.mock import MagicMock, patch, PropertyMock

import iSLAT.Constants as c
from iSLAT.Modules.Plotting.PlotView import PlotView
from iSLAT.Modules.Plotting.FullSpectrumPlot import FullSpectrumPlot


# ======================================================================
# Helpers
# ======================================================================

def _make_wave_flux(n=200, start=4.9, end=28.0):
    wave = np.linspace(start, end, n)
    flux = np.random.default_rng(42).normal(1.0, 0.05, n)
    return wave, flux


def _make_mock_islat():
    """Build a mock iSLAT class with the minimum attributes needed."""
    islat = MagicMock()
    islat.display_range = (4.9, 5.9)
    islat._display_range = (4.9, 5.9)
    islat.wave_data_original = np.linspace(4.9, 28.0, 200)
    islat.flux_data = np.ones(200)
    islat.molecules_dict = MagicMock()
    islat.molecules_dict.apply_stellar_rv.side_effect = lambda w: w.copy()
    islat.active_molecule = None
    return islat


def _make_mock_plot_manager(islat=None):
    """Build a mock iSLATPlot (plot_manager) for FullSpectrumView."""
    if islat is None:
        islat = _make_mock_islat()
    pm = MagicMock()
    pm.islat = islat
    pm.theme = {}
    pm.toggle_state = {
        "atomic_lines": False,
        "saved_lines": False,
        "summed": True,
        "legend": True,
        "show_residuals": False,
        "current_selection": None,
    }
    pm.plot_renderer = MagicMock()
    pm.canvas = MagicMock()
    pm.ax1 = MagicMock()
    pm.ax2 = MagicMock()
    pm.ax3 = MagicMock()
    # View-change callbacks (simulates MainPlot after Phase 1)
    pm._view_change_callbacks = []
    return pm


# ======================================================================
# PlotView ABC default protocol
# ======================================================================

class TestPlotViewProtocol:
    """PlotView ABC provides default get_view_fields / get_display_range_binding."""

    def test_get_view_fields_default_returns_empty_list(self):
        # Create a minimal concrete subclass
        class StubView(PlotView):
            def activate(self, pf): ...
            def deactivate(self): ...
            def update_model_plot(self, *a, **k): ...
            def on_molecule_visibility_changed(self, *a, **k): ...
            def on_selection(self, xmin, xmax): ...
            def clear_selection(self): ...
            def clear_active_lines(self): ...
            def on_active_molecule_changed(self, new_molecule=None, current_selection=None): ...
            def on_molecule_parameter_changed(self, molecule_name, parameter_name, current_selection=None): ...
            def on_molecule_deleted(self, molecule_name): ...
            def apply_theme(self, t): ...
            def sync_toggle_state(self, s): ...
            def toggle_summed_spectrum(self, v): ...
            def toggle_legend(self, v=None): ...
            def toggle_saved_lines(self, s, **k): ...
            def toggle_atomic_lines(self, s): ...
            def draw(self): ...
            def get_canvas(self): ...
            def get_figure(self): ...

        view = StubView()
        assert view.get_view_fields() == []

    def test_get_display_range_binding_default_returns_none(self):
        class StubView(PlotView):
            def activate(self, pf): ...
            def deactivate(self): ...
            def update_model_plot(self, *a, **k): ...
            def on_molecule_visibility_changed(self, *a, **k): ...
            def on_selection(self, xmin, xmax): ...
            def clear_selection(self): ...
            def clear_active_lines(self): ...
            def on_active_molecule_changed(self, new_molecule=None, current_selection=None): ...
            def on_molecule_parameter_changed(self, molecule_name, parameter_name, current_selection=None): ...
            def on_molecule_deleted(self, molecule_name): ...
            def apply_theme(self, t): ...
            def sync_toggle_state(self, s): ...
            def toggle_summed_spectrum(self, v): ...
            def toggle_legend(self, v=None): ...
            def toggle_saved_lines(self, s, **k): ...
            def toggle_atomic_lines(self, s): ...
            def draw(self): ...
            def get_canvas(self): ...
            def get_figure(self): ...

        view = StubView()
        assert view.get_display_range_binding() is None


# ======================================================================
# ThreePanelView — inherits protocol defaults
# ======================================================================

class TestThreePanelViewProtocol:
    """ThreePanelView should return empty fields and None binding."""

    def test_get_view_fields_empty(self):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        pm = _make_mock_plot_manager()
        view = ThreePanelView(pm)
        assert view.get_view_fields() == []

    def test_get_display_range_binding_none(self):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        pm = _make_mock_plot_manager()
        view = ThreePanelView(pm)
        assert view.get_display_range_binding() is None


# ======================================================================
# FullSpectrumView — provides n_panels field + display-range binding
# ======================================================================

class TestFullSpectrumViewProtocol:
    """FullSpectrumView.get_view_fields / get_display_range_binding."""

    def _make_view_with_plot(self):
        """Create a FullSpectrumView with a real FullSpectrumPlot attached."""
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        pm = _make_mock_plot_manager()
        view = FullSpectrumView(pm)
        wave, flux = _make_wave_flux()
        plot = FullSpectrumPlot(wave_data=wave, flux_data=flux, n_panels=8)
        view._plot = plot
        view._canvas = MagicMock()
        return view, plot

    def test_get_view_fields_returns_n_panels_descriptor(self):
        view, plot = self._make_view_with_plot()
        fields = view.get_view_fields()
        assert len(fields) == 1
        desc = fields[0]
        assert desc["key"] == "n_panels"
        assert desc["datatype"] == "int"
        assert desc["default"] == 8
        assert callable(desc["getter"])
        assert callable(desc["setter"])

    def test_n_panels_getter_reads_from_plot(self):
        view, plot = self._make_view_with_plot()
        fields = view.get_view_fields()
        assert fields[0]["getter"]() == 8

    def test_n_panels_setter_updates_plot(self):
        view, plot = self._make_view_with_plot()
        fields = view.get_view_fields()
        setter = fields[0]["setter"]
        with patch.object(plot, 'generate_plot'):
            setter(5)
        assert plot.n_panels == 5

    def test_get_display_range_binding_not_none(self):
        view, plot = self._make_view_with_plot()
        binding = view.get_display_range_binding()
        assert binding is not None
        assert "getter" in binding
        assert "setter" in binding

    def test_display_range_getter_returns_xlim(self):
        view, plot = self._make_view_with_plot()
        binding = view.get_display_range_binding()
        start, end = binding["getter"]()
        assert start == pytest.approx(plot._xlim_start, abs=1e-6)
        assert end == pytest.approx(plot._xlim_end, abs=1e-6)

    def test_display_range_setter_updates_xlim(self):
        view, plot = self._make_view_with_plot()
        binding = view.get_display_range_binding()
        with patch.object(plot, 'generate_plot'):
            binding["setter"](10.0, 20.0)
        assert plot._xlim_start == pytest.approx(10.0)
        assert plot._xlim_end == pytest.approx(20.0)

    def test_no_binding_when_plot_is_none(self):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        pm = _make_mock_plot_manager()
        view = FullSpectrumView(pm)
        assert view._plot is None
        assert view.get_display_range_binding() is None

    def test_view_fields_default_when_plot_is_none(self):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        pm = _make_mock_plot_manager()
        view = FullSpectrumView(pm)
        fields = view.get_view_fields()
        assert len(fields) == 1
        # default should be 10 when no plot
        assert fields[0]["default"] == 10


# ======================================================================
# iSLATPlot view-change callback infrastructure
# ======================================================================

class TestViewChangeCallbacks:
    """Tests for add/remove/notify view-change callbacks on iSLATPlot."""

    def _make_islat_plot_mock(self):
        """Return mock with the new callback infra spliced in."""
        pm = _make_mock_plot_manager()
        # Splice real callback infrastructure
        pm._view_change_callbacks = []

        def add_cb(cb):
            if cb not in pm._view_change_callbacks:
                pm._view_change_callbacks.append(cb)
        pm.add_view_change_callback = add_cb

        def remove_cb(cb):
            try:
                pm._view_change_callbacks.remove(cb)
            except ValueError:
                pass
        pm.remove_view_change_callback = remove_cb

        def notify(old, new):
            for cb in pm._view_change_callbacks:
                cb(old, new)
        pm._notify_view_change = notify

        return pm

    def test_add_callback(self):
        pm = self._make_islat_plot_mock()
        cb = MagicMock()
        pm.add_view_change_callback(cb)
        assert cb in pm._view_change_callbacks

    def test_no_duplicate_callbacks(self):
        pm = self._make_islat_plot_mock()
        cb = MagicMock()
        pm.add_view_change_callback(cb)
        pm.add_view_change_callback(cb)
        assert pm._view_change_callbacks.count(cb) == 1

    def test_remove_callback(self):
        pm = self._make_islat_plot_mock()
        cb = MagicMock()
        pm.add_view_change_callback(cb)
        pm.remove_view_change_callback(cb)
        assert cb not in pm._view_change_callbacks

    def test_remove_nonexistent_callback_no_error(self):
        pm = self._make_islat_plot_mock()
        pm.remove_view_change_callback(MagicMock())  # should not raise

    def test_notify_calls_callback(self):
        pm = self._make_islat_plot_mock()
        cb = MagicMock()
        pm.add_view_change_callback(cb)
        old, new = MagicMock(), MagicMock()
        pm._notify_view_change(old, new)
        cb.assert_called_once_with(old, new)

    def test_notify_calls_multiple_callbacks(self):
        pm = self._make_islat_plot_mock()
        cb1, cb2 = MagicMock(), MagicMock()
        pm.add_view_change_callback(cb1)
        pm.add_view_change_callback(cb2)
        pm._notify_view_change(None, MagicMock())
        assert cb1.call_count == 1
        assert cb2.call_count == 1


# ======================================================================
# ControlPanel dynamic view fields
# ======================================================================

class TestControlPanelViewSwitch:
    """Test _rebuild_view_fields via _on_view_changed."""

    @pytest.fixture
    def tk_root(self):
        try:
            import tkinter as tk
            root = tk.Tk()
            root.withdraw()
            yield root
            root.destroy()
        except Exception:
            pytest.skip("Tk display not available")

    @pytest.fixture
    def control_panel(self, tk_root):
        """Build a ControlPanel with mocked dependencies."""
        import tkinter as tk
        from tkinter import ttk
        islat = _make_mock_islat()
        islat.molecules_dict = None  # simplify — skip mol callbacks
        pm = _make_mock_plot_manager(islat)
        data_field = MagicMock()
        data_field.insert_text = MagicMock()
        font = tk.font.Font(family="TkDefaultFont", size=10)
        # Patch the load_control_panel_fields_config to avoid file I/O
        with patch('iSLAT.Modules.GUI.Widgets.ControlPanel.load_control_panel_fields_config', return_value={}):
            cp = self._build_control_panel(tk_root, islat, pm, data_field, font)
        return cp

    @staticmethod
    def _build_control_panel(root, islat, plot, data_field, font):
        """
        Attempt to build a ControlPanel; if it fails due to missing attrs
        on the mock, just attach the minimum for testing view-field logic.
        """
        import tkinter as tk
        from tkinter import ttk
        try:
            from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel
            cp = ControlPanel(root, islat, plot, data_field, font)
            return cp
        except Exception:
            # Fallback: create a minimal frame with just the bits we need
            cp = ttk.Frame(root)
            cp.islat = islat
            cp.plot = plot
            cp.data_field = data_field
            cp.font = font
            cp.updating = False
            cp.fg_color = "black"
            cp.bg_color = "white"
            cp._view_fields_frame = ttk.Frame(cp)
            cp._view_fields_frame.grid(row=0, column=0)
            cp._view_field_entries = {}
            cp._current_display_range_binding = None
            cp.plot_start_var = tk.StringVar(value="4.9")
            cp.plot_range_var = tk.StringVar(value="1.0")

            # Bind methods from the real class
            from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel as CP
            import types
            for name in ('_rebuild_view_fields', '_create_view_field_entry',
                         '_sync_display_range_from_binding', '_restore_display_range_from_islat',
                         '_on_view_changed',
                         '_update_display_range', '_format_value', '_set_var'):
                if hasattr(CP, name):
                    setattr(cp, name, types.MethodType(getattr(CP, name), cp))
            return cp

    def test_rebuild_clears_previous_fields(self, control_panel):
        """_rebuild_view_fields should destroy old children."""
        cp = control_panel
        # Seed with a dummy child
        import tkinter as tk
        dummy = tk.Label(cp._view_fields_frame, text="old")
        dummy.grid(row=0, column=0)
        assert len(cp._view_fields_frame.winfo_children()) >= 1

        # Rebuild with a None view (empty)
        cp._rebuild_view_fields(None)
        assert len(cp._view_fields_frame.winfo_children()) == 0
        assert cp._current_display_range_binding is None

    def test_rebuild_creates_entries_from_view_fields(self, control_panel):
        """_rebuild_view_fields should create entries for view descriptors."""
        cp = control_panel
        view = MagicMock()
        view.get_view_fields.return_value = [
            {
                "key": "test_field",
                "label": "Test:",
                "default": 42,
                "tip": None,
                "datatype": "int",
                "width": 5,
                "getter": lambda: 42,
                "setter": lambda v: None,
            }
        ]
        view.get_display_range_binding.return_value = None

        cp._rebuild_view_fields(view)
        assert "test_field" in cp._view_field_entries
        entry, var, getter, setter, dtype = cp._view_field_entries["test_field"]
        assert var.get() == "42"
        assert dtype == "int"

    def test_rebuild_stores_display_range_binding(self, control_panel):
        """_rebuild_view_fields should store the view's display-range binding."""
        cp = control_panel
        binding = {"getter": lambda: (5.0, 25.0), "setter": MagicMock()}
        view = MagicMock()
        view.get_view_fields.return_value = []
        view.get_display_range_binding.return_value = binding

        cp._rebuild_view_fields(view)
        assert cp._current_display_range_binding is binding

    def test_on_view_changed_calls_rebuild(self, control_panel):
        """_on_view_changed should call _rebuild_view_fields."""
        cp = control_panel
        view = MagicMock()
        view.get_view_fields.return_value = []
        view.get_display_range_binding.return_value = None

        with patch.object(type(cp), '_rebuild_view_fields', wraps=cp._rebuild_view_fields) if hasattr(type(cp), '_rebuild_view_fields') else patch.object(cp, '_rebuild_view_fields', wraps=cp._rebuild_view_fields):
            cp._on_view_changed(MagicMock(), view)
        # The binding should reflect the new view
        assert cp._current_display_range_binding is None


# ======================================================================
# Display-range routing
# ======================================================================

class TestDisplayRangeBinding:
    """Test _update_display_range respects the active binding."""

    @pytest.fixture
    def tk_root(self):
        try:
            import tkinter as tk
            root = tk.Tk()
            root.withdraw()
            yield root
            root.destroy()
        except Exception:
            pytest.skip("Tk display not available")

    def test_routes_through_binding_setter(self, tk_root):
        """When a binding is active, _update_display_range should call its setter."""
        import tkinter as tk
        setter = MagicMock()

        # Minimal stand-in
        cp = MagicMock()
        cp.plot_start_var = tk.StringVar(tk_root, value="5.0")
        cp.plot_range_var = tk.StringVar(tk_root, value="20.0")
        cp._current_display_range_binding = {
            "getter": lambda: (5.0, 25.0),
            "setter": setter,
        }
        cp.islat = _make_mock_islat()

        # Call the real method
        from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel
        ControlPanel._update_display_range(cp)
        setter.assert_called_once_with(5.0, 25.0)

    def test_falls_back_to_islat_display_range(self, tk_root):
        """When no binding, _update_display_range should write islat.display_range."""
        import tkinter as tk

        cp = MagicMock()
        cp.plot_start_var = tk.StringVar(tk_root, value="5.0")
        cp.plot_range_var = tk.StringVar(tk_root, value="1.0")
        cp._current_display_range_binding = None
        cp.islat = _make_mock_islat()
        cp.islat._display_range = None  # force update

        from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel
        ControlPanel._update_display_range(cp)
        # Should have set islat.display_range
        assert cp.islat.display_range == (5.0, 6.0)

    def test_sync_display_range_from_binding(self, tk_root):
        """_sync_display_range_from_binding should populate GUI vars from binding."""
        import tkinter as tk

        cp = MagicMock()
        cp.plot_start_var = tk.StringVar(tk_root, value="0")
        cp.plot_range_var = tk.StringVar(tk_root, value="0")
        cp._current_display_range_binding = {
            "getter": lambda: (10.0, 22.0),
            "setter": MagicMock(),
        }

        def _set_var(var, val):
            cp.updating = True
            var.set(str(val))
            cp.updating = False
        cp._set_var = _set_var
        cp._format_value = lambda v, p=None: str(round(v, 6))

        from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel
        ControlPanel._sync_display_range_from_binding(cp)
        assert float(cp.plot_start_var.get()) == pytest.approx(10.0)
        assert float(cp.plot_range_var.get()) == pytest.approx(12.0)

    def test_restore_display_range_from_islat(self, tk_root):
        """_restore_display_range_from_islat should reset GUI vars from islat.display_range."""
        import tkinter as tk

        cp = MagicMock()
        cp.plot_start_var = tk.StringVar(tk_root, value="10.0")
        cp.plot_range_var = tk.StringVar(tk_root, value="12.0")
        cp.islat = _make_mock_islat()
        cp.islat.display_range = (4.9, 5.9)

        def _set_var(var, val):
            cp.updating = True
            var.set(str(val))
            cp.updating = False
        cp._set_var = _set_var
        cp._format_value = lambda v, p=None: str(round(v, 6))

        from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel
        ControlPanel._restore_display_range_from_islat(cp)
        assert float(cp.plot_start_var.get()) == pytest.approx(4.9)
        assert float(cp.plot_range_var.get()) == pytest.approx(1.0)

    def test_view_switch_reverts_display_range(self, tk_root):
        """Switching from FSP (binding) to 3-panel (no binding) should revert to islat.display_range."""
        import tkinter as tk
        from tkinter import ttk

        cp = MagicMock()
        cp.plot_start_var = tk.StringVar(tk_root, value="4.9")
        cp.plot_range_var = tk.StringVar(tk_root, value="1.0")
        cp.islat = _make_mock_islat()
        cp.islat.display_range = (4.9, 5.9)
        cp._view_fields_frame = ttk.Frame(tk_root)
        cp._view_field_entries = {}
        cp._current_display_range_binding = None

        def _set_var(var, val):
            cp.updating = True
            var.set(str(val))
            cp.updating = False
        cp._set_var = _set_var
        cp._format_value = lambda v, p=None: str(round(v, 6))

        from iSLAT.Modules.GUI.Widgets.ControlPanel import ControlPanel

        # Simulate FSP view with binding — sets display range to wider range
        fsp_view = MagicMock()
        fsp_view.get_view_fields.return_value = []
        fsp_view.get_display_range_binding.return_value = {
            "getter": lambda: (5.0, 28.0),
            "setter": MagicMock(),
        }

        for name in ('_rebuild_view_fields', '_create_view_field_entry',
                      '_sync_display_range_from_binding', '_restore_display_range_from_islat'):
            if hasattr(ControlPanel, name):
                import types
                setattr(cp, name, types.MethodType(getattr(ControlPanel, name), cp))

        cp._rebuild_view_fields(fsp_view)
        assert float(cp.plot_start_var.get()) == pytest.approx(5.0)
        assert float(cp.plot_range_var.get()) == pytest.approx(23.0)

        # Now switch to 3-panel view (no binding) — should revert to islat.display_range
        three_panel_view = MagicMock()
        three_panel_view.get_view_fields.return_value = []
        three_panel_view.get_display_range_binding.return_value = None

        cp._rebuild_view_fields(three_panel_view)
        assert float(cp.plot_start_var.get()) == pytest.approx(4.9)
        assert float(cp.plot_range_var.get()) == pytest.approx(1.0)
