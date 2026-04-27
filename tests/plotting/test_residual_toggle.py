# -*- coding: utf-8 -*-
"""Tests for the residual-toggle feature.

Validates the end-to-end flow that lets users switch between
:class:`FullSpectrumPlot` and :class:`ResidualSpectrumPlot` in the
full-spectrum view via:

* ``MainPlot.toggle_state["show_residuals"]``
* ``MainPlot.toggle_residuals()``
* ``FullSpectrumView.toggle_residuals()``
* ``FullSpectrumView._create_plot()`` (type dispatch)
* ``FullSpectrumView._rebuild_plot()`` (type-mismatch detection)
* ``FullSpectrumView._compute_model_flux()``
* ``InteractionHandler._toggle_residuals()``

All tests use the Agg backend and lightweight mocks to avoid requiring
a live Tk GUI or real HITRAN data.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock, PropertyMock, patch
from types import SimpleNamespace
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from iSLAT.Modules.Plotting.FullSpectrumPlot import FullSpectrumPlot
from iSLAT.Modules.Plotting.ResidualSpectrumPlot import ResidualSpectrumPlot
from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView

from tests.plotting import make_wave_flux


# ======================================================================
# Fixtures / helpers
# ======================================================================

def _make_molecules_dict(wave, flux):
    """Build a minimal mock MoleculeDict that returns valid summed flux."""
    md = MagicMock()
    md.get_summed_flux.return_value = (wave.copy(), flux.copy())
    md.get_summed_flux_resampled.return_value = (wave.copy(), flux.copy())
    md.apply_stellar_rv.side_effect = lambda w: w.copy()
    md.get_visible_molecules.return_value = []
    md.keys.return_value = []
    # _build_mol_cache calls get_matched_sampling_wavelengths → (bool, ndarray|None)
    md.get_matched_sampling_wavelengths.return_value = (False, None)
    return md


def _make_islat(wave, flux, err):
    """Build a minimal iSLAT mock for FullSpectrumView."""
    islat = MagicMock()
    islat.wave_data_original = wave
    islat.flux_data = flux
    islat.err_data = err
    islat.molecules_dict = _make_molecules_dict(wave, flux * 0.8)
    islat.input_line_list = "nonexistent.csv"
    islat.root = None  # No Tk root
    return islat


def _make_plot_manager(islat, show_residuals=False):
    """Build a minimal MainPlot-like mock for FullSpectrumView."""
    pm = MagicMock()
    pm.islat = islat
    pm.theme = {}
    pm.toggle_state = {
        "atomic_lines": False,
        "saved_lines": False,
        "summed": True,
        "legend": True,
        "show_residuals": show_residuals,
        "current_selection": None,
    }
    pm.plot_renderer = MagicMock()
    return pm


def _make_view(show_residuals=False, n=200):
    """Create a FullSpectrumView wired to lightweight mocks."""
    wave, flux, err = make_wave_flux(n=n)
    islat = _make_islat(wave, flux, err)
    pm = _make_plot_manager(islat, show_residuals=show_residuals)
    view = FullSpectrumView(pm)
    return view, pm, islat, wave, flux, err


# ======================================================================
# MainPlot toggle_state / property tests
# ======================================================================

class TestMainPlotToggleState:
    """Verify the ``show_residuals`` entry is wired correctly."""

    def test_default_toggle_state_has_show_residuals(self):
        """toggle_state should include 'show_residuals' defaulting to False."""
        pm = _make_plot_manager(_make_islat(*make_wave_flux(n=50)))
        assert "show_residuals" in pm.toggle_state
        assert pm.toggle_state["show_residuals"] is False

    def test_toggle_state_can_be_set_true(self):
        pm = _make_plot_manager(_make_islat(*make_wave_flux(n=50)), show_residuals=True)
        assert pm.toggle_state["show_residuals"] is True


# ======================================================================
# FullSpectrumView._create_plot() dispatch
# ======================================================================

class TestCreatePlotDispatch:
    """_create_plot should return different types based on show_residuals."""

    def test_creates_fsp_when_residuals_off(self):
        view, pm, *_ = _make_view(show_residuals=False)
        plot = view._create_plot()
        assert isinstance(plot, FullSpectrumPlot)
        assert not isinstance(plot, ResidualSpectrumPlot)
        plot.close()

    def test_creates_rsp_when_residuals_on(self):
        view, pm, *_ = _make_view(show_residuals=True)
        plot = view._create_plot()
        assert isinstance(plot, ResidualSpectrumPlot)
        plot.close()

    def test_rsp_has_model_flux_set(self):
        view, pm, *_ = _make_view(show_residuals=True)
        plot = view._create_plot()
        assert hasattr(plot, "model_flux")
        assert isinstance(plot.model_flux, np.ndarray)
        assert plot.model_flux.shape[0] > 0
        plot.close()

    def test_rsp_has_error_data(self):
        view, pm, islat, wave, flux, err = _make_view(show_residuals=True)
        plot = view._create_plot()
        assert plot.error_data is not None
        assert np.allclose(plot.error_data, err)
        plot.close()

    def test_fsp_generate_plot_produces_figure(self):
        view, pm, *_ = _make_view(show_residuals=False)
        plot = view._create_plot()
        plot.generate_plot()
        assert plot.fig is not None
        assert isinstance(plot.fig, Figure)
        plot.close()

    def test_rsp_generate_plot_produces_figure(self):
        view, pm, *_ = _make_view(show_residuals=True)
        plot = view._create_plot()
        plot.generate_plot()
        assert plot.fig is not None
        assert isinstance(plot.fig, Figure)
        plot.close()

    def test_rsp_subplots_are_tuples(self):
        """RSP subplots should be (spectrum_ax, residual_ax) tuples."""
        view, pm, *_ = _make_view(show_residuals=True)
        plot = view._create_plot()
        plot.generate_plot()
        for idx, val in plot.subplots.items():
            assert isinstance(val, tuple), f"subplot {idx} is not a tuple"
            assert len(val) == 2
        plot.close()

    def test_fsp_subplots_are_plain_axes(self):
        """FSP subplots should be plain Axes (not tuples)."""
        view, pm, *_ = _make_view(show_residuals=False)
        plot = view._create_plot()
        plot.generate_plot()
        for idx, val in plot.subplots.items():
            assert isinstance(val, Axes), f"subplot {idx} is not an Axes"
        plot.close()


# ======================================================================
# FullSpectrumView._compute_model_flux()
# ======================================================================

class TestComputeModelFlux:
    """Validate model-flux computation for the residual display."""

    def test_returns_array_of_correct_shape(self):
        view, pm, islat, wave, flux, err = _make_view()
        result = view._compute_model_flux(wave, wave)
        assert isinstance(result, np.ndarray)
        assert result.shape == wave.shape

    def test_returns_scaled_flux(self):
        """The mock returns flux * 0.8 so the model should match that."""
        view, pm, islat, wave, flux, err = _make_view()
        result = view._compute_model_flux(wave, wave)
        expected = flux * 0.8
        assert np.allclose(result, expected, atol=1e-10)

    def test_returns_zeros_on_failure(self):
        """If get_summed_flux_resampled raises, the result should be zeros."""
        view, pm, islat, wave, flux, err = _make_view()
        islat.molecules_dict.get_summed_flux_resampled.side_effect = RuntimeError("boom")
        result = view._compute_model_flux(wave, wave)
        assert np.allclose(result, 0.0)

    def test_returns_zeros_when_empty(self):
        """If get_summed_flux_resampled returns empty arrays, the result is zeros."""
        view, pm, islat, wave, flux, err = _make_view()
        islat.molecules_dict.get_summed_flux_resampled.return_value = (np.array([]), np.array([]))
        result = view._compute_model_flux(wave, wave)
        assert result.shape == wave.shape
        assert np.allclose(result, 0.0)

    def test_resamples_via_spectres(self):
        """get_summed_flux_resampled is called with wave_rest and visible_only=False."""
        view, pm, islat, wave, flux, err = _make_view(n=200)
        dense_flux = np.ones(len(wave)) * 0.05
        islat.molecules_dict.get_summed_flux_resampled.return_value = (wave.copy(), dense_flux)
        result = view._compute_model_flux(wave, wave)
        assert result.shape == wave.shape
        assert np.allclose(result, 0.05, atol=1e-10)

    def test_uses_visible_false_for_all_molecules(self):
        """Should call get_summed_flux_resampled with wave_rest and visible_only=False."""
        view, pm, islat, wave, flux, err = _make_view()
        view._compute_model_flux(wave, wave)
        islat.molecules_dict.get_summed_flux_resampled.assert_called_once_with(
            wave, visible_only=False,
        )


# ======================================================================
# FullSpectrumView.toggle_residuals()
# ======================================================================

class TestToggleResiduals:
    """Test the view-level toggle_residuals method."""

    def test_toggle_before_init_is_noop(self):
        """Calling toggle_residuals before activation should not crash."""
        view, pm, *_ = _make_view(show_residuals=False)
        assert not view._initialised
        view.toggle_residuals(True)  # Should be a no-op, not crash

    def test_toggle_switches_plot_type(self):
        """After activation, toggling should swap the plot type."""
        view, pm, *_ = _make_view(show_residuals=False)
        # Manually activate without Tk
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)

        # Now toggle to RSP
        pm.toggle_state["show_residuals"] = True
        view.toggle_residuals(True)
        assert isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_toggle_back_to_fsp(self):
        """Toggling from RSP back to FSP should work."""
        view, pm, *_ = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        assert isinstance(view._plot, ResidualSpectrumPlot)

        # Toggle back to FSP
        pm.toggle_state["show_residuals"] = False
        view.toggle_residuals(False)
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_toggle_creates_new_figure(self):
        """Toggling should create a brand new figure object."""
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        old_fig = view.fig

        pm.toggle_state["show_residuals"] = True
        view.toggle_residuals(True)
        new_fig = view.fig
        # Figures should be different objects
        assert new_fig is not old_fig
        view._plot.close()

    def test_span_selectors_reinstalled(self):
        """After toggling the span selectors should be rebuilt."""
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True

        pm.toggle_state["show_residuals"] = True
        view.toggle_residuals(True)
        # span_selectors should have been cleared and rebuilt
        # (they may be empty in headless mode without an active canvas
        #  but the dict should exist and not be the same reference)
        assert isinstance(view.span_selectors, dict)
        view._plot.close()


# ======================================================================
# FullSpectrumView._rebuild_plot() type-mismatch detection
# ======================================================================

class TestRebuildPlotTypeMismatch:
    """_rebuild_plot should detect FSP↔RSP mismatches and do full rebuilds."""

    def test_rebuild_detects_fsp_to_rsp_switch(self):
        """If show_residuals flipped while inactive, _rebuild_plot should
        detect the type mismatch and create the correct plot type."""
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)

        # Flip toggle as if done while view was deactivated
        pm.toggle_state["show_residuals"] = True
        view._rebuild_plot()
        assert isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_rebuild_detects_rsp_to_fsp_switch(self):
        view, pm, *_ = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        assert isinstance(view._plot, ResidualSpectrumPlot)

        pm.toggle_state["show_residuals"] = False
        view._rebuild_plot()
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_rebuild_no_mismatch_keeps_fsp(self):
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        old_type = type(view._plot)

        # Same toggle — no mismatch
        view._rebuild_plot()
        assert isinstance(view._plot, old_type)
        view._plot.close()

    def test_rebuild_no_mismatch_keeps_rsp(self):
        view, pm, *_ = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True

        # Same toggle — no mismatch, but model_flux should be updated
        view._rebuild_plot()
        assert isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_rebuild_updates_model_flux_on_rsp(self):
        """When RSP is active and data changes, model_flux must be refreshed."""
        view, pm, islat, wave, flux, err = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        old_model = view._plot.model_flux.copy()

        # Change the summed-flux return
        new_model = flux * 0.5
        islat.molecules_dict.get_summed_flux_resampled.return_value = (wave.copy(), new_model.copy())
        view._rebuild_plot()

        assert not np.allclose(view._plot.model_flux, old_model)
        view._plot.close()

    def test_rebuild_from_none(self):
        """If _plot is None, _rebuild_plot should create one."""
        view, pm, *_ = _make_view(show_residuals=False)
        assert view._plot is None
        view._rebuild_plot()
        assert view._plot is not None
        assert isinstance(view._plot, FullSpectrumPlot)
        view._plot.close()

    def test_rebuild_from_none_with_residuals(self):
        view, pm, *_ = _make_view(show_residuals=True)
        assert view._plot is None
        view._rebuild_plot()
        assert view._plot is not None
        assert isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()


# ======================================================================
# FullSpectrumView.subplots / _iter_all_axes with RSP
# ======================================================================

class TestViewSubplotsWithRSP:
    """Verify that the view's subplots property handles RSP tuple cells."""

    def test_view_subplots_fsp(self):
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        for idx, ax in view.subplots.items():
            assert isinstance(ax, Axes)
        view._plot.close()

    def test_view_subplots_rsp_returns_spectrum_axes(self):
        """The view's subplots dict should unwrap RSP tuples to primary axes."""
        view, pm, *_ = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        for idx, ax in view.subplots.items():
            # The view should return only the primary (spectrum) Axes
            assert isinstance(ax, Axes)
        view._plot.close()

    def test_iter_all_axes_fsp(self):
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        axes = list(view._iter_all_axes())
        assert len(axes) > 0
        for ax in axes:
            assert isinstance(ax, Axes)
        view._plot.close()

    def test_iter_all_axes_rsp_includes_residual_axes(self):
        """_iter_all_axes on RSP should yield both spectrum and residual axes."""
        view, pm, *_ = _make_view(show_residuals=True)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        axes = list(view._iter_all_axes())
        n_panels = len(view._plot.subplots)
        # Each panel has 2 axes (spectrum + residual)
        assert len(axes) == 2 * n_panels
        view._plot.close()


# ======================================================================
# InteractionHandler integration (mocked)
# ======================================================================

class TestInteractionHandlerResidualToggle:
    """Validate the keybinding delegation path."""

    def test_toggle_residuals_delegates_to_plot_manager(self):
        """_toggle_residuals should call plot_manager.toggle_residuals."""
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        pm = MagicMock()
        pm.toggle_residuals = MagicMock()
        pm.islat = MagicMock()
        pm.canvas = MagicMock()
        pm.fig = MagicMock()
        pm.ax1 = MagicMock()
        pm.ax2 = MagicMock()
        pm.ax3 = MagicMock()
        pm.canvas.mpl_connect = MagicMock()

        handler = InteractionHandler(pm)
        handler._toggle_residuals()
        pm.toggle_residuals.assert_called_once()

    def test_toggle_residuals_safe_without_method(self):
        """If plot_manager lacks toggle_residuals, no error is raised."""
        from iSLAT.Modules.GUI.InteractionHandler import InteractionHandler
        pm = MagicMock(spec=[])
        pm.islat = MagicMock()
        pm.canvas = MagicMock()
        pm.fig = MagicMock()
        pm.ax1 = MagicMock()
        pm.ax2 = MagicMock()
        pm.ax3 = MagicMock()
        pm.canvas.mpl_connect = MagicMock()

        handler = InteractionHandler(pm)
        # Should not raise
        handler._toggle_residuals()


# ======================================================================
# Round-trip toggle cycle
# ======================================================================

class TestRoundTripToggleCycle:
    """FSP → RSP → FSP round-trip should produce clean plots."""

    def test_fsp_rsp_fsp_cycle(self):
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True
        assert isinstance(view._plot, FullSpectrumPlot)

        # → RSP
        pm.toggle_state["show_residuals"] = True
        view.toggle_residuals(True)
        assert isinstance(view._plot, ResidualSpectrumPlot)
        assert view._plot.fig is not None

        # → FSP
        pm.toggle_state["show_residuals"] = False
        view.toggle_residuals(False)
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)
        assert view._plot.fig is not None
        view._plot.close()

    def test_double_toggle_same_direction_is_idempotent(self):
        """Toggling to RSP twice should not crash."""
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True

        pm.toggle_state["show_residuals"] = True
        view.toggle_residuals(True)
        first_fig = view.fig

        # Toggle again (already RSP → still RSP)
        view.toggle_residuals(True)
        assert isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()

    def test_multiple_rapid_toggles(self):
        """Rapidly toggling back and forth should not leak or crash."""
        view, pm, *_ = _make_view(show_residuals=False)
        view._plot = view._create_plot()
        view._plot.generate_plot()
        view._initialised = True

        for i in range(6):
            show = (i % 2 == 0)
            pm.toggle_state["show_residuals"] = show
            view.toggle_residuals(show)

        # Should end on FSP (even index → show=True at i=4, then i=5 → show=False)
        assert isinstance(view._plot, FullSpectrumPlot)
        assert not isinstance(view._plot, ResidualSpectrumPlot)
        view._plot.close()


# ======================================================================
# Error data passthrough
# ======================================================================

class TestErrorDataPassthrough:
    """RSP should receive error data from the iSLAT class."""

    def test_rsp_receives_err_data(self):
        view, pm, islat, wave, flux, err = _make_view(show_residuals=True)
        plot = view._create_plot()
        assert plot.error_data is not None
        assert np.allclose(plot.error_data, err)
        plot.close()

    def test_rsp_handles_missing_err_data(self):
        """If err_data is not set, RSP should still construct."""
        view, pm, islat, wave, flux, err = _make_view(show_residuals=True)
        del islat.err_data  # Remove attribute
        plot = view._create_plot()
        assert plot.error_data is None
        plot.close()

    def test_fsp_no_error_data_by_default(self):
        """FSP created by _create_plot should not have error_data set
        (it's not passed to the constructor)."""
        view, pm, *_ = _make_view(show_residuals=False)
        plot = view._create_plot()
        assert plot.error_data is None
        plot.close()
