# -*- coding: utf-8 -*-
"""Tests for the SpectralPanel ABC (via ConcreteSpectralPanel stub)."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.axes import Axes

from iSLAT.Modules.Plotting import GapMode, XScaling

from tests.plotting import (
    ConcreteSpectralPanel,
    make_wave_flux,
    make_atomic_lines,
    make_line_list,
)


class TestSpectralPanelProperties:
    def test_xlim_property(self):
        panel = ConcreteSpectralPanel(
            np.linspace(0, 10, 100), np.ones(100), 2.0, 8.0,
        )
        assert panel.xlim == (2.0, 8.0)

    def test_get_panel_mask(self):
        wave = np.linspace(0, 10, 101)
        panel = ConcreteSpectralPanel(wave, np.ones(101), 3.0, 7.0)
        mask = panel.get_panel_mask()
        assert mask.sum() > 0
        assert np.all(wave[mask] >= 3.0)
        assert np.all(wave[mask] <= 7.0)

    def test_get_panel_data(self):
        wave = np.linspace(0, 10, 101)
        flux = np.ones(101) * 0.5
        err = np.ones(101) * 0.01
        panel = ConcreteSpectralPanel(wave, flux, 2.0, 5.0, error_data=err)
        pw, pf, pe = panel.get_panel_data()
        assert len(pw) == len(pf) == len(pe)
        assert np.all(pw >= 2.0)
        assert np.all(pw <= 5.0)

    def test_get_panel_data_no_error(self):
        wave = np.linspace(0, 10, 50)
        panel = ConcreteSpectralPanel(wave, np.ones(50), 0.0, 10.0)
        pw, pf, pe = panel.get_panel_data()
        assert pe is None

    def test_compute_ylim(self):
        wave = np.linspace(0, 10, 100)
        flux = np.random.default_rng(0).uniform(0.1, 0.5, 100)
        panel = ConcreteSpectralPanel(wave, flux, 0.0, 10.0)
        ymin, ymax = panel.compute_ylim()
        assert ymin < 0
        assert ymax > np.max(flux)

    def test_panel_has_data_true(self):
        wave = np.linspace(0, 10, 100)
        panel = ConcreteSpectralPanel(wave, np.ones(100), 0.0, 10.0)
        assert panel.panel_has_data() is True

    def test_panel_has_data_false_nan(self):
        wave = np.linspace(0, 10, 100)
        flux = np.full(100, np.nan)
        panel = ConcreteSpectralPanel(wave, flux, 0.0, 10.0)
        assert panel.panel_has_data() is False

    def test_panel_has_data_false_empty_range(self):
        wave = np.linspace(0, 10, 100)
        panel = ConcreteSpectralPanel(wave, np.ones(100), 20.0, 30.0)
        assert panel.panel_has_data() is False

    def test_gap_mode_default(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        assert panel.gap_mode is GapMode.CONNECT

    def test_gap_mode_from_string(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            gap_mode="skip",
        )
        assert panel.gap_mode is GapMode.SKIP

    def test_x_scaling_default(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        assert panel.x_scaling is XScaling.WAVELENGTH

    def test_x_scaling_from_string(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            x_scaling="data_density",
        )
        assert panel.x_scaling is XScaling.DATA_DENSITY

    def test_set_range(self):
        wave = np.linspace(0, 20, 200)
        panel = ConcreteSpectralPanel(wave, np.ones(200), 0.0, 10.0)
        panel.set_range(5.0, 15.0)
        assert panel.xmin == 5.0
        assert panel.xmax == 15.0

    def test_resolve_axes_external(self):
        fig, ax = plt.subplots()
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0, ax=ax,
        )
        resolved = panel._resolve_axes()
        assert resolved is ax
        plt.close(fig)

    def test_resolve_axes_creates_own(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        ax = panel._resolve_axes()
        assert ax is not None
        assert isinstance(ax, Axes)
        panel.close()


class TestGapDetection:
    def test_detect_gaps_nan(self):
        wave = np.linspace(10, 20, 100)
        flux = np.ones(100)
        flux[30:50] = np.nan
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        gaps = panel.detect_gaps()
        assert len(gaps) >= 1
        g_start, g_end = gaps[0]
        assert g_start < wave[40]
        assert g_end > wave[40]

    def test_detect_gaps_large_jump(self):
        """A wavelength jump > 5 % of the span should be detected."""
        wave = np.concatenate([
            np.linspace(10, 12, 50),
            np.linspace(18, 20, 50),
        ])
        flux = np.ones(100)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        gaps = panel.detect_gaps()
        assert len(gaps) >= 1
        g_start, g_end = gaps[0]
        assert g_start <= 12.0
        assert g_end >= 18.0

    def test_detect_gaps_no_gap(self):
        wave = np.linspace(10, 20, 200)
        flux = np.ones(200)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        gaps = panel.detect_gaps()
        assert len(gaps) == 0

    def test_detect_gaps_short_data(self):
        panel = ConcreteSpectralPanel(
            np.array([5.0]), np.array([1.0]), 5.0, 5.0,
        )
        assert panel.detect_gaps() == []

    def test_detect_gaps_custom_threshold(self):
        wave = np.concatenate([
            np.linspace(10, 14, 40),
            np.linspace(16, 20, 40),
        ])
        flux = np.ones(80)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        gaps_default = panel.detect_gaps()
        assert len(gaps_default) >= 1


class TestPanelDataWithGaps:
    def test_connect_mode_unchanged(self):
        """CONNECT mode should return the same data as get_panel_data."""
        wave = np.linspace(10, 20, 100)
        flux = np.ones(100)
        flux[40:60] = np.nan
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.CONNECT,
        )
        w1, f1, _ = panel.get_panel_data()
        w2, f2, _ = panel.get_panel_data_with_gaps()
        np.testing.assert_array_equal(w1, w2)
        np.testing.assert_array_equal(f1, f2)

    def test_skip_mode_inserts_nan(self):
        wave = np.linspace(10, 20, 100)
        flux = np.ones(100)
        flux[40:60] = np.nan
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        w_orig, _, _ = panel.get_panel_data()
        w_gap, f_gap, _ = panel.get_panel_data_with_gaps()
        assert len(w_gap) >= len(w_orig)
        nan_count = np.sum(np.isnan(f_gap))
        assert nan_count > np.sum(np.isnan(flux[(wave >= 10) & (wave <= 20)]))


class TestGapIndicators:
    def test_no_gaps(self):
        wave = np.linspace(10, 20, 200)
        flux = np.ones(200)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        gap_texts = [t for t in panel.ax.texts if hasattr(t, panel._GAP_TAG)]
        assert len(gap_texts) == 0
        panel.close()

    def test_internal_gap(self):
        """Internal gaps get white-out, break marks, and annotation."""
        wave = np.concatenate([
            np.linspace(10, 13, 100),
            np.linspace(17, 20, 100),
        ])
        flux = np.ones(200)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        gap_texts = [t for t in panel.ax.texts if hasattr(t, panel._GAP_TAG)]
        assert len(gap_texts) >= 1
        panel.close()

    def test_edge_gap_tightens_xlim(self):
        """Edge gaps should tighten the axes xlim.

        An 'edge gap' is one whose start or end coincides with the
        panel boundary.  Create a small cluster of NaN points at the
        left edge (10–12), then real data from 12–20.  The NaN region
        creates a detected gap that touches the left boundary,
        triggering tightening.
        """
        wave = np.linspace(10, 20, 200)
        flux = np.ones(200)
        # NaN-out the left 20% to create a left-edge gap
        flux[wave < 12.0] = np.nan
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        vis_lo, vis_hi = panel.ax.get_xlim()
        # The left edge should have been tightened past 10.0
        assert vis_lo > 10.5, (
            f"Expected left edge tightened past 10.5, got {vis_lo}"
        )
        panel.close()

    def test_remove_gap_indicators(self):
        wave = np.concatenate([
            np.linspace(10, 13, 100),
            np.linspace(17, 20, 100),
        ])
        flux = np.ones(200)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        gap_texts = [t for t in panel.ax.texts if hasattr(t, panel._GAP_TAG)]
        assert len(gap_texts) >= 1
        panel.remove_gap_indicators()
        gap_texts_after = [t for t in panel.ax.texts if hasattr(t, panel._GAP_TAG)]
        assert len(gap_texts_after) == 0
        panel.close()

    def test_explicit_gaps(self):
        wave = np.linspace(10, 20, 200)
        flux = np.ones(200)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        explicit_gaps = [(13.0, 17.0)]
        panel.draw_gap_indicators(gaps=explicit_gaps)
        gap_texts = [t for t in panel.ax.texts if hasattr(t, panel._GAP_TAG)]
        assert len(gap_texts) >= 1
        panel.close()


class TestAnnotationsAfterGapTightening:
    def test_plot_atomic_lines_uses_visible_xlim(self):
        """After gap tightening, annotations should use the visible xlim."""
        wave = np.concatenate([
            np.linspace(10, 10.5, 10),
            np.linspace(17, 20, 100),
        ])
        flux = np.ones(110)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        atomic_df = make_atomic_lines([12.0, 18.0])
        panel.plot_atomic_lines(atomic_df)
        tagged = [t for t in panel.ax.texts if hasattr(t, "_islat_atomic_line")]
        wavelengths_annotated = [t.get_position()[0] for t in tagged]
        assert any(17.5 < w < 19.0 for w in wavelengths_annotated)
        panel.close()

    def test_plot_saved_lines_uses_visible_xlim(self):
        wave = np.concatenate([
            np.linspace(10, 10.5, 10),
            np.linspace(17, 20, 100),
        ])
        flux = np.ones(110)
        panel = ConcreteSpectralPanel(
            wave, flux, 10.0, 20.0, gap_mode=GapMode.SKIP,
        )
        panel.generate_plot()
        panel.draw_gap_indicators()
        line_df = make_line_list([12.0, 18.0])
        panel.plot_saved_lines(line_df)
        tagged = [t for t in panel.ax.texts if hasattr(t, "_islat_saved_line")]
        assert len(tagged) >= 1
        panel.close()
