# -*- coding: utf-8 -*-
"""Tests for StackedSpectralPanel panel-layout computation."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

from iSLAT.Modules.Plotting import FullSpectrumPlot

from tests.plotting import make_wave_flux


class TestPanelLayoutWavelength:
    """WAVELENGTH-mode panel-edge computation."""

    def _make_fsp(self, n_panels=5, **kwargs):
        wave, flux, _ = make_wave_flux(n=200)
        return FullSpectrumPlot(wave, flux, n_panels=n_panels, **kwargs)

    def test_n_panels_exact_count(self):
        """WAVELENGTH mode should produce exactly n_panels panels."""
        for n in (1, 3, 5, 7, 10, 13):
            fsp = self._make_fsp(n_panels=n)
            assert len(fsp._panel_edges) == n, (
                f"Expected {n} panels, got {len(fsp._panel_edges)}"
            )
            assert len(fsp._panel_ends) == n
            fsp.close()

    def test_panel_edges_cover_full_range(self):
        fsp = self._make_fsp(n_panels=5)
        assert fsp._panel_edges[0] == pytest.approx(fsp._xlim_start)
        assert fsp._panel_ends[-1] == pytest.approx(fsp._xlim_end)
        fsp.close()

    def test_panel_edges_contiguous(self):
        """Each panel end should equal the next panel start."""
        fsp = self._make_fsp(n_panels=7)
        for i in range(len(fsp._panel_edges) - 1):
            assert fsp._panel_ends[i] == pytest.approx(fsp._panel_edges[i + 1])
        fsp.close()

    def test_step_override(self):
        wave, flux, _ = make_wave_flux(wmin=10, wmax=20, n=200)
        fsp = FullSpectrumPlot(wave, flux, step=2.5)
        assert fsp._step == pytest.approx(2.5)
        assert len(fsp._panel_edges) == 4
        fsp.close()

    def test_xlim_range_restricts(self):
        wave, flux, _ = make_wave_flux(wmin=10, wmax=20, n=200)
        fsp = FullSpectrumPlot(
            wave, flux, n_panels=5, xlim_range=(12.0, 18.0),
        )
        assert fsp._xlim_start == pytest.approx(12.0)
        assert fsp._xlim_end == pytest.approx(18.0)
        fsp.close()


class TestPanelLayoutDataDensity:
    """DATA_DENSITY-mode panel-edge computation."""

    def _make_fsp(self, n_panels=5, **kwargs):
        wave, flux, _ = make_wave_flux(n=200)
        return FullSpectrumPlot(
            wave, flux, n_panels=n_panels, x_scaling="data_density", **kwargs,
        )

    def test_panel_count(self):
        fsp = self._make_fsp(n_panels=5)
        assert len(fsp._panel_edges) == 5
        assert len(fsp._panel_ends) == 5
        fsp.close()

    def test_boundaries_cover_range(self):
        fsp = self._make_fsp(n_panels=5)
        assert fsp._panel_edges[0] == pytest.approx(fsp._xlim_start)
        assert fsp._panel_ends[-1] == pytest.approx(fsp._xlim_end)
        fsp.close()

    def test_contiguous(self):
        fsp = self._make_fsp(n_panels=4)
        for i in range(len(fsp._panel_edges) - 1):
            assert fsp._panel_ends[i] == pytest.approx(fsp._panel_edges[i + 1])
        fsp.close()

    def test_step_is_none(self):
        fsp = self._make_fsp(n_panels=5)
        assert fsp._step is None
        fsp.close()


class TestActivePanelEdges:
    """Gap-skip panel filtering."""

    def _make_fsp(self, n_panels=5, **kwargs):
        wave, flux, _ = make_wave_flux(n=200)
        return FullSpectrumPlot(wave, flux, n_panels=n_panels, **kwargs)

    def test_all_active_when_connect(self):
        fsp = self._make_fsp(n_panels=5)
        idx, edges = fsp._active_panel_edges()
        assert len(idx) == 5

    def test_filters_empty_cells(self):
        wave, flux, _ = make_wave_flux(n=200, gap=(13.0, 17.0))
        fsp = FullSpectrumPlot(
            wave, flux, n_panels=5, gap_mode="skip",
        )
        idx, edges = fsp._active_panel_edges()
        assert len(idx) <= 5
        assert len(idx) >= 2

    def test_cell_has_data_true(self):
        fsp = self._make_fsp(n_panels=5)
        assert fsp._cell_has_data(10.0, 20.0) is True

    def test_cell_has_data_false(self):
        fsp = self._make_fsp(n_panels=5)
        assert fsp._cell_has_data(100.0, 200.0) is False
