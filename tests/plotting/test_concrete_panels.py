# -*- coding: utf-8 -*-
"""Tests for the concrete SpectrumPanel and ResidualPanel classes."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from iSLAT.Modules.Plotting import GapMode, SpectrumPanel, ResidualPanel

from tests.plotting import make_wave_flux


class TestSpectrumPanel:
    def test_basic_render(self):
        wave, flux, err = make_wave_flux(n=100)
        panel = SpectrumPanel(
            wave, flux, 10.0, 20.0,
            error_data=err,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert panel.ax is not None
        assert len(panel.ax.lines) > 0
        plt.close("all")

    def test_with_mol_cache(self):
        wave, flux, err = make_wave_flux(n=100)
        mol_wave = np.linspace(10, 20, 100)
        mol_flux = 0.02 * np.exp(-((mol_wave - 15) ** 2) / 2)
        cache = [(mol_wave, mol_flux, "blue", "H2O", "h2o")]
        panel = SpectrumPanel(
            wave, flux, 10.0, 20.0,
            mol_cache=cache,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        mol_lines = [l for l in panel.ax.lines if hasattr(l, "_molecule_name")]
        assert len(mol_lines) == 1
        assert mol_lines[0]._molecule_name == "h2o"
        plt.close("all")

    def test_with_summed_spectrum(self):
        wave, flux, _ = make_wave_flux(n=100)
        summed_w = np.linspace(10, 20, 100)
        summed_f = 0.05 * np.ones(100)
        panel = SpectrumPanel(
            wave, flux, 10.0, 20.0,
            summed_wave=summed_w,
            summed_flux=summed_f,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        summed_coll = [c for c in panel.ax.collections if hasattr(c, "_islat_summed")]
        assert len(summed_coll) == 1
        plt.close("all")

    def test_gap_mode_skip(self):
        wave, flux, _ = make_wave_flux(n=200, gap=(14.0, 16.0))
        panel = SpectrumPanel(
            wave, flux, 10.0, 20.0,
            gap_mode=GapMode.SKIP,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        plt.close("all")


class TestResidualPanel:
    def test_basic_render(self):
        wave, flux, err = make_wave_flux(n=100)
        model = flux + np.random.default_rng(1).normal(0, 0.002, 100)
        panel = ResidualPanel(
            wave, flux, 10.0, 20.0,
            model_flux_adj=model,
            error_data=err,
            error_adj=err,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert panel.ax is not None
        plt.close("all")

    def test_with_excluded_ranges(self):
        wave, flux, err = make_wave_flux(n=100)
        model = flux.copy()
        panel = ResidualPanel(
            wave, flux, 10.0, 20.0,
            model_flux_adj=model,
            error_data=err,
            error_adj=err,
            excluded_ranges=[(13.0, 14.0)],
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert len(panel.ax.patches) >= 1
        plt.close("all")

    def test_with_noise_floor(self):
        wave, flux, err = make_wave_flux(n=100)
        model = flux.copy()
        err_adj = np.sqrt(err ** 2 + 0.01 ** 2)
        panel = ResidualPanel(
            wave, flux, 10.0, 20.0,
            model_flux_adj=model,
            error_data=err,
            error_adj=err_adj,
            has_noise_floor=True,
            is_first_row=True,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert len(panel.ax.collections) >= 2
        plt.close("all")

    def test_compute_ylim(self):
        wave, flux, err = make_wave_flux(n=100)
        model = flux + 0.01
        panel = ResidualPanel(
            wave, flux, 10.0, 20.0,
            model_flux_adj=model,
            error_data=err,
            error_adj=err,
            ax=plt.subplots()[1],
        )
        ymin, ymax = panel.compute_ylim()
        assert ymin < 0
        assert ymax > 0
        assert abs(ymin) == pytest.approx(abs(ymax))
        plt.close("all")

    def test_gap_mode_skip_filters_nan_residuals(self):
        wave, flux, err = make_wave_flux(n=100, gap=(14.0, 16.0))
        model = np.where(np.isnan(flux), np.nan, flux + 0.001)
        panel = ResidualPanel(
            wave, flux, 10.0, 20.0,
            model_flux_adj=model,
            gap_mode=GapMode.SKIP,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        plt.close("all")
