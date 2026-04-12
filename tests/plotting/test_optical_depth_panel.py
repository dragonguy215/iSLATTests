# -*- coding: utf-8 -*-
"""Tests for OpticalDepthPanel and OpticalDepthSpectrumPlot."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.axes import Axes

from iSLAT.Modules.Plotting import OpticalDepthPanel, OpticalDepthSpectrumPlot

from tests.plotting import make_wave_flux


def _make_tau_cache(wave, n_molecules=1):
    """Build a synthetic mol_tau_cache for testing."""
    colors = ["blue", "red", "green", "orange", "purple"]
    cache = []
    for i in range(n_molecules):
        tau = 0.1 * (i + 1) * np.abs(np.sin(
            2 * np.pi * (wave - wave[0]) / (wave[-1] - wave[0]) * (i + 1)
        ))
        cache.append((
            wave, tau,
            colors[i % len(colors)],
            f"Mol{i}",
            f"mol_{i}",
        ))
    return cache


class TestOpticalDepthPanel:
    """Tests for the single-panel OpticalDepthPanel."""

    def _make_panel(self, display_mode="profile", log_scale=False,
                    n_molecules=1, **kwargs):
        wave, flux, err = make_wave_flux(n=200)
        tau_cache = _make_tau_cache(wave, n_molecules=n_molecules)
        return OpticalDepthPanel(
            wave, flux,
            xmin=float(wave[0]), xmax=float(wave[-1]),
            error_data=err,
            mol_tau_cache=tau_cache,
            display_mode=display_mode,
            log_scale=log_scale,
            **kwargs,
        )

    def test_generate_profile_mode(self):
        panel = self._make_panel(display_mode="profile")
        panel.generate_plot()
        assert panel.ax is not None
        panel.close()

    def test_generate_stem_mode(self):
        """Stem mode with no molecules still renders without error."""
        wave, flux, err = make_wave_flux(n=200)
        panel = OpticalDepthPanel(
            wave, flux,
            xmin=float(wave[0]), xmax=float(wave[-1]),
            display_mode="stem",
        )
        panel.generate_plot()
        assert panel.ax is not None
        panel.close()

    def test_multiple_molecules_profile(self):
        panel = self._make_panel(n_molecules=3)
        panel.generate_plot()
        # Should have fills rendered (check that axes has collections)
        assert len(panel.ax.collections) > 0 or len(panel.ax.lines) > 0
        panel.close()

    def test_stacked_fill_false(self):
        panel = self._make_panel(n_molecules=2, stacked_fill=False)
        panel.generate_plot()
        # Should have line artists instead of fills
        assert panel.ax is not None
        panel.close()

    def test_show_total_line(self):
        panel = self._make_panel(n_molecules=2, show_total=True)
        panel.generate_plot()
        assert panel.ax is not None
        panel.close()

    def test_log_scale(self):
        panel = self._make_panel(log_scale=True)
        panel.generate_plot()
        assert panel.ax.get_yscale() == "log"
        panel.close()

    def test_linear_scale(self):
        panel = self._make_panel(log_scale=False)
        panel.generate_plot()
        assert panel.ax.get_yscale() == "linear"
        panel.close()

    def test_compute_ylim_default(self):
        panel = self._make_panel()
        ymin, ymax = panel.compute_ylim()
        assert ymin == 0.0
        assert ymax > 0.0

    def test_compute_ylim_log(self):
        panel = self._make_panel(log_scale=True)
        ymin, ymax = panel.compute_ylim()
        assert ymin > 0.0  # Floor for log scale
        assert ymax > ymin

    def test_compute_ylim_empty_cache(self):
        wave, flux, err = make_wave_flux(n=200)
        panel = OpticalDepthPanel(
            wave, flux,
            xmin=float(wave[0]), xmax=float(wave[-1]),
            mol_tau_cache=[],
        )
        ymin, ymax = panel.compute_ylim()
        assert ymax > 0  # Should provide a sensible default

    def test_ylabel_set(self):
        panel = self._make_panel()
        panel.generate_plot()
        ylabel = panel.ax.get_ylabel()
        assert "\u03c4" in ylabel or "Optical" in ylabel
        panel.close()


class TestOpticalDepthSpectrumPlot:
    """Tests for the multi-panel stacked OpticalDepthSpectrumPlot."""

    def _make_odsp(self, n_panels=3, **kwargs):
        wave, flux, err = make_wave_flux(n=200)
        return OpticalDepthSpectrumPlot(
            wave, flux, error_data=err, n_panels=n_panels, **kwargs,
        )

    def test_generate_plot(self):
        odsp = self._make_odsp()
        odsp.generate_plot()
        assert odsp.fig is not None
        assert len(odsp.panels) == 3
        for idx, panels in odsp.panels.items():
            assert len(panels) == 1
            assert isinstance(panels[0], OpticalDepthPanel)
        odsp.close()

    def test_subplots_populated(self):
        odsp = self._make_odsp()
        odsp.generate_plot()
        assert len(odsp.subplots) == 3
        for idx in range(3):
            assert isinstance(odsp.subplots[idx], Axes)
        odsp.close()

    def test_panel_count(self):
        for n in (2, 5, 8):
            odsp = self._make_odsp(n_panels=n)
            odsp.generate_plot()
            assert len(odsp.panels) == n
            odsp.close()

    def test_log_scale_propagated(self):
        odsp = self._make_odsp(log_scale=True)
        odsp.generate_plot()
        for panels in odsp.panels.values():
            assert panels[0].ax.get_yscale() == "log"
        odsp.close()

    def test_display_mode_propagated(self):
        odsp = self._make_odsp(display_mode="stem")
        odsp.generate_plot()
        for panels in odsp.panels.values():
            assert panels[0].display_mode == "stem"
        odsp.close()

    def test_save_and_close(self, tmp_path):
        odsp = self._make_odsp()
        odsp.generate_plot()
        out = odsp.save(tmp_path / "odsp.png", dpi=72)
        assert out.exists()
        odsp.close()
        assert odsp.fig is None

    def test_uniform_ylim(self):
        odsp = self._make_odsp(n_panels=3, uniform_ylim=True)
        odsp.generate_plot()
        ylims = [
            list(odsp.panels.values())[i][0].ax.get_ylim()
            for i in odsp.panels
        ]
        for yl in ylims[1:]:
            assert yl == pytest.approx(ylims[0], abs=1e-6)
        odsp.close()

    def test_composite_stacking(self):
        """OpticalDepthSpectrumPlot + FullSpectrumPlot via stack_with."""
        from iSLAT.Modules.Plotting import FullSpectrumPlot
        wave, flux, err = make_wave_flux(n=200)
        fsp = FullSpectrumPlot(wave, flux, error_data=err, n_panels=3)
        odsp = OpticalDepthSpectrumPlot(
            wave, flux, error_data=err, n_panels=3,
        )
        composite = fsp + odsp
        composite.generate_plot()
        assert composite.fig is not None
        # Should have panels from both sources
        assert len(composite.panels) >= 3
        composite.close()
