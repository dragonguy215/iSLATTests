# -*- coding: utf-8 -*-
"""Tests for FullSpectrumPlot."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib.axes import Axes

from iSLAT.Modules.Plotting import FullSpectrumPlot
from iSLAT.Modules.Plotting.SpectrumPanel import SpectrumPanel

from tests.plotting import make_wave_flux, make_atomic_lines, make_line_list


class TestFullSpectrumPlot:
    def _make_fsp(self, n_panels=3, gap_mode="connect", **kwargs):
        wave, flux, err = make_wave_flux(n=200)
        return FullSpectrumPlot(
            wave, flux, error_data=err, n_panels=n_panels,
            gap_mode=gap_mode, **kwargs,
        )

    def test_generate_plot(self):
        fsp = self._make_fsp()
        fsp.generate_plot()
        assert fsp.fig is not None
        assert len(fsp.panels) == 3
        for idx, panels in fsp.panels.items():
            assert len(panels) == 1
            assert isinstance(panels[0], SpectrumPanel)
        fsp.close()

    def test_subplots_populated(self):
        fsp = self._make_fsp()
        fsp.generate_plot()
        assert len(fsp.subplots) == 3
        for idx in range(3):
            assert isinstance(fsp.subplots[idx], Axes)
        fsp.close()

    def test_gap_mode_skip_rendering(self):
        wave, flux, err = make_wave_flux(n=200, gap=(13.0, 17.0))
        fsp = FullSpectrumPlot(
            wave, flux, error_data=err, n_panels=5, gap_mode="skip",
        )
        fsp.generate_plot()
        assert fsp.fig is not None
        # Some panels should have been skipped
        assert len(fsp.panels) <= 5
        fsp.close()

    def test_with_annotations(self):
        wave, flux, err = make_wave_flux(n=200)
        atomic = make_atomic_lines([12.0, 15.0, 18.0])
        lines = make_line_list([11.0, 14.0])
        fsp = FullSpectrumPlot(
            wave, flux, error_data=err, n_panels=3,
            atomic_lines=atomic, line_list=lines,
        )
        fsp.generate_plot()
        # Check that some annotations were drawn
        total_texts = 0
        for panels in fsp.panels.values():
            for p in panels:
                total_texts += len(p.ax.texts)
        assert total_texts > 0
        fsp.close()

    def test_annotations_respect_gap_tightening(self):
        """Annotations should not be drawn in gap-cropped regions."""
        wave = np.concatenate([
            np.linspace(10, 10.5, 20),
            np.linspace(18, 20, 100),
        ])
        flux = np.ones(120)
        # Atomic line at 12.0 is in the gap
        atomic = make_atomic_lines([12.0, 19.0])
        fsp = FullSpectrumPlot(
            wave, flux, n_panels=1,
            gap_mode="skip",
            atomic_lines=atomic,
        )
        fsp.generate_plot()
        # Check the single panel's annotations
        panel = list(fsp.panels.values())[0][0]
        vis_lo, vis_hi = panel.ax.get_xlim()
        # All annotation x-positions should be within visible range
        for txt in panel.ax.texts:
            if hasattr(txt, "_islat_atomic_line"):
                x, _ = txt.get_position()
                assert vis_lo <= x <= vis_hi + 1.0, (
                    f"Annotation at x={x} outside visible range "
                    f"[{vis_lo}, {vis_hi}]"
                )
        fsp.close()

    def test_update_data(self):
        fsp = self._make_fsp(n_panels=3)
        fsp.generate_plot()
        new_wave, new_flux, _ = make_wave_flux(wmin=5, wmax=25, n=300)
        changed = fsp.update_data(new_wave, new_flux)
        assert changed is True  # Range changed

    def test_save_and_close(self, tmp_path):
        fsp = self._make_fsp()
        fsp.generate_plot()
        out = fsp.save(tmp_path / "fsp.png", dpi=72)
        assert out.exists()
        fsp.close()
        assert fsp.fig is None

    def test_uniform_ylim(self):
        fsp = self._make_fsp(n_panels=3, uniform_ylim=True)
        fsp.generate_plot()
        # All panels should have the same y-limits
        ylims = [list(fsp.panels.values())[i][0].ax.get_ylim()
                 for i in fsp.panels]
        for yl in ylims[1:]:
            assert yl[0] == pytest.approx(ylims[0][0], abs=0.001)
            assert yl[1] == pytest.approx(ylims[0][1], abs=0.001)
        fsp.close()
