# -*- coding: utf-8 -*-
"""Tests for BasePlot class methods and infrastructure."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure as MplFigure

from iSLAT.Modules.Plotting import BasePlot, DEFAULT_THEME

from tests.plotting import ConcreteSpectralPanel, make_wave_flux, make_atomic_lines, make_line_list


class TestBasePlotTheme:
    """Theme helper tests."""

    def test_default_theme_keys(self):
        assert "foreground" in DEFAULT_THEME
        assert "background" in DEFAULT_THEME
        assert "summed_spectra_color" in DEFAULT_THEME

    def test_get_theme_value_returns_from_theme(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            theme={"foreground": "red"},
        )
        assert panel._get_theme_value("foreground") == "red"

    def test_get_theme_value_default(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        assert panel._get_theme_value("nonexistent_key", "fallback") == "fallback"


class TestBasePlotFigureLifecycle:
    """Figure creation, ownership, and teardown."""

    def test_ensure_figure_creates_figure(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        assert panel.fig is None
        panel._ensure_figure()
        assert panel.fig is not None
        assert isinstance(panel.fig, MplFigure)
        panel.close()

    def test_ensure_figure_respects_figsize(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            figsize=(6, 3),
        )
        panel._ensure_figure()
        w, h = panel.fig.get_size_inches()
        assert abs(w - 6.0) < 0.01
        assert abs(h - 3.0) < 0.01
        panel.close()

    def test_close_releases_figure(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel._ensure_figure()
        panel.close()
        assert panel.fig is None

    def test_get_figure_returns_fig(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel._ensure_figure()
        assert panel.get_figure() is panel.fig
        panel.close()

    def test_external_figure_not_closed(self):
        fig = MplFigure(figsize=(4, 2))
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            fig=fig,
        )
        assert panel.fig is fig
        assert not panel._owns_figure
        panel.close()

    def test_apply_theme_sets_colors(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
            theme={"foreground": "#111111", "background": "#eeeeee",
                    "graph_fill_color": "#dddddd"},
        )
        panel.generate_plot()
        panel.apply_theme_to_figure()
        assert panel.fig.get_facecolor() is not None
        panel.close()


class TestClearTaggedArtists:
    def test_removes_lines(self):
        fig, ax = plt.subplots()
        line, = ax.plot([0, 1], [0, 1])
        line._test_tag = True
        assert len(ax.lines) == 1
        BasePlot._clear_tagged_artists(ax, "_test_tag", lines=True)
        assert len(ax.lines) == 0
        plt.close(fig)

    def test_removes_texts(self):
        fig, ax = plt.subplots()
        txt = ax.text(0.5, 0.5, "hello")
        txt._test_tag = True
        assert len(ax.texts) == 1
        BasePlot._clear_tagged_artists(ax, "_test_tag", texts=True)
        assert len(ax.texts) == 0
        plt.close(fig)

    def test_ignores_untagged(self):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])  # no tag
        BasePlot._clear_tagged_artists(ax, "_test_tag", lines=True)
        assert len(ax.lines) == 1
        plt.close(fig)


class TestBasePlotRenderingHelpers:
    """_plot_observed_spectrum, _plot_summed_spectrum, annotations."""

    def test_plot_observed_spectrum(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel.generate_plot()
        ax = panel.ax
        n_before = len(ax.lines)
        panel._plot_observed_spectrum(ax, np.arange(5.0), np.ones(5))
        assert len(ax.lines) > n_before
        panel.close()

    def test_plot_observed_spectrum_with_error(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel.generate_plot()
        ax = panel.ax
        panel._plot_observed_spectrum(
            ax, np.arange(5.0), np.ones(5), error_data=np.full(5, 0.01),
        )
        panel.close()

    def test_plot_observed_spectrum_deduplicate(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel.generate_plot()
        ax = panel.ax
        panel._plot_observed_spectrum(
            ax, np.arange(5.0), np.ones(5), deduplicate=True,
        )
        tagged_before = sum(1 for l in ax.lines if hasattr(l, '_islat_observed'))
        panel._plot_observed_spectrum(
            ax, np.arange(5.0), np.ones(5), deduplicate=True,
        )
        tagged_after = sum(1 for l in ax.lines if hasattr(l, '_islat_observed'))
        assert tagged_after == tagged_before
        panel.close()

    def test_plot_summed_spectrum(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel.generate_plot()
        ax = panel.ax
        panel._plot_summed_spectrum(
            ax, np.arange(5.0), np.array([0, 0.1, 0.2, 0.1, 0]),
        )
        tagged = [c for c in ax.collections if hasattr(c, '_islat_summed')]
        assert len(tagged) == 1
        panel.close()

    def test_plot_summed_spectrum_empty_data(self):
        panel = ConcreteSpectralPanel(
            np.arange(5.0), np.ones(5), 0.0, 4.0,
        )
        panel.generate_plot()
        panel._plot_summed_spectrum(panel.ax, np.array([]), np.array([]))
        panel.close()

    def test_plot_line_annotations(self):
        fig, ax = plt.subplots()
        ax.set_xlim(10, 20)
        ax.set_ylim(0, 1)
        df = make_line_list([12.0, 15.0, 25.0])
        BasePlot._plot_line_annotations(ax, df, (10, 20), 0, 1, tag="_test")
        tagged_texts = [t for t in ax.texts if hasattr(t, "_test")]
        assert len(tagged_texts) == 2
        plt.close(fig)

    def test_plot_atomic_lines(self):
        fig, ax = plt.subplots()
        ax.set_xlim(10, 20)
        ax.set_ylim(0, 1)
        df = make_atomic_lines([11.0, 15.0, 22.0])
        BasePlot._plot_atomic_lines(ax, df, xr=(10, 20), tag="_atag")
        tagged = [t for t in ax.texts if hasattr(t, "_atag")]
        assert len(tagged) == 2
        plt.close(fig)

    def test_plot_saved_line_markers(self):
        fig, ax = plt.subplots()
        df = pd.DataFrame({"lam": [5.0, 10.0], "xmin": [4.5, 9.5], "xmax": [5.5, 10.5]})
        BasePlot._plot_saved_line_markers(ax, df, tag="_stag")
        tagged_lines = [l for l in ax.lines if hasattr(l, "_stag")]
        assert len(tagged_lines) == 6
        plt.close(fig)


class TestBuildMoleculeLegend:
    def test_empty(self):
        fig, ax = plt.subplots()
        BasePlot.build_molecule_legend(ax, [], [])
        assert ax.get_legend() is None
        plt.close(fig)

    def test_creates_legend(self):
        fig, ax = plt.subplots()
        BasePlot.build_molecule_legend(
            ax, ["H2O", "CO"], ["blue", "red"],
            use_figure_transform=False,
        )
        leg = ax.get_legend()
        assert leg is not None
        texts = leg.get_texts()
        assert len(texts) == 2
        assert texts[0].get_text() == "H2O"
        plt.close(fig)


class TestBasePlotIO:
    """Save and IPython rich-display."""

    def test_save(self, tmp_path):
        wave, flux, _ = make_wave_flux(n=50)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        panel.generate_plot()
        out = panel.save(tmp_path / "test.png", dpi=72)
        assert out.exists()
        assert out.stat().st_size > 0
        panel.close()

    def test_repr_png(self):
        wave, flux, _ = make_wave_flux(n=50)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        panel.generate_plot()
        png = panel._repr_png_()
        assert png is not None
        assert len(png) > 100
        panel.close()

    def test_repr_html(self):
        wave, flux, _ = make_wave_flux(n=50)
        panel = ConcreteSpectralPanel(wave, flux, 10.0, 20.0)
        panel.generate_plot()
        html = panel._repr_html_()
        assert html is not None
        assert "<img" in html
        panel.close()
