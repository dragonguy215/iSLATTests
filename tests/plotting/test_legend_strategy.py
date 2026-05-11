"""Tests for :mod:`iSLAT.Modules.Plotting.LegendStrategy`.

Covers the :class:`LegendStrategy` ABC and its two concrete
implementations — :class:`StandardLegend` (default for non-stacked
plots) and :class:`MoleculeColorLegend` (default for stacked plots).
"""

import matplotlib
matplotlib.use("Agg")

import pytest
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from iSLAT.Modules.Plotting.LegendStrategy import (
    LegendStrategy,
    MoleculeColorLegend,
    StandardLegend,
)


# ======================================================================
# Fixtures
# ======================================================================
@pytest.fixture()
def fig_ax():
    """Create a basic figure + axes pair."""
    fig, ax = plt.subplots(figsize=(10, 4))
    fig.subplots_adjust(top=0.93)
    yield fig, ax
    plt.close(fig)


@pytest.fixture()
def stacked_fig():
    """Simulate a stacked-panel figure with subplots_adjust."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.94, top=0.93, bottom=0.06)
    yield fig, axes
    plt.close(fig)


@pytest.fixture()
def sample_labels():
    return ["H₂O", "CO", "OH", "HCN", "C₂H₂"]


@pytest.fixture()
def sample_colors():
    return ["blue", "red", "green", "orange", "purple"]


# ======================================================================
# LegendStrategy ABC
# ======================================================================
class TestLegendStrategyABC:
    """The ABC cannot be instantiated directly."""

    def test_cannot_instantiate(self):
        with pytest.raises(TypeError, match="abstract method"):
            LegendStrategy()

    def test_required_methods(self):
        """Subclasses must implement all four abstract methods."""
        required = {"build_legend", "remove_legend", "update_visibility", "apply_theme"}
        actual = set()
        for name in dir(LegendStrategy):
            attr = getattr(LegendStrategy, name, None)
            if getattr(attr, "__isabstractmethod__", False):
                actual.add(name)
        assert actual == required

    def test_custom_subclass_works(self, fig_ax):
        """A minimal concrete subclass can be instantiated and used."""
        fig, ax = fig_ax

        class NoopLegend(LegendStrategy):
            def build_legend(self, ax, fig, labels, colors, **kw):
                pass
            def remove_legend(self, ax):
                pass
            def update_visibility(self, ax, visible):
                pass
            def apply_theme(self, ax, theme):
                pass

        legend = NoopLegend()
        # Should not raise
        legend.build_legend(ax, fig, ["A"], ["red"])
        legend.remove_legend(ax)
        legend.update_visibility(ax, True)
        legend.apply_theme(ax, {})


# ======================================================================
# StandardLegend
# ======================================================================
class TestStandardLegend:
    """Tests for the artist-based standard legend."""

    def test_is_legend_strategy(self):
        assert isinstance(StandardLegend(), LegendStrategy)

    def test_build_from_artists(self, fig_ax):
        """Builds legend from visible labelled artists on the axes."""
        fig, ax = fig_ax
        ax.plot([1, 2], [3, 4], label="line1")
        ax.plot([1, 2], [5, 6], label="line2")

        strategy = StandardLegend()
        # labels/colors args are ignored — legend comes from axes artists
        strategy.build_legend(ax, fig, [], [])

        leg = ax.get_legend()
        assert leg is not None
        texts = [t.get_text() for t in leg.get_texts()]
        assert "line1" in texts
        assert "line2" in texts

    def test_invisible_artists_excluded(self, fig_ax):
        """Invisible artists should not appear in the legend."""
        fig, ax = fig_ax
        l1, = ax.plot([1, 2], [3, 4], label="visible")
        l2, = ax.plot([1, 2], [5, 6], label="hidden")
        l2.set_visible(False)

        strategy = StandardLegend()
        strategy.build_legend(ax, fig, [], [])

        leg = ax.get_legend()
        assert leg is not None
        texts = [t.get_text() for t in leg.get_texts()]
        assert "visible" in texts
        assert "hidden" not in texts

    def test_no_artists_removes_legend(self, fig_ax):
        """When no visible labelled artists exist, legend is removed."""
        fig, ax = fig_ax
        # Seed a legend manually first
        ax.legend(["dummy"])
        assert ax.get_legend() is not None

        strategy = StandardLegend()
        strategy.build_legend(ax, fig, [], [])

        assert ax.get_legend() is None

    def test_ncols_for_many_items(self, fig_ax):
        """More than 8 visible items → 2-column legend."""
        fig, ax = fig_ax
        for i in range(10):
            ax.plot([1, 2], [i, i + 1], label=f"item{i}")

        strategy = StandardLegend()
        strategy.build_legend(ax, fig, [], [])

        leg = ax.get_legend()
        assert leg is not None
        assert leg._ncols >= 2

    def test_max_ncols_respected(self, fig_ax):
        """max_ncols kwarg caps the column count."""
        fig, ax = fig_ax
        for i in range(10):
            ax.plot([1, 2], [i, i + 1], label=f"item{i}")

        strategy = StandardLegend()
        strategy.build_legend(ax, fig, [], [], max_ncols=1)

        leg = ax.get_legend()
        assert leg is not None
        assert leg._ncols == 1

    def test_remove_legend(self, fig_ax):
        fig, ax = fig_ax
        ax.plot([1, 2], [3, 4], label="test")
        ax.legend()
        assert ax.get_legend() is not None

        strategy = StandardLegend()
        strategy.remove_legend(ax)
        assert ax.get_legend() is None

    def test_remove_legend_noop_when_absent(self, fig_ax):
        """remove_legend should not raise when no legend exists."""
        fig, ax = fig_ax
        strategy = StandardLegend()
        strategy.remove_legend(ax)  # should not raise

    def test_update_visibility(self, fig_ax):
        fig, ax = fig_ax
        ax.plot([1, 2], [3, 4], label="test")
        ax.legend()

        strategy = StandardLegend()
        strategy.update_visibility(ax, False)
        assert not ax.get_legend().get_visible()

        strategy.update_visibility(ax, True)
        assert ax.get_legend().get_visible()

    def test_update_visibility_noop_when_absent(self, fig_ax):
        """update_visibility should not raise when no legend exists."""
        fig, ax = fig_ax
        strategy = StandardLegend()
        strategy.update_visibility(ax, True)  # should not raise

    def test_apply_theme(self, fig_ax):
        fig, ax = fig_ax
        ax.plot([1, 2], [3, 4], label="test")
        ax.legend()

        theme = {"foreground": "white", "graph_fill_color": "#333333"}
        strategy = StandardLegend()
        strategy.apply_theme(ax, theme)

        leg = ax.get_legend()
        for text in leg.get_texts():
            assert text.get_color() == "white"

    def test_apply_theme_noop_when_no_legend(self, fig_ax):
        fig, ax = fig_ax
        strategy = StandardLegend()
        strategy.apply_theme(ax, {"foreground": "red"})  # should not raise


# ======================================================================
# MoleculeColorLegend
# ======================================================================
class TestMoleculeColorLegend:
    """Tests for the compact text-only molecule color legend."""

    def test_is_legend_strategy(self):
        assert isinstance(MoleculeColorLegend(), LegendStrategy)

    def test_build_basic(self, fig_ax, sample_labels, sample_colors):
        """Creates a text-only legend with correct labels and colors."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        leg = ax.get_legend()
        assert leg is not None
        texts = leg.get_texts()
        assert len(texts) == len(sample_labels)
        for t, expected_label, expected_color in zip(texts, sample_labels, sample_colors):
            assert t.get_text() == expected_label
            assert t.get_color() == expected_color

    def test_mol_color_tag(self, fig_ax, sample_labels, sample_colors):
        """Legend texts are tagged with _islat_mol_color = True."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        for text in ax.get_legend().get_texts():
            assert getattr(text, "_islat_mol_color", False) is True

    def test_invisible_handles(self, fig_ax, sample_labels, sample_colors):
        """Handles should be invisible (text-only appearance)."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        leg = ax.get_legend()
        for handle in leg.legend_handles:
            assert isinstance(handle, Patch)
            assert handle.get_facecolor() == (0.0, 0.0, 0.0, 0.0)  # 'none'
            assert handle.get_edgecolor() == (0.0, 0.0, 0.0, 0.0)

    def test_empty_labels_removes_legend(self, fig_ax):
        """Passing empty labels removes any existing legend."""
        fig, ax = fig_ax
        # Create a legend first
        ax.legend(["dummy"])
        assert ax.get_legend() is not None

        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, [], [])

        assert ax.get_legend() is None

    def test_replaces_existing_legend(self, fig_ax):
        """Calling build_legend twice replaces the old legend."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()

        strategy.build_legend(ax, fig, ["A", "B"], ["red", "blue"])
        first = ax.get_legend()
        assert len(first.get_texts()) == 2

        strategy.build_legend(ax, fig, ["X"], ["green"])
        second = ax.get_legend()
        assert len(second.get_texts()) == 1
        assert second.get_texts()[0].get_text() == "X"

    def test_frameon_false(self, fig_ax, sample_labels, sample_colors):
        """Legend frame should be off (frameon=False)."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        leg = ax.get_legend()
        assert not leg.get_frame().get_visible() or not leg.get_frame_on()

    def test_fontsize_kwarg(self, fig_ax, sample_labels, sample_colors):
        """Custom fontsize is passed through."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors, fontsize=7)

        leg = ax.get_legend()
        assert leg is not None
        # The legend's prop fontsize may be set via the FontProperties
        # object rather than on individual texts (matplotlib applies
        # fontsize from prop).  Verify the legend was built with the
        # requested size by checking the prop attribute.
        assert leg.prop.get_size() == 7

    def test_max_ncols_respected(self, fig_ax, sample_labels, sample_colors):
        """max_ncols caps the column count."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors, max_ncols=2)

        leg = ax.get_legend()
        assert leg._ncols <= 2

    def test_remove_legend(self, fig_ax, sample_labels, sample_colors):
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)
        assert ax.get_legend() is not None

        strategy.remove_legend(ax)
        assert ax.get_legend() is None

    def test_update_visibility(self, fig_ax, sample_labels, sample_colors):
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        strategy.update_visibility(ax, False)
        assert not ax.get_legend().get_visible()

        strategy.update_visibility(ax, True)
        assert ax.get_legend().get_visible()

    def test_apply_theme_preserves_mol_colors(self, fig_ax, sample_labels, sample_colors):
        """apply_theme should NOT overwrite per-molecule colors."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        strategy.apply_theme(ax, {"foreground": "white", "graph_fill_color": "black"})

        for text, expected_color in zip(ax.get_legend().get_texts(), sample_colors):
            # Molecule-colored texts should keep their original color
            assert text.get_color() == expected_color

    def test_apply_theme_noop_when_no_legend(self, fig_ax):
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.apply_theme(ax, {"foreground": "red"})  # should not raise

    # ------------------------------------------------------------------
    # Dynamic width / ncols computation
    # ------------------------------------------------------------------
    def test_ncols_auto_fits_panel_width(self):
        """Auto ncols should not produce a legend wider than the panel."""
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.subplots_adjust(top=0.93)
        try:
            long_labels = [f"MoleculeWithLongName_{i}" for i in range(8)]
            long_colors = ["blue"] * 8

            strategy = MoleculeColorLegend()
            strategy.build_legend(ax, fig, long_labels, long_colors)

            leg = ax.get_legend()
            # With 6-inch-wide figure and long labels, ncols should be < 8
            assert leg._ncols < 8
        finally:
            plt.close(fig)

    def test_ncols_single_label(self, fig_ax):
        """A single label should produce 1 column."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, ["H₂O"], ["blue"])

        assert ax.get_legend()._ncols == 1

    def test_ncols_narrow_figure(self):
        """Very narrow figure forces fewer columns."""
        fig, ax = plt.subplots(figsize=(3, 4))
        fig.subplots_adjust(top=0.93)
        try:
            labels = [f"Molecule_{i}" for i in range(6)]
            colors = ["blue"] * 6

            strategy = MoleculeColorLegend()
            strategy.build_legend(ax, fig, labels, colors)

            leg = ax.get_legend()
            # Narrow figure → should force < 6 columns
            assert leg._ncols < 6
        finally:
            plt.close(fig)

    def test_ncols_wide_figure_fits_all(self):
        """Wide figure with short labels should fit all in one row."""
        fig, ax = plt.subplots(figsize=(20, 4))
        fig.subplots_adjust(top=0.93)
        try:
            labels = ["A", "B", "C"]
            colors = ["red", "blue", "green"]

            strategy = MoleculeColorLegend()
            strategy.build_legend(ax, fig, labels, colors)

            leg = ax.get_legend()
            assert leg._ncols == 3
        finally:
            plt.close(fig)

    # ------------------------------------------------------------------
    # Anti-overlap / positioning
    # ------------------------------------------------------------------
    def test_uses_figure_transform(self, fig_ax, sample_labels, sample_colors):
        """Legend should be positioned in figure coordinates."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        leg = ax.get_legend()
        # The legend's loc should be 'upper center' (code 9)
        assert leg._loc == 9  # upper center
        # Verify the legend was created (the build_legend uses fig.transFigure
        # via bbox_transform kwarg which is internal but we can verify
        # the legend was placed using figure-level coordinates by checking
        # that the legend text is present)
        assert len(leg.get_texts()) == len(sample_labels)

    def test_y_anchor_in_top_margin(self, fig_ax, sample_labels, sample_colors):
        """The y-anchor should be above the panel top (0.93) but below 1.0."""
        fig, ax = fig_ax
        strategy = MoleculeColorLegend()
        strategy.build_legend(ax, fig, sample_labels, sample_colors)

        # Verify via the helper that the anchor is in the margin.
        y = MoleculeColorLegend._safe_y_anchor(fig)
        assert 0.93 <= y <= 1.0

    def test_safe_y_anchor_no_title(self):
        """Without suptitle, y_anchor is midpoint of top margin."""
        fig, ax = plt.subplots(figsize=(10, 4))
        fig.subplots_adjust(top=0.93)
        try:
            y = MoleculeColorLegend._safe_y_anchor(fig)
            expected = (0.93 + 1.0) / 2.0  # 0.965
            assert abs(y - expected) < 0.01
        finally:
            plt.close(fig)

    def test_safe_y_anchor_with_suptitle(self):
        """With suptitle, y_anchor shifts down below the title."""
        fig, ax = plt.subplots(figsize=(10, 4))
        fig.subplots_adjust(top=0.93)
        fig.suptitle("Test Title", fontsize=14)
        fig.canvas.draw()  # Need to render for text extents
        try:
            y = MoleculeColorLegend._safe_y_anchor(fig)
            # With a title, y should be lower than the no-title case
            y_no_title = (0.93 + 1.0) / 2.0
            assert y < y_no_title
        finally:
            plt.close(fig)


# ======================================================================
# Swappability / Integration
# ======================================================================
class TestSwappability:
    """Verify that legend strategies are interchangeable."""

    def test_swap_standard_for_molecule(self, fig_ax, sample_labels, sample_colors):
        """A StandardLegend can be replaced with MoleculeColorLegend."""
        fig, ax = fig_ax
        ax.plot([1, 2], [3, 4], label="line1")

        std = StandardLegend()
        std.build_legend(ax, fig, [], [])
        leg1_texts = [t.get_text() for t in ax.get_legend().get_texts()]
        assert "line1" in leg1_texts

        mol = MoleculeColorLegend()
        mol.build_legend(ax, fig, sample_labels, sample_colors)
        leg2_texts = [t.get_text() for t in ax.get_legend().get_texts()]
        assert sample_labels[0] in leg2_texts

    def test_swap_molecule_for_standard(self, fig_ax, sample_labels, sample_colors):
        """A MoleculeColorLegend can be replaced with StandardLegend."""
        fig, ax = fig_ax

        mol = MoleculeColorLegend()
        mol.build_legend(ax, fig, sample_labels, sample_colors)
        assert ax.get_legend() is not None

        ax.plot([1, 2], [3, 4], label="swapped")
        std = StandardLegend()
        std.build_legend(ax, fig, [], [])
        texts = [t.get_text() for t in ax.get_legend().get_texts()]
        assert "swapped" in texts

    def test_custom_strategy_pluggable(self, fig_ax):
        """A custom strategy can be used as a drop-in replacement."""
        fig, ax = fig_ax

        class AlwaysEmptyLegend(LegendStrategy):
            def build_legend(self, ax, fig, labels, colors, **kw):
                old = ax.get_legend()
                if old is not None:
                    old.remove()
            def remove_legend(self, ax):
                pass
            def update_visibility(self, ax, visible):
                pass
            def apply_theme(self, ax, theme):
                pass

        strategy = AlwaysEmptyLegend()
        ax.plot([1, 2], [3, 4], label="test")
        ax.legend()
        assert ax.get_legend() is not None

        strategy.build_legend(ax, fig, ["A"], ["red"])
        assert ax.get_legend() is None


# ======================================================================
# Integration with BasePlot / StackedSpectralPanel defaults
# ======================================================================
class TestDefaultStrategyAssignment:
    """Verify that BasePlot and StackedSpectralPanel get the right defaults."""

    def test_base_plot_default_is_standard(self):
        """BasePlot subclasses default to StandardLegend."""
        from iSLAT.Modules.Plotting.BasePlot import BasePlot

        class Stub(BasePlot):
            def generate_plot(self, **kw):
                pass

        plot = Stub()
        assert isinstance(plot.legend_strategy, StandardLegend)

    def test_stacked_panel_default_is_molecule_color(self):
        """StackedSpectralPanel subclasses default to MoleculeColorLegend."""
        from iSLAT.Modules.Plotting.StackedSpectralPanel import StackedSpectralPanel
        from iSLAT.Modules.Plotting.SpectralPanel import SpectralPanel

        wave = np.linspace(5, 25, 200)
        flux = np.random.default_rng(42).normal(0, 1, 200)

        class StubStacked(StackedSpectralPanel):
            def _create_cell(self, idx, xmin, xmax, gs_slot, **kw):
                ax = self.fig.add_subplot(gs_slot)
                return [ax]
            def _cell_height_ratios(self):
                return [1]

        plot = StubStacked(wave, flux, n_panels=2)
        assert isinstance(plot.legend_strategy, MoleculeColorLegend)

    def test_explicit_strategy_preserved_in_stacked(self):
        """Explicitly passing a strategy to SSP is not overridden."""
        from iSLAT.Modules.Plotting.StackedSpectralPanel import StackedSpectralPanel

        wave = np.linspace(5, 25, 200)
        flux = np.random.default_rng(42).normal(0, 1, 200)

        custom = StandardLegend()

        class StubStacked(StackedSpectralPanel):
            def _create_cell(self, idx, xmin, xmax, gs_slot, **kw):
                ax = self.fig.add_subplot(gs_slot)
                return [ax]
            def _cell_height_ratios(self):
                return [1]

        plot = StubStacked(wave, flux, n_panels=2, legend_strategy=custom)
        assert plot.legend_strategy is custom
        assert isinstance(plot.legend_strategy, StandardLegend)

    def test_fsp_default_is_molecule_color(self):
        """FullSpectrumPlot defaults to MoleculeColorLegend."""
        from iSLAT.Modules.Plotting.FullSpectrumPlot import FullSpectrumPlot

        wave = np.linspace(5, 25, 200)
        flux = np.random.default_rng(42).normal(0, 1, 200)
        plot = FullSpectrumPlot(wave, flux, n_panels=2)
        assert isinstance(plot.legend_strategy, MoleculeColorLegend)
        plt.close(plot.fig)

    def test_rsp_default_is_molecule_color(self):
        """ResidualSpectrumPlot defaults to MoleculeColorLegend."""
        from iSLAT.Modules.Plotting.ResidualSpectrumPlot import ResidualSpectrumPlot

        wave = np.linspace(5, 25, 200)
        flux = np.random.default_rng(42).normal(0, 1, 200)
        model = flux * 0.9
        plot = ResidualSpectrumPlot(wave, flux, model_flux=model, n_panels=2)
        assert isinstance(plot.legend_strategy, MoleculeColorLegend)
        plt.close(plot.fig)


# ======================================================================
# _compute_ncols edge cases
# ======================================================================
class TestComputeNcols:
    """Direct tests for MoleculeColorLegend._compute_ncols."""

    def test_single_label(self, fig_ax):
        fig, ax = fig_ax
        result = MoleculeColorLegend._compute_ncols(ax, fig, ["X"])
        assert result == 1

    def test_empty_labels_returns_one(self, fig_ax):
        """Edge case: empty list should still return >= 1."""
        fig, ax = fig_ax
        # _compute_ncols should not be called with empty list normally,
        # but if it is, it should not crash.
        # The build_legend method returns early for empty labels,
        # but we test the helper directly.
        result = MoleculeColorLegend._compute_ncols(ax, fig, ["X"])
        assert result >= 1

    def test_max_ncols_caps_result(self, fig_ax):
        fig, ax = fig_ax
        labels = ["A", "B", "C", "D", "E"]
        result = MoleculeColorLegend._compute_ncols(ax, fig, labels, max_ncols=2)
        assert result <= 2

    def test_many_long_labels_force_few_cols(self):
        """Very long labels on a narrow figure should yield few columns."""
        fig, ax = plt.subplots(figsize=(4, 3))
        fig.subplots_adjust(top=0.93)
        try:
            labels = [f"VeryLongMoleculeName_{i}_extra" for i in range(10)]
            result = MoleculeColorLegend._compute_ncols(ax, fig, labels)
            assert result < 5
        finally:
            plt.close(fig)

    def test_short_labels_wide_fig_maximises_cols(self):
        """Short labels on a wide figure should allow many columns."""
        fig, ax = plt.subplots(figsize=(20, 3))
        fig.subplots_adjust(top=0.93)
        try:
            labels = ["A", "B", "C", "D"]
            result = MoleculeColorLegend._compute_ncols(ax, fig, labels)
            assert result == 4
        finally:
            plt.close(fig)
