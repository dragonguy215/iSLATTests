# -*- coding: utf-8 -*-
"""Tests for DisplayConfig — centralised high-DPI / rendering configuration."""

import matplotlib
matplotlib.use("Agg")

import platform
import pytest

from iSLAT.Modules.GUI.DisplayConfig import (
    _DisplayConfig,
    _build_config,
    apply_matplotlib_defaults,
    detect_scale_factor,
    display_config,
    _detect_macos_scale_factor,
    _detect_linux_scale_factor,
)


class TestDetectScaleFactor:
    """Platform scale-factor detection."""

    def test_returns_float(self):
        result = detect_scale_factor()
        assert isinstance(result, float)

    def test_returns_positive(self):
        assert detect_scale_factor() > 0

    def test_returns_reasonable_range(self):
        """Scale factor should be between 0.5 and 4.0 on any real display."""
        sf = detect_scale_factor()
        assert 0.5 <= sf <= 4.0

    def test_platform_specific_dispatch(self):
        """detect_scale_factor should dispatch to the correct platform function."""
        system = platform.system()
        sf = detect_scale_factor()
        if system == "Darwin":
            assert sf in (1.0, 2.0)  # macOS: standard or Retina
        else:
            assert sf >= 1.0

    def test_macos_detector_returns_1_or_2(self):
        """The macOS detector should return 1.0 or 2.0."""
        result = _detect_macos_scale_factor()
        assert result in (1.0, 2.0)

    def test_linux_detector_fallback(self):
        """Linux detector should fall back to 1.0 when xrdb is unavailable."""
        # On macOS/CI this will hit the fallback
        result = _detect_linux_scale_factor()
        assert isinstance(result, float)
        assert result >= 1.0


class TestBuildConfig:
    """_build_config factory and _DisplayConfig dataclass."""

    def test_returns_display_config(self):
        cfg = _build_config()
        assert isinstance(cfg, _DisplayConfig)

    def test_frozen_dataclass(self):
        cfg = _build_config()
        with pytest.raises(AttributeError):
            cfg.scale_factor = 99  # type: ignore[misc]

    def test_figure_dpi_reasonable(self):
        cfg = _build_config()
        assert 72 <= cfg.figure_dpi <= 200

    def test_savefig_dpi_high_quality(self):
        cfg = _build_config()
        assert cfg.savefig_dpi >= 200

    def test_agg_chunksize_large(self):
        cfg = _build_config()
        assert cfg.agg_chunksize >= 10_000

    def test_antialiasing_enabled(self):
        cfg = _build_config()
        assert cfg.text_antialiased is True
        assert cfg.lines_antialiased is True

    def test_is_hidpi_flag(self):
        cfg = _build_config()
        assert cfg.is_hidpi == (cfg.scale_factor >= 1.5)

    def test_tk_scaling_none_on_macos(self):
        """On macOS, tk_scaling should be None (Tk handles Retina natively)."""
        cfg = _build_config()
        if platform.system() == "Darwin":
            assert cfg.tk_scaling is None


class TestSingleton:
    """Module-level display_config is the eagerly-created singleton."""

    def test_singleton_exists(self):
        assert display_config is not None
        assert isinstance(display_config, _DisplayConfig)

    def test_singleton_matches_fresh_build(self):
        """Singleton values should be consistent with a fresh build."""
        fresh = _build_config()
        assert display_config.scale_factor == fresh.scale_factor
        assert display_config.figure_dpi == fresh.figure_dpi
        assert display_config.savefig_dpi == fresh.savefig_dpi


class TestApplyMatplotlibDefaults:
    """apply_matplotlib_defaults sets rcParams correctly."""

    def test_sets_figure_dpi(self):
        apply_matplotlib_defaults()
        assert matplotlib.rcParams["figure.dpi"] == display_config.figure_dpi

    def test_sets_savefig_dpi(self):
        apply_matplotlib_defaults()
        assert matplotlib.rcParams["savefig.dpi"] == display_config.savefig_dpi

    def test_sets_agg_chunksize(self):
        apply_matplotlib_defaults()
        assert matplotlib.rcParams["agg.path.chunksize"] == display_config.agg_chunksize

    def test_sets_font_type(self):
        apply_matplotlib_defaults()
        assert matplotlib.rcParams["pdf.fonttype"] == 42
        assert matplotlib.rcParams["ps.fonttype"] == 42

    def test_idempotent(self):
        """Calling apply_matplotlib_defaults twice should not change state."""
        apply_matplotlib_defaults()
        dpi1 = matplotlib.rcParams["figure.dpi"]
        apply_matplotlib_defaults()
        dpi2 = matplotlib.rcParams["figure.dpi"]
        assert dpi1 == dpi2

    def test_path_simplify_on(self):
        apply_matplotlib_defaults()
        assert matplotlib.rcParams["path.simplify"] is True
