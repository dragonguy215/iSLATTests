# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.LineSaveService.

Tests line info creation, formatting, and selection bounds logic.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock

from iSLAT.Modules.DataProcessing.LineSaveService import LineSaveService


# ============================================================================
# Construction
# ============================================================================

class TestLineSaveServiceInit:
    """Tests for LineSaveService initialization."""

    def test_can_instantiate(self):
        svc = LineSaveService()
        assert svc is not None


# ============================================================================
# create_default_line_info
# ============================================================================

class TestCreateDefaultLineInfo:
    """Tests for default line info creation."""

    def test_basic_creation(self):
        svc = LineSaveService()
        info = svc.create_default_line_info(center_wave=15.0, line_flux=1e-15, line_err=1e-17)

        assert info['lam'] == 15.0
        assert info['wavelength'] == 15.0
        assert info['flux'] == 1e-15
        assert info['intensity'] == 1e-15
        assert info['up_lev'] == 'Unknown'
        assert info['low_lev'] == 'Unknown'
        assert info['tau'] == 0.0

    def test_default_err_zero(self):
        svc = LineSaveService()
        info = svc.create_default_line_info(center_wave=10.0, line_flux=5e-16)
        # line_err defaults to 0.0; result should still work
        assert info['lam'] == 10.0

    def test_values_propagated(self):
        svc = LineSaveService()
        info = svc.create_default_line_info(center_wave=7.5, line_flux=3.14e-15, line_err=1e-17)
        assert info['inten'] == 3.14e-15
        assert info['e'] == 0.0
        assert info['a'] == 0.0
        assert info['g'] == 1.0


# ============================================================================
# format_line_for_save
# ============================================================================

class TestFormatLineForSave:
    """Tests for formatting line info for file output."""

    def test_basic_format(self):
        svc = LineSaveService()
        line_info = {
            'lam': 15.0,
            'up_lev': '0|10',
            'low_lev': '0|9',
            'tau': 0.5,
            'inten': 1e-15,
            'a_stein': 0.01,
            'e_up': 4586.0,
            'e_low': 3379.0,
            'g_up': 21,
            'g_low': 19,
        }
        result = svc.format_line_for_save(line_info, 'H2O', xmin=14.5, xmax=15.5)

        assert result['species'] == 'H2O'
        assert result['lam'] == 15.0
        assert result['lev_up'] == '0|10'
        assert result['lev_low'] == '0|9'
        assert result['xmin'] == 14.5
        assert result['xmax'] == 15.5
        assert result['a_stein'] == 0.01
        assert result['e_up'] == 4586.0
        assert result['g_up'] == 21

    def test_missing_keys_use_defaults(self):
        svc = LineSaveService()
        line_info = {'lam': 12.0}
        result = svc.format_line_for_save(line_info, 'CO', xmin=11.0, xmax=13.0)

        assert result['species'] == 'CO'
        assert result['lam'] == 12.0
        assert result['lev_up'] == ''
        assert result['lev_low'] == ''
        assert result['tau'] == 0.0
        assert result['a_stein'] == 0.0
        assert result['g_up'] == 1.0
        assert result['g_low'] == 1.0


# ============================================================================
# get_selection_bounds
# ============================================================================

class TestGetSelectionBounds:
    """Tests for selection bounds logic."""

    def test_uses_selected_wave_when_available(self):
        svc = LineSaveService()
        selected_wave = np.array([12.0, 12.5, 13.0, 13.5, 14.0])
        xmin, xmax = svc.get_selection_bounds(
            selected_wave, current_selection=(10.0, 20.0), line_wavelength=13.0
        )
        assert xmin == 12.0
        assert xmax == 14.0

    def test_falls_back_to_current_selection(self):
        svc = LineSaveService()
        xmin, xmax = svc.get_selection_bounds(
            selected_wave=None, current_selection=(11.0, 19.0), line_wavelength=15.0
        )
        assert xmin == 11.0
        assert xmax == 19.0

    def test_empty_array_falls_back(self):
        svc = LineSaveService()
        xmin, xmax = svc.get_selection_bounds(
            selected_wave=np.array([]), current_selection=(10.0, 20.0), line_wavelength=15.0
        )
        assert xmin == 10.0
        assert xmax == 20.0

    def test_single_element_wave(self):
        svc = LineSaveService()
        selected_wave = np.array([15.0])
        xmin, xmax = svc.get_selection_bounds(
            selected_wave, current_selection=(10.0, 20.0), line_wavelength=15.0
        )
        assert xmin == 15.0
        # xmax uses fallback for single-element array
        assert xmax == pytest.approx(15.0 + 0.01)


# ============================================================================
# extract_line_info_from_selection
# ============================================================================

class TestExtractLineInfoFromSelection:
    """Tests for extracting line info from plot selection."""

    def test_no_selection_returns_error(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = None

        result, error = svc.extract_line_info_from_selection(main_plot)
        assert result is None
        assert "No region selected" in error

    def test_no_current_selection_attr(self):
        svc = LineSaveService()
        main_plot = MagicMock(spec=[])  # no attributes at all

        result, error = svc.extract_line_info_from_selection(main_plot)
        assert result is None
        assert error != ""

    def test_invalid_save_type_returns_error(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = (10.0, 20.0)

        result, error = svc.extract_line_info_from_selection(main_plot, save_type="invalid")
        assert result is None
        assert "Invalid save type" in error

    def test_strongest_no_line_found(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = (10.0, 20.0)
        main_plot.find_strongest_line_from_data.return_value = None

        result, error = svc.extract_line_info_from_selection(main_plot, save_type="strongest")
        assert result is None
        assert "No valid line" in error

    def test_strongest_returns_line_info(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = (14.0, 16.0)
        expected_info = {'lam': 15.0, 'up_lev': '0|10', 'low_lev': '0|9'}
        main_plot.find_strongest_line_from_data.return_value = expected_info

        result, error = svc.extract_line_info_from_selection(main_plot, save_type="strongest")
        assert result is expected_info
        assert error == ""

    def test_selected_with_selected_line(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = (14.0, 16.0)
        expected_info = {'lam': 15.0, 'inten': 1e-15}
        main_plot.selected_line = expected_info

        result, error = svc.extract_line_info_from_selection(main_plot, save_type="selected")
        assert result is expected_info
        assert error == ""

    def test_selected_fallback_creates_default(self):
        svc = LineSaveService()
        main_plot = MagicMock()
        main_plot.current_selection = (14.0, 16.0)
        main_plot.selected_line = None
        main_plot.flux_integral.return_value = (1e-15, 1e-17)

        wave = np.linspace(10, 20, 500)
        flux = np.ones_like(wave) * 1e-16
        err = np.ones_like(wave) * 1e-18

        result, error = svc.extract_line_info_from_selection(
            main_plot, save_type="selected",
            wave_data=wave, flux_data=flux, err_data=err,
        )
        assert result is not None
        assert error == ""
        assert result['lam'] == 15.0  # midpoint of (14, 16)
