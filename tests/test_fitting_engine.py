# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.FittingEngine.

Tests Gaussian fitting, parameter extraction, and result formatting.
"""

import pytest
import numpy as np

import iSLAT.Constants as c
from iSLAT.Modules.DataProcessing.FittingEngine import FittingEngine


# ============================================================================
# Helpers
# ============================================================================

def _gaussian(x, amplitude, center, sigma):
    """Pure Gaussian for generating test data."""
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


# ============================================================================
# Construction
# ============================================================================

class TestFittingEngineInit:
    """Tests for FittingEngine initialization."""

    def test_default_state(self):
        fe = FittingEngine()
        assert fe.last_fit_result is None
        assert fe.last_fit_params is None
        assert isinstance(FittingEngine.VERBOSE_FIT_OUTPUT, bool)


# ============================================================================
# Single Gaussian fit
# ============================================================================

class TestSingleGaussianFit:
    """Tests for _fit_single_gaussian via fit_gaussian_line."""

    def test_fits_clean_gaussian(self):
        """Should recover parameters from a noise-free Gaussian."""
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 200)
        sigma = 0.05
        amplitude = 10.0
        y = _gaussian(x, amplitude, 15.0, sigma)

        result, fitted_wave, fitted_flux = fe.fit_gaussian_line(
            wave_data=x, flux_data=y, xmin=14.5, xmax=15.5
        )

        assert result is not None
        assert result.success
        assert fitted_wave is not None
        assert fitted_flux is not None
        assert len(fitted_flux) == len(x)

        # Recovered center should be near 15.0
        assert result.params['center'].value == pytest.approx(15.0, abs=0.01)

    def test_fits_noisy_gaussian(self):
        """Should fit reasonably despite moderate noise."""
        rng = np.random.default_rng(42)
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 300)
        y = _gaussian(x, 10.0, 15.0, 0.05) + rng.normal(0, 0.3, len(x))

        result, fitted_wave, fitted_flux = fe.fit_gaussian_line(
            wave_data=x, flux_data=y, xmin=14.5, xmax=15.5
        )

        assert result is not None
        assert result.params['center'].value == pytest.approx(15.0, abs=0.05)

    def test_stores_last_fit(self):
        """After fitting, last_fit_result and last_fit_params should be populated."""
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 200)
        y = _gaussian(x, 10.0, 15.0, 0.05)

        fe.fit_gaussian_line(wave_data=x, flux_data=y)

        assert fe.last_fit_result is not None
        assert fe.last_fit_params is not None

    def test_with_error_weights(self):
        """Providing error data should produce weighted fit without crashing."""
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 200)
        y = _gaussian(x, 10.0, 15.0, 0.05)
        err = np.full_like(y, 0.1)

        result, _, _ = fe.fit_gaussian_line(
            wave_data=x, flux_data=y, xmin=14.5, xmax=15.5, err_data=err
        )

        assert result is not None
        assert result.success

    def test_flat_data_still_returns(self):
        """Flat data (no line) should still return without crashing."""
        fe = FittingEngine()
        x = np.linspace(10, 20, 200)
        y = np.full_like(x, 1e-16)

        result, fitted_wave, fitted_flux = fe.fit_gaussian_line(
            wave_data=x, flux_data=y
        )

        # May or may not "succeed" but should not crash
        # result could be None or a result object
        assert fitted_wave is None or len(fitted_wave) == len(x)


# ============================================================================
# extract_line_parameters
# ============================================================================

class TestExtractLineParameters:
    """Tests for extracting parameters from last fit."""

    def test_no_fit_returns_empty(self):
        fe = FittingEngine()
        params = fe.extract_line_parameters()
        assert params == {}

    def test_single_gaussian_extraction(self):
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 300)
        y = _gaussian(x, 10.0, 15.0, 0.05)
        fe.fit_gaussian_line(wave_data=x, flux_data=y, xmin=14.5, xmax=15.5)

        params = fe.extract_line_parameters(rest_wavelength=15.0)

        assert 'center' in params
        assert 'fwhm' in params
        assert 'area' in params
        assert 'doppler_shift' in params
        assert params['center'] == pytest.approx(15.0, abs=0.01)

    def test_doppler_shift_calculation(self):
        """Doppler shift should be (center - rest) / rest * c_km/s."""
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 300)
        y = _gaussian(x, 10.0, 15.0, 0.05)
        fe.fit_gaussian_line(wave_data=x, flux_data=y, xmin=14.5, xmax=15.5)

        params = fe.extract_line_parameters(rest_wavelength=15.0)
        center = params['center']
        expected_doppler = (center - 15.0) / 15.0 * c.SPEED_OF_LIGHT_KMS
        assert params['doppler_shift'] == pytest.approx(expected_doppler, abs=1.0)


# ============================================================================
# _extract_component_parameters
# ============================================================================

class TestExtractComponentParameters:
    """Tests for the component parameter extraction helper."""

    def test_extracts_from_single_fit(self):
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 300)
        y = _gaussian(x, 10.0, 15.0, 0.05)
        fe.fit_gaussian_line(wave_data=x, flux_data=y, xmin=14.5, xmax=15.5)

        comp = fe._extract_component_parameters(
            fe.last_fit_params, prefix='', rest_wavelength=15.0
        )

        assert 'center' in comp
        assert 'sigma' in comp
        assert 'fwhm' in comp
        assert 'area' in comp
        assert 'fit_detected' in comp
        assert comp['fit_detected'] is True or comp['fit_detected'] is False or isinstance(comp['fit_detected'], (bool, np.bool_))


# ============================================================================
# format_fit_results_for_csv
# ============================================================================

class TestFormatFitResultsForCSV:
    """Tests for the CSV result formatter."""

    def test_successful_fit_format(self):
        fe = FittingEngine()
        x = np.linspace(14.5, 15.5, 300)
        y = _gaussian(x, 10.0, 15.0, 0.05)
        result, _, _ = fe.fit_gaussian_line(wave_data=x, flux_data=y, xmin=14.5, xmax=15.5)

        line_info = {
            'species': 'H2O', 'lev_up': '0|10', 'lev_low': '0|9',
            'lam': 15.0, 'a_stein': 0.01, 'e_up': 4586.0,
            'e_low': 3379.0, 'g_up': 21, 'g_low': 19,
        }

        entry = fe.format_fit_results_for_csv(
            fit_result=result,
            flux_data_integral=1e-15,
            err_data_integral=1e-17,
            xmin=14.5, xmax=15.5,
            rest_wavelength=15.0,
            line_info=line_info,
        )

        assert entry['species'] == 'H2O'
        assert entry['xmin'] == 14.5
        assert entry['xmax'] == 15.5
        assert 'Flux_data' in entry
        assert 'Err_data' in entry
        assert 'Line_SN' in entry
        assert 'Flux_fit' in entry
        assert entry['Fit_success'] is True

    def test_failed_fit_format(self):
        fe = FittingEngine()

        line_info = {
            'species': 'CO', 'lev_up': '0|2', 'lev_low': '0|1',
            'lam': 12.0, 'a_stein': 0.01, 'e_up': 2000.0,
            'e_low': 1000.0, 'g_up': 5, 'g_low': 3,
        }

        entry = fe.format_fit_results_for_csv(
            fit_result=None,
            flux_data_integral=5e-16,
            err_data_integral=2e-17,
            xmin=11.5, xmax=12.5,
            rest_wavelength=12.0,
            line_info=line_info,
        )

        assert entry['species'] == 'CO'
        assert entry['Fit_success'] is False
        assert np.isnan(entry['Flux_fit'])
        assert np.isnan(entry['FWHM_fit'])

    def test_detection_flag_below_threshold(self):
        """Line below detection threshold should have Line_det=False."""
        fe = FittingEngine()
        line_info = {'species': 'OH', 'lam': 10.0}
        # flux < sig_det_lim * error → not detected
        entry = fe.format_fit_results_for_csv(
            fit_result=None,
            flux_data_integral=1e-18,
            err_data_integral=1e-17,
            xmin=9.5, xmax=10.5,
            rest_wavelength=10.0,
            line_info=line_info,
            sig_det_lim=2,
        )
        assert entry['Line_det'] is False

    def test_detection_flag_above_threshold(self):
        """Line well above detection threshold should have Line_det=True."""
        fe = FittingEngine()
        line_info = {'species': 'OH', 'lam': 10.0}
        entry = fe.format_fit_results_for_csv(
            fit_result=None,
            flux_data_integral=1e-14,
            err_data_integral=1e-17,
            xmin=9.5, xmax=10.5,
            rest_wavelength=10.0,
            line_info=line_info,
            sig_det_lim=2,
        )
        assert entry['Line_det'] is True


# ============================================================================
# _estimate_components_from_user_selection
# ============================================================================

class TestEstimateComponents:
    """Tests for component estimation from line lists."""

    def test_single_strong_line(self):
        """A single strong line should give 1 component."""
        fe = FittingEngine()
        x = np.linspace(14, 16, 500)
        y = _gaussian(x, 10.0, 15.0, 0.05)

        mock_line = type('Line', (), {'lam': 15.0})()
        lines_with_intensity = [(mock_line, 10.0, 0.0)]

        n, centers = fe._estimate_components_from_user_selection(
            x, y, xmin=14.0, xmax=16.0,
            lines_with_intensity=lines_with_intensity,
            line_threshold=0.03,
        )

        assert n == 1
        assert len(centers) == 1
        assert abs(centers[0] - 15.0) < 0.01

    def test_two_strong_lines(self):
        """Two lines above threshold should give 2 components."""
        fe = FittingEngine()
        x = np.linspace(14, 16, 500)
        y = _gaussian(x, 10.0, 14.5, 0.05) + _gaussian(x, 8.0, 15.5, 0.05)

        line1 = type('Line', (), {'lam': 14.5})()
        line2 = type('Line', (), {'lam': 15.5})()
        lines = [(line1, 10.0, 0.0), (line2, 8.0, 0.0)]

        n, centers = fe._estimate_components_from_user_selection(
            x, y, xmin=14.0, xmax=16.0,
            lines_with_intensity=lines,
            line_threshold=0.03,
        )

        assert n == 2
        assert len(centers) == 2

    def test_weak_line_filtered_out(self):
        """A line below the threshold fraction should be filtered."""
        fe = FittingEngine()
        x = np.linspace(14, 16, 500)
        y = _gaussian(x, 10.0, 15.0, 0.05)

        line_strong = type('Line', (), {'lam': 15.0})()
        line_weak = type('Line', (), {'lam': 15.5})()
        # Second line has intensity < 3% of max (10.0)
        lines = [(line_strong, 10.0, 0.0), (line_weak, 0.01, 0.0)]

        n, centers = fe._estimate_components_from_user_selection(
            x, y, xmin=14.0, xmax=16.0,
            lines_with_intensity=lines,
            line_threshold=0.03,
        )

        assert n == 1
