# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.LineAnalyzer.

Tests line detection, characterization, flux integration,
and rotation diagram calculations.
"""

import pytest
import numpy as np

import iSLAT.Constants as c
from iSLAT.Modules.DataProcessing.LineAnalyzer import LineAnalyzer


# ============================================================================
# Construction
# ============================================================================

class TestLineAnalyzerInit:
    """Tests for LineAnalyzer initialization."""

    def test_default_state(self):
        la = LineAnalyzer()
        assert la.detected_lines == []
        assert la.line_measurements == {}
        assert la.min_snr == 0.1
        assert la.min_line_width > 0
        assert la.max_line_width > la.min_line_width

    def test_class_level_settings(self):
        """Class-level parallel settings should be accessible."""
        assert isinstance(LineAnalyzer.PARALLEL_LINE_FITTING, bool)
        assert LineAnalyzer.LINE_FITTING_MAX_WORKERS is None or isinstance(
            LineAnalyzer.LINE_FITTING_MAX_WORKERS, int
        )


# ============================================================================
# find_single_lines
# ============================================================================

class TestFindSingleLines:
    """Tests for the find_single_lines method."""

    @staticmethod
    def _make_gaussian(center, sigma, amplitude, wavs):
        """Generate a Gaussian emission line profile."""
        return amplitude * np.exp(-0.5 * ((wavs - center) / sigma) ** 2)

    def test_single_isolated_emission(self):
        """A single strong isolated Gaussian should be detected."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 1000)
        flux = self._make_gaussian(15.0, 0.05, 10.0, wavs)

        lines = la.find_single_lines(wavs, flux, specsep=0.5, line_threshold=0.05)

        assert len(lines) >= 1
        # The detected line should be near 15.0
        detected_wavelengths = [l['wavelength'] for l in lines]
        assert any(abs(w - 15.0) < 0.1 for w in detected_wavelengths)

    def test_no_lines_in_flat_spectrum(self):
        """A perfectly flat spectrum should produce no detected lines."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 500)
        flux = np.ones_like(wavs) * 0.001  # flat, low-level

        lines = la.find_single_lines(wavs, flux, specsep=0.5, line_threshold=0.5)

        # Should find 0 (or very few spurious) lines
        assert len(lines) == 0

    def test_blended_lines_filtered_by_isolation(self):
        """Two closely spaced lines should fail the isolation test."""
        la = LineAnalyzer()
        wavs = np.linspace(14, 16, 2000)
        # Two Gaussians separated by 0.1 microns (within specsep)
        flux = (self._make_gaussian(15.0, 0.02, 10.0, wavs)
                + self._make_gaussian(15.05, 0.02, 8.0, wavs))

        lines = la.find_single_lines(
            wavs, flux, specsep=0.2, line_threshold=0.05, isolation_threshold=0.1
        )

        # Both should be rejected as non-isolated (they're within specsep of each other)
        # Exact count depends on thresholds, but should be fewer than 2
        detected_wavelengths = [l['wavelength'] for l in lines]
        # At most one (or zero) should survive isolation filtering
        assert len(detected_wavelengths) <= 1

    def test_widely_separated_lines_both_detected(self):
        """Two well-separated lines should both be detected."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 2000)
        flux = (self._make_gaussian(12.0, 0.05, 10.0, wavs)
                + self._make_gaussian(18.0, 0.05, 8.0, wavs))

        lines = la.find_single_lines(
            wavs, flux, specsep=1.0, line_threshold=0.05, isolation_threshold=0.1
        )

        assert len(lines) >= 2
        wavelengths = sorted([l['wavelength'] for l in lines])
        assert any(abs(w - 12.0) < 0.2 for w in wavelengths)
        assert any(abs(w - 18.0) < 0.2 for w in wavelengths)

    def test_results_sorted_by_wavelength(self):
        """Detected lines should be sorted by wavelength."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 2000)
        flux = (self._make_gaussian(18.0, 0.05, 10.0, wavs)
                + self._make_gaussian(12.0, 0.05, 8.0, wavs)
                + self._make_gaussian(15.0, 0.05, 6.0, wavs))

        lines = la.find_single_lines(
            wavs, flux, specsep=0.5, line_threshold=0.05, isolation_threshold=0.1
        )

        if len(lines) > 1:
            wavelengths = [l['wavelength'] for l in lines]
            assert wavelengths == sorted(wavelengths)

    def test_detected_lines_stored_as_attribute(self):
        """find_single_lines should also set self.detected_lines."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 1000)
        flux = self._make_gaussian(15.0, 0.05, 10.0, wavs)

        result = la.find_single_lines(wavs, flux, specsep=0.5, line_threshold=0.05)
        assert la.detected_lines is result

    def test_line_info_has_required_keys(self):
        """Each detected line dict should have wavelength, peak_flux, peak_index."""
        la = LineAnalyzer()
        wavs = np.linspace(10, 20, 1000)
        flux = self._make_gaussian(15.0, 0.05, 10.0, wavs)

        lines = la.find_single_lines(wavs, flux, specsep=0.5, line_threshold=0.05)

        if lines:
            line = lines[0]
            assert 'wavelength' in line
            assert 'peak_flux' in line
            assert 'peak_index' in line


# ============================================================================
# flux_integral (standalone function in spectral_utils)
# ============================================================================

from iSLAT.Modules.DataProcessing.spectral_utils import flux_integral


class TestFluxIntegral:
    """Tests for flux_integral from spectral_utils."""

    def test_zero_flux_returns_zero(self):
        """Zero flux should give zero integral."""
        lam = np.linspace(10, 20, 500)
        flux = np.zeros_like(lam)
        err = np.ones_like(lam) * 1e-18

        result, err_result = flux_integral(lam, flux, err, 12.0, 18.0)
        assert result == pytest.approx(0.0, abs=1e-30)

    def test_no_overlap_returns_zero(self):
        """If lam_min/lam_max don't overlap data, should return (0, 0)."""
        lam = np.linspace(10, 20, 500)
        flux = np.ones_like(lam)
        err = np.ones_like(lam) * 0.01

        result, err_result = flux_integral(lam, flux, err, 25.0, 30.0)
        assert result == 0.0
        assert err_result == 0.0

    def test_integral_sign(self):
        """Positive flux should produce a non-zero flux integral (converts to erg/s/cm^2)."""
        lam = np.linspace(10, 20, 1000)
        flux = np.ones_like(lam) * 1.0  # 1 Jy constant
        err = np.ones_like(lam) * 0.1

        result, err_result = flux_integral(lam, flux, err, 12.0, 18.0)
        assert result != 0.0

    def test_null_error_returns_zero_error(self):
        """When err=None, error result should be 0."""
        lam = np.linspace(10, 20, 500)
        flux = np.ones_like(lam)

        result, err_result = flux_integral(lam, flux, None, 12.0, 18.0)
        assert err_result == 0.0

    def test_single_point_returns_zero(self):
        """Range with only one data point should return (0, 0)."""
        lam = np.array([15.0])
        flux = np.array([1.0])
        err = np.array([0.1])

        result, err_result = flux_integral(lam, flux, err, 14.0, 16.0)
        assert result == 0.0
        assert err_result == 0.0

    def test_wider_range_larger_integral(self):
        """Wider wavelength range should give a larger integral magnitude."""
        lam = np.linspace(10, 20, 1000)
        flux = np.ones_like(lam) * 1.0
        err = np.ones_like(lam) * 0.1

        narrow_result, _ = flux_integral(lam, flux, err, 14.0, 16.0)
        wide_result, _ = flux_integral(lam, flux, err, 12.0, 18.0)
        assert abs(wide_result) > abs(narrow_result)


# ============================================================================
# add_rotation_diagram_values
# ============================================================================

class TestRotationDiagram:
    """Tests for rotation diagram calculations."""

    def test_valid_entry_gets_rd_y(self):
        """Entry with valid a_stein, g_up, Flux_fit should get a finite RD_y."""
        la = LineAnalyzer()
        entries = [{
            'a_stein': 1.05e-2,
            'g_up': 21,
            'Flux_fit': 1e-15,
            'lam': 12.407,
        }]
        la.add_rotation_diagram_values(entries)
        assert 'RD_y' in entries[0]
        assert np.isfinite(entries[0]['RD_y'])

    def test_zero_a_stein_gives_nan(self):
        """Zero Einstein-A should give NaN RD_y."""
        la = LineAnalyzer()
        entries = [{'a_stein': 0, 'g_up': 21, 'Flux_fit': 1e-15, 'lam': 12.0}]
        la.add_rotation_diagram_values(entries)
        assert np.isnan(entries[0]['RD_y'])

    def test_zero_g_up_gives_nan(self):
        """Zero degeneracy should give NaN RD_y."""
        la = LineAnalyzer()
        entries = [{'a_stein': 0.01, 'g_up': 0, 'Flux_fit': 1e-15, 'lam': 12.0}]
        la.add_rotation_diagram_values(entries)
        assert np.isnan(entries[0]['RD_y'])

    def test_zero_flux_gives_nan(self):
        """Zero flux should give NaN RD_y."""
        la = LineAnalyzer()
        entries = [{'a_stein': 0.01, 'g_up': 21, 'Flux_fit': 0, 'lam': 12.0}]
        la.add_rotation_diagram_values(entries)
        assert np.isnan(entries[0]['RD_y'])

    def test_missing_keys_gives_nan(self):
        """Entry missing required keys should get NaN."""
        la = LineAnalyzer()
        entries = [{'lam': 12.0}]
        la.add_rotation_diagram_values(entries)
        assert np.isnan(entries[0]['RD_y'])

    def test_multiple_entries(self):
        """Multiple entries should each get their own RD_y."""
        la = LineAnalyzer()
        entries = [
            {'a_stein': 0.01, 'g_up': 21, 'Flux_fit': 1e-15, 'lam': 12.0},
            {'a_stein': 0.02, 'g_up': 13, 'Flux_fit': 2e-15, 'lam': 16.0},
            {'a_stein': 0, 'g_up': 5, 'Flux_fit': 1e-16, 'lam': 8.0},
        ]
        la.add_rotation_diagram_values(entries)
        assert np.isfinite(entries[0]['RD_y'])
        assert np.isfinite(entries[1]['RD_y'])
        assert np.isnan(entries[2]['RD_y'])

    def test_rd_y_is_rounded(self):
        """RD_y should be rounded to 3 decimal places."""
        la = LineAnalyzer()
        entries = [{'a_stein': 0.01, 'g_up': 21, 'Flux_fit': 1e-15, 'lam': 12.0}]
        la.add_rotation_diagram_values(entries)
        rd_y = entries[0]['RD_y']
        # Check it's rounded to 3 decimals
        assert rd_y == pytest.approx(round(rd_y, 3), abs=1e-10)
