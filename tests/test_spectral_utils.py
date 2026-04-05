# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.spectral_utils.

Tests the flux-conserving spectral resampling functions: make_bins and spectres.
"""

import pytest
import numpy as np

from iSLAT.Modules.DataProcessing.spectral_utils import make_bins, spectres


# ============================================================================
# make_bins
# ============================================================================

class TestMakeBins:
    """Tests for the make_bins helper."""

    def test_uniform_grid_edges(self):
        """Uniform spacing should give symmetric half-step edges."""
        wavs = np.array([10.0, 11.0, 12.0, 13.0])
        edges, widths = make_bins(wavs)

        assert len(edges) == len(wavs) + 1
        assert len(widths) == len(wavs)
        np.testing.assert_allclose(edges, [9.5, 10.5, 11.5, 12.5, 13.5])
        np.testing.assert_allclose(widths, [1.0, 1.0, 1.0, 1.0])

    def test_non_uniform_grid(self):
        """Non-uniform spacing: edges should be midpoints of neighbours."""
        wavs = np.array([1.0, 2.0, 4.0])
        edges, widths = make_bins(wavs)

        assert edges[0] == pytest.approx(0.5)          # 1.0 - (2.0-1.0)/2
        assert edges[1] == pytest.approx(1.5)          # (1.0+2.0)/2
        assert edges[2] == pytest.approx(3.0)          # (2.0+4.0)/2
        assert edges[3] == pytest.approx(5.0)          # 4.0 + (4.0-2.0)/2
        assert len(widths) == 3

    def test_two_element_grid(self):
        """Minimum valid grid with two wavelengths."""
        wavs = np.array([5.0, 7.0])
        edges, widths = make_bins(wavs)

        assert len(edges) == 3
        assert len(widths) == 2
        assert edges[0] == pytest.approx(4.0)
        assert edges[1] == pytest.approx(6.0)
        assert edges[2] == pytest.approx(8.0)

    def test_bin_widths_sum_to_total_range(self):
        """Sum of all bin widths should equal total bin range (edge-to-edge)."""
        wavs = np.linspace(5.0, 25.0, 100)
        edges, widths = make_bins(wavs)
        np.testing.assert_allclose(widths.sum(), edges[-1] - edges[0], rtol=1e-12)

    def test_edges_monotonically_increasing(self):
        """Bin edges must be strictly increasing for a sorted input grid."""
        wavs = np.linspace(4.5, 28.0, 500)
        edges, _ = make_bins(wavs)
        assert np.all(np.diff(edges) > 0)

    def test_widths_positive(self):
        """All bin widths must be positive."""
        wavs = np.array([1.0, 1.5, 3.0, 7.0, 8.0])
        _, widths = make_bins(wavs)
        assert np.all(widths > 0)


# ============================================================================
# spectres — identity / constant tests
# ============================================================================

class TestSpectresBasic:
    """Basic correctness tests for spectres."""

    def test_identity_resampling(self):
        """Resampling onto the same grid should return identical fluxes."""
        wavs = np.linspace(10, 20, 100)
        flux = np.sin(wavs)
        result = spectres(wavs, wavs, flux)
        np.testing.assert_allclose(result, flux, atol=1e-10)

    def test_constant_spectrum_preserved(self):
        """A flat spectrum should remain flat after resampling."""
        old = np.linspace(10, 20, 200)
        new = np.linspace(11, 19, 50)
        flat = np.full_like(old, 42.0)
        result = spectres(new, old, flat)
        np.testing.assert_allclose(result, 42.0, atol=1e-10)

    def test_linear_spectrum_preserved(self):
        """A linear spectrum should be well preserved."""
        old = np.linspace(10, 20, 500)
        new = np.linspace(11, 19, 80)
        linear = 3.0 * old + 7.0
        result = spectres(new, old, linear)
        expected = 3.0 * new + 7.0
        np.testing.assert_allclose(result, expected, rtol=1e-3)

    def test_coarser_resampling(self):
        """Resampling to a coarser grid should still conserve a constant."""
        old = np.linspace(5, 25, 1000)
        new = np.linspace(6, 24, 20)
        flux = np.full_like(old, 100.0)
        result = spectres(new, old, flux)
        np.testing.assert_allclose(result, 100.0, atol=1e-8)

    def test_finer_resampling(self):
        """Resampling to a finer grid should preserve a constant."""
        old = np.linspace(10, 20, 50)
        new = np.linspace(11, 19, 500)
        flux = np.full_like(old, 7.5)
        result = spectres(new, old, flux)
        np.testing.assert_allclose(result, 7.5, atol=1e-8)


# ============================================================================
# spectres — fill / out-of-range
# ============================================================================

class TestSpestresFill:
    """Tests for out-of-range behavior and fill values."""

    def test_fill_default_zero(self):
        """Out-of-range bins should default to 0.0."""
        old = np.linspace(10, 20, 100)
        new = np.array([5.0, 15.0, 25.0])
        flux = np.ones_like(old)
        result = spectres(new, old, flux)
        assert result[0] == 0.0
        assert result[2] == 0.0

    def test_fill_custom_value(self):
        """Custom fill value should appear for out-of-range bins."""
        old = np.linspace(10, 20, 100)
        new = np.array([5.0, 15.0, 25.0])
        flux = np.ones_like(old)
        result = spectres(new, old, flux, fill=-999.0)
        assert result[0] == -999.0
        assert result[2] == -999.0
        np.testing.assert_allclose(result[1], 1.0, atol=1e-10)

    def test_all_out_of_range(self):
        """If all new wavelengths are outside old grid, all filled."""
        old = np.linspace(10, 20, 100)
        new = np.array([1.0, 2.0, 3.0])
        flux = np.ones_like(old)
        result = spectres(new, old, flux, fill=-1.0)
        np.testing.assert_array_equal(result, -1.0)

    def test_fill_nan(self):
        """NaN fill value should produce NaN for out-of-range bins."""
        old = np.linspace(10, 20, 100)
        new = np.array([5.0, 15.0])
        flux = np.ones_like(old)
        result = spectres(new, old, flux, fill=np.nan)
        assert np.isnan(result[0])
        np.testing.assert_allclose(result[1], 1.0, atol=1e-10)

    def test_verbose_warning(self):
        """verbose=True should emit RuntimeWarning when fill is used."""
        old = np.linspace(10, 20, 100)
        new = np.array([1.0, 2.0])
        flux = np.ones_like(old)
        with pytest.warns(RuntimeWarning, match="spectres"):
            spectres(new, old, flux, verbose=True)

    def test_partial_overlap_left(self):
        """Only the left portion of new grid is out of range."""
        old = np.linspace(10, 20, 200)
        new = np.array([8.0, 9.0, 12.0, 15.0])
        flux = np.full_like(old, 5.0)
        result = spectres(new, old, flux, fill=0.0)
        assert result[0] == 0.0
        assert result[1] == 0.0
        np.testing.assert_allclose(result[2], 5.0, atol=1e-8)
        np.testing.assert_allclose(result[3], 5.0, atol=1e-8)


# ============================================================================
# spectres — flux conservation
# ============================================================================

class TestSpectresConservation:
    """Tests verifying flux-conservation properties."""

    def test_total_flux_conserved_gaussian(self):
        """Total integrated flux of a Gaussian should be roughly conserved."""
        old = np.linspace(10, 20, 1000)
        new = np.linspace(10.5, 19.5, 200)
        # Gaussian centered at 15 with sigma=1
        flux = np.exp(-0.5 * ((old - 15.0) / 1.0) ** 2)
        result = spectres(new, old, flux)

        # Integrated flux comparison (trapezoidal rule)
        old_integral = np.trapezoid(flux, old)
        # Only compare over the overlapping range
        old_mask = (old >= new[0]) & (old <= new[-1])
        old_trimmed_integral = np.trapezoid(flux[old_mask], old[old_mask])
        new_integral = np.trapezoid(result, new)
        np.testing.assert_allclose(new_integral, old_trimmed_integral, rtol=0.02)

    def test_delta_function_peak_preserved(self):
        """A narrow spike's peak region should still be non-zero after resampling."""
        old = np.linspace(10, 20, 500)
        flux = np.zeros_like(old)
        # Put a spike near 15.0
        spike_idx = np.argmin(np.abs(old - 15.0))
        flux[spike_idx] = 100.0
        new = np.linspace(10, 20, 200)
        result = spectres(new, old, flux)
        # The spike should appear somewhere in the result
        assert np.max(result) > 0

    def test_step_function(self):
        """A step function should produce a smooth transition when resampled."""
        old = np.linspace(10, 20, 1000)
        flux = np.where(old >= 15.0, 10.0, 0.0)
        new = np.linspace(10.5, 19.5, 100)
        result = spectres(new, old, flux)
        # Below 15 should be ~0, above should be ~10
        low_mask = new < 14.5
        high_mask = new > 15.5
        np.testing.assert_allclose(result[low_mask], 0.0, atol=1e-8)
        np.testing.assert_allclose(result[high_mask], 10.0, atol=0.5)


# ============================================================================
# spectres — return type and shape
# ============================================================================

class TestSpectresShape:
    """Shape and dtype tests for spectres output."""

    def test_output_length_matches_new_grid(self):
        """Output should always have the same length as new_wavs."""
        old = np.linspace(10, 20, 300)
        for n in [2, 5, 50, 300, 1000]:  # make_bins requires ≥ 2 points
            new = np.linspace(11, 19, n)
            flux = np.ones_like(old)
            result = spectres(new, old, flux)
            assert result.shape == (n,), f"Mismatch for n={n}"

    def test_output_dtype_float64(self):
        """Output should always be float64."""
        old = np.linspace(10, 20, 50)
        new = np.linspace(11, 19, 20)
        flux = np.ones_like(old, dtype=np.float32)
        result = spectres(new, old, flux)
        assert result.dtype == np.float64
