# -*- coding: utf-8 -*-
"""Unit tests for the Spectrum class."""

import pytest
import numpy as np
import pandas as pd


class TestSpectrum:
    """Tests for the Spectrum class."""

    def _make_spectrum(self, lam_min=10.0, lam_max=20.0, dlambda=0.01,
                       R=3000.0, distance=160.0):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        return Spectrum(
            lam_min=lam_min, lam_max=lam_max, dlambda=dlambda,
            R=R, distance=distance,
        )

    def test_init(self):
        spec = self._make_spectrum()
        assert spec._lam_min == 10.0
        assert spec._lam_max == 20.0
        assert spec._dlambda == 0.01
        assert spec._R == 3000.0
        assert spec._distance == 160.0

    def test_init_invalid_range(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        with pytest.raises(ValueError, match='lam_min must be < lam_max'):
            Spectrum(lam_min=20.0, lam_max=10.0, dlambda=0.01, R=3000, distance=160)

    def test_init_equal_range(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        with pytest.raises(ValueError, match='lam_min must be < lam_max'):
            Spectrum(lam_min=15.0, lam_max=15.0, dlambda=0.01, R=3000, distance=160)

    def test_wavelength_grid(self):
        spec = self._make_spectrum()
        grid = spec.lamgrid
        assert isinstance(grid, np.ndarray)
        assert grid[0] == 10.0
        assert grid[-1] == 20.0
        assert len(grid) > 0
        # Check monotonically increasing
        assert np.all(np.diff(grid) > 0)

    def test_flux_empty_spectrum(self):
        spec = self._make_spectrum()
        flux = spec.flux
        assert isinstance(flux, np.ndarray)
        assert len(flux) == len(spec.lamgrid)
        # No intensities added, flux should be zeros
        np.testing.assert_array_equal(flux, np.zeros_like(spec.lamgrid))

    def test_flux_jy_empty_spectrum(self):
        spec = self._make_spectrum()
        flux_jy = spec.flux_jy
        assert isinstance(flux_jy, np.ndarray)
        assert len(flux_jy) == len(spec.lamgrid)

    def test_reset(self):
        spec = self._make_spectrum()
        # Access flux to fill cache
        _ = spec.flux
        spec.reset()
        assert spec._flux is None
        assert spec._flux_jy is None
        assert spec._flux_valid is False
        assert len(spec._I_arrays) == 0
        assert len(spec._lam_arrays) == 0

    def test_cache_stats(self):
        spec = self._make_spectrum()
        stats = spec.get_cache_stats()
        assert 'hits' in stats
        assert 'misses' in stats
        assert 'invalidations' in stats

    def test_flux_caching(self):
        spec = self._make_spectrum()
        f1 = spec.flux
        f2 = spec.flux
        # Second access should be a cache hit
        np.testing.assert_array_equal(f1, f2)
        assert spec._flux_valid is True

    def test_wavelength_range_property(self):
        spec = self._make_spectrum(lam_min=5.0, lam_max=25.0)
        wr = spec.wavelength_range
        assert wr == (5.0, 25.0)

    def test_wavelength_range_setter(self):
        spec = self._make_spectrum()
        spec._flux_valid = True
        spec.wavelength_range = (8.0, 18.0)
        assert spec.wavelength_range == (8.0, 18.0)
        # Should have invalidated the cache
        assert spec._flux_valid is False

    def test_components_empty(self):
        spec = self._make_spectrum()
        assert spec.components == []

    def test_R_func_property(self):
        spec = self._make_spectrum()
        assert spec.R_func is None
        rfunc = lambda w: 3000 * np.ones_like(w)
        spec.R_func = rfunc
        assert spec.R_func is rfunc

    def test_R_func_invalidates_cache(self):
        spec = self._make_spectrum()
        _ = spec.flux
        spec.R_func = lambda w: 2500 * np.ones_like(w)
        assert spec._flux_valid is False

    def test_get_table(self):
        spec = self._make_spectrum()
        table = spec.get_table
        assert isinstance(table, pd.DataFrame)
        assert 'lam' in table.columns
        assert 'flux' in table.columns
        assert 'flux_jy' in table.columns
        assert len(table) == len(spec.lamgrid)

    def test_init_with_wavelength_range(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        spec = Spectrum(
            lam_min=10.0, lam_max=20.0, dlambda=0.01,
            R=3000, distance=160.0,
            wavelength_range=(12.0, 18.0),
        )
        assert spec.wavelength_range == (12.0, 18.0)

    def test_init_with_R_func(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        rfunc = lambda w: 2000 + 100 * w
        spec = Spectrum(
            lam_min=10.0, lam_max=20.0, dlambda=0.01,
            R=3000, distance=160.0, R_func=rfunc,
        )
        assert spec.R_func is rfunc

    def test_very_small_grid(self):
        """Spectrum with a very small wavelength range should still work."""
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        spec = Spectrum(lam_min=10.0, lam_max=10.1, dlambda=0.01,
                        R=3000, distance=160.0)
        assert len(spec.lamgrid) >= 2
