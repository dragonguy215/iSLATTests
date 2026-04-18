# -*- coding: utf-8 -*-
"""Unit tests for the Intensity class."""

import pytest
import numpy as np

import iSLAT.Constants as c


class TestIntensity:
    """Tests for the Intensity class."""

    def _make_mock_line_list(self, n_lines=10):
        """Create a mock MoleculeLineList for Intensity tests."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        # Generate realistic line data
        np.random.seed(42)
        lines_data = []
        for i in range(n_lines):
            lam = 10.0 + i * 1.0  # wavelengths from 10 to 10+n
            freq = c.SPEED_OF_LIGHT_MICRONS / lam
            lines_data.append({
                'nr': i + 1,
                'lev_up': f'0|{2*i+2}',
                'lev_low': f'0|{2*i+1}',
                'lam': lam,
                'freq': freq,
                'a_stein': np.random.uniform(1e-3, 5e-2),
                'e_up': 1000.0 + i * 500,
                'e_low': 500.0 + i * 500,
                'g_up': 2 * i + 3,
                'g_low': 2 * i + 1,
            })
        mll = MoleculeLineList(molecule_id='TestMol', lines_data=lines_data)
        # Set up a dummy partition function
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015  # H2O
        return mll

    def test_init(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        assert inten.molecule is mll
        assert inten.intensity is None
        assert inten.tau is None
        assert inten.t_kin is None
        assert inten.n_mol is None
        assert inten.dv is None

    def test_init_with_wavelength_range(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll, wavelength_range=(12.0, 18.0))
        assert inten.wavelength_range == (12.0, 18.0)

    def test_blackbody_scalar(self):
        """Test the static blackbody function with scalar temperature."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        nu = np.array([1e13, 2e13, 3e13])
        T = np.array([500.0])
        bb = Intensity._bb(nu, T)
        assert isinstance(bb, np.ndarray)
        assert bb.shape == (3,)
        # Blackbody should be positive
        assert np.all(bb > 0)

    def test_blackbody_batch(self):
        """Test the blackbody function with multiple temperatures."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        nu = np.array([1e13, 2e13, 3e13])
        T = np.array([300.0, 500.0, 1000.0])
        bb = Intensity._bb(nu, T)
        assert bb.shape == (3, 3)  # (n_temps, n_freqs)
        assert np.all(bb > 0)

    def test_blackbody_increases_with_temperature(self):
        """Blackbody intensity should increase with temperature."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        nu = np.array([2e13])
        bb_low = Intensity._bb(nu, np.array([300.0]))
        bb_high = Intensity._bb(nu, np.array([1000.0]))
        assert bb_high > bb_low

    def test_blackbody_rayleigh_jeans(self):
        """In the Rayleigh-Jeans regime (low freq / high T), B ~ T * nu^2."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        # Very low frequency for given T ensures x = h*nu/(k*T) << 1
        nu = np.array([1e8])  # Very low frequency
        T = np.array([5000.0])
        bb = Intensity._bb(nu, T)
        assert bb.shape == (1,)
        assert bb[0] > 0

    def test_blackbody_very_high_temperature(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        nu = np.array([2e13])
        T = np.array([1e6])
        bb = Intensity._bb(nu, T)
        assert np.all(np.isfinite(bb))
        assert bb[0] > 0

    def test_blackbody_very_low_temperature(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        nu = np.array([2e13])
        T = np.array([10.0])
        bb = Intensity._bb(nu, T)
        assert np.all(np.isfinite(bb))
        assert bb[0] > 0

    def test_class_constants(self):
        """Verify pre-computed class-level constants exist and are reasonable."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        assert Intensity._SQRT_LN2_INV > 0
        assert Intensity._C3_OVER_8PI > 0
        assert Intensity._TWO_SQRT_LN2 > 0
        assert Intensity._TAU_THIN > 0

    def test_calc_intensity_requires_all_params(self):
        """calc_intensity should raise if any parameter is None."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        with pytest.raises(ValueError, match="t_kin, n_mol, and dv must all be provided"):
            inten.calc_intensity(t_kin=500, n_mol=None, dv=1.0)
        with pytest.raises(ValueError, match="t_kin, n_mol, and dv must all be provided"):
            inten.calc_intensity(t_kin=None, n_mol=1e18, dv=1.0)
        with pytest.raises(ValueError, match="t_kin, n_mol, and dv must all be provided"):
            inten.calc_intensity(t_kin=500, n_mol=1e18, dv=None)

    def test_calc_intensity_basic(self):
        """Basic intensity calculation should produce valid arrays."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        assert inten.intensity is not None
        assert inten.tau is not None
        assert inten.t_kin == 500.0
        assert inten.n_mol == 1e18
        assert inten.dv == 1.0
        assert len(inten.intensity) == 10
        assert len(inten.tau) == 10
        # Intensities should be non-negative
        assert np.all(inten.intensity >= 0)
        # Tau should be non-negative
        assert np.all(inten.tau >= 0)

    def test_calc_intensity_caching(self):
        """Second call with same params should use cache."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        i1 = inten.intensity.copy()
        # Second call with same params — should be cached
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        np.testing.assert_array_equal(inten.intensity, i1)

    def test_calc_intensity_different_params(self):
        """Different parameters should produce different results."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        i1 = inten.intensity.copy()
        inten.calc_intensity(t_kin=1000.0, n_mol=1e18, dv=1.0)
        i2 = inten.intensity.copy()
        assert not np.array_equal(i1, i2)

    def test_invalidate_cache(self):
        """invalidate_cache should clear stored intensity and tau."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        assert inten.intensity is not None
        inten.invalidate_cache()
        assert inten.intensity is None
        assert inten.tau is None
        assert inten._cache_valid is False

    def test_wavelength_range_setter_invalidates_cache(self):
        """Changing wavelength range should invalidate all caches."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        inten.wavelength_range = (12.0, 16.0)
        assert inten._cache_valid is False
        assert inten._sorted_idx is None

    def test_calc_intensity_radex_method(self):
        """RADEX method should also produce valid results."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0, method='radex')
        assert inten.intensity is not None
        assert np.all(inten.intensity >= 0)

    def test_calc_intensity_batch(self):
        """Batch intensity calculation should return correct shape."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        t_arr = np.array([300, 500, 700])
        n_arr = np.array([1e17, 1e18, 1e19])
        dv_arr = np.array([1.0, 1.0, 1.0])
        result = inten.calc_intensity_batch(t_arr, n_arr, dv_arr)
        assert result.shape[0] == 3   # 3 conditions
        assert result.shape[1] == 10  # 10 lines
        assert np.all(result >= 0)

    def test_repr(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_mock_line_list()
        inten = Intensity(mll)
        r = repr(inten)
        assert 'Intensity' in r
        assert 'TestMol' in r

    def test_gauss_quad_initialization(self):
        """Gaussian quadrature should initialize correctly."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        Intensity._initialize_gauss_quad()
        assert Intensity._GAUSS_QUAD_INITIALIZED is True
        assert Intensity._GAUSS_QUAD_X is not None
        assert Intensity._GAUSS_QUAD_W is not None
        assert len(Intensity._GAUSS_QUAD_X) == 20
        assert len(Intensity._GAUSS_QUAD_W) == 20
