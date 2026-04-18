# -*- coding: utf-8 -*-
"""Unit tests for MoleculeLineList."""

import pytest
import numpy as np


class TestMoleculeLineList:
    """Tests for MoleculeLineList."""

    def test_init_from_data(self, sample_lines_data):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        mll = MoleculeLineList(molecule_id='H2O', lines_data=sample_lines_data)
        assert mll.molecule_id == 'H2O'
        assert mll.num_lines == 5
        assert len(mll) == 5

    def test_init_empty(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        mll = MoleculeLineList(molecule_id='empty')
        assert mll.molecule_id == 'empty'
        assert mll.num_lines == 0

    def test_init_empty_lines_data(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        mll = MoleculeLineList(molecule_id='empty', lines_data=[])
        assert mll.num_lines == 0

    def test_get_wavelengths(self, molecule_line_list):
        wavelengths = molecule_line_list.get_wavelengths()
        assert isinstance(wavelengths, np.ndarray)
        assert len(wavelengths) == 5
        np.testing.assert_allclose(sorted(wavelengths), [6.5, 12.407, 14.95, 17.221, 22.5])

    def test_get_frequencies(self, molecule_line_list):
        freqs = molecule_line_list.get_frequencies()
        assert isinstance(freqs, np.ndarray)
        assert len(freqs) == 5
        assert all(f > 0 for f in freqs)

    def test_get_einstein_coefficients(self, molecule_line_list):
        a_stein = molecule_line_list.get_einstein_coefficients()
        assert isinstance(a_stein, np.ndarray)
        assert len(a_stein) == 5
        assert all(a > 0 for a in a_stein)

    def test_get_upper_energies(self, molecule_line_list):
        e_up = molecule_line_list.get_upper_energies()
        assert isinstance(e_up, np.ndarray)
        assert len(e_up) == 5

    def test_get_lower_energies(self, molecule_line_list):
        e_low = molecule_line_list.get_lower_energies()
        assert isinstance(e_low, np.ndarray)
        assert len(e_low) == 5

    def test_get_upper_weights(self, molecule_line_list):
        g_up = molecule_line_list.get_upper_weights()
        assert isinstance(g_up, np.ndarray)
        assert len(g_up) == 5

    def test_get_lower_weights(self, molecule_line_list):
        g_low = molecule_line_list.get_lower_weights()
        assert isinstance(g_low, np.ndarray)
        assert len(g_low) == 5

    def test_wavelength_range_filtering(self, molecule_line_list_with_range):
        """Lines outside range (10, 20) should be excluded."""
        mll = molecule_line_list_with_range
        n = mll.num_lines
        # Lines at 12.407, 14.95, 17.221 are in [10, 20]; 6.5 and 22.5 are out
        assert n == 3

        wavelengths = mll.get_wavelengths()
        assert len(wavelengths) == 3
        assert all(10.0 <= w <= 20.0 for w in wavelengths)

    def test_wavelength_range_change_invalidates_cache(self, molecule_line_list):
        mll = molecule_line_list
        _ = mll.get_wavelengths()  # populate cache
        assert mll.num_lines == 5

        mll.wavelength_range = (10.0, 20.0)
        assert mll.num_lines == 3

        mll.wavelength_range = None
        assert mll.num_lines == 5

    def test_get_lines_in_range(self, molecule_line_list):
        lines = molecule_line_list.get_lines_in_range(10.0, 15.0)
        assert len(lines) == 2  # 12.407 and 14.95

    def test_get_lines_in_range_empty(self, molecule_line_list):
        lines = molecule_line_list.get_lines_in_range(30.0, 40.0)
        assert len(lines) == 0

    def test_lines_as_namedtuple(self, molecule_line_list):
        nt = molecule_line_list.lines_as_namedtuple
        assert hasattr(nt, 'lam')
        assert hasattr(nt, 'freq')
        assert hasattr(nt, 'a_stein')
        assert hasattr(nt, 'e_up')
        assert hasattr(nt, 'e_low')
        assert hasattr(nt, 'g_up')
        assert hasattr(nt, 'g_low')
        assert hasattr(nt, 'nr')
        assert hasattr(nt, 'lev_up')
        assert hasattr(nt, 'lev_low')
        assert len(nt.lam) == 5

    def test_name_property(self, molecule_line_list):
        assert molecule_line_list.name == 'H2O'

    def test_fname_property(self, molecule_line_list):
        assert molecule_line_list.fname is None
        molecule_line_list.fname = 'test.par'
        assert molecule_line_list.fname == 'test.par'

    def test_get_ndarray(self, molecule_line_list):
        arr = molecule_line_list.get_ndarray()
        assert isinstance(arr, np.ndarray)
        assert arr.shape[0] == 5  # 5 lines
        assert arr.shape[1] == 10  # 10 fields per line

    def test_get_pandas_table(self, molecule_line_list):
        import pandas as pd
        df = molecule_line_list.get_pandas_table()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 5
        assert 'lam' in df.columns
        assert 'freq' in df.columns

    def test_get_ndarray_of_attribute(self, molecule_line_list):
        lam = molecule_line_list.get_ndarray_of_attribute('lam')
        assert isinstance(lam, np.ndarray)
        assert len(lam) == 5
        np.testing.assert_allclose(sorted(lam), [6.5, 12.407, 14.95, 17.221, 22.5])

    def test_caching_behavior(self, molecule_line_list):
        """Second call should return cached result."""
        w1 = molecule_line_list.get_wavelengths()
        w2 = molecule_line_list.get_wavelengths()
        # Should be the same array object due to caching
        np.testing.assert_array_equal(w1, w2)

    def test_invalidate_caches(self, molecule_line_list):
        """After invalidation, cache should be rebuilt."""
        _ = molecule_line_list.get_wavelengths()
        molecule_line_list._invalidate_caches()
        assert molecule_line_list._wavelengths_cache is None
        # Should rebuild on next access
        w = molecule_line_list.get_wavelengths()
        assert len(w) == 5

    def test_single_line(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        mll = MoleculeLineList(molecule_id='single', lines_data=[{
            'nr': 1, 'lev_up': '0|2', 'lev_low': '0|1',
            'lam': 15.0, 'freq': 2e13, 'a_stein': 0.01,
            'e_up': 2000, 'e_low': 1000, 'g_up': 5, 'g_low': 3,
        }])
        assert mll.num_lines == 1
        assert mll.get_wavelengths()[0] == 15.0

    def test_duplicate_wavelengths(self):
        """MoleculeLineList should handle duplicate wavelengths."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        lines_data = [
            {'nr': 1, 'lev_up': '0|2', 'lev_low': '0|1',
             'lam': 15.0, 'freq': 2e13, 'a_stein': 0.01,
             'e_up': 2000, 'e_low': 1000, 'g_up': 5, 'g_low': 3},
            {'nr': 2, 'lev_up': '0|4', 'lev_low': '0|3',
             'lam': 15.0, 'freq': 2e13, 'a_stein': 0.02,
             'e_up': 3000, 'e_low': 2000, 'g_up': 9, 'g_low': 7},
        ]
        mll = MoleculeLineList(molecule_id='dup', lines_data=lines_data)
        assert mll.num_lines == 2
        wavelengths = mll.get_wavelengths()
        assert np.all(wavelengths == 15.0)
