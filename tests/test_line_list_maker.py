# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.LineListMaker.

Tests the chainable builder API for filtering and exporting line lists.
"""

import pytest
import numpy as np
import warnings
from pathlib import Path

import pandas as pd

from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker, _ensure_dataframe


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_lines_data():
    """Five sample lines spanning 6.5 – 22.5 µm."""
    return [
        {'nr': 1, 'lev_up': '0|10', 'lev_low': '0|9',
         'lam': 12.407, 'freq': 2.416e13, 'a_stein': 1.05e-2,
         'e_up': 4586.4, 'e_low': 3379.0, 'g_up': 21, 'g_low': 19},
        {'nr': 2, 'lev_up': '0|8', 'lev_low': '0|7',
         'lam': 14.950, 'freq': 2.005e13, 'a_stein': 2.10e-2,
         'e_up': 3200.0, 'e_low': 2500.0, 'g_up': 17, 'g_low': 15},
        {'nr': 3, 'lev_up': '0|6', 'lev_low': '0|5',
         'lam': 17.221, 'freq': 1.741e13, 'a_stein': 5.00e-3,
         'e_up': 2100.0, 'e_low': 1600.0, 'g_up': 13, 'g_low': 11},
        {'nr': 4, 'lev_up': '0|20', 'lev_low': '0|19',
         'lam': 6.500, 'freq': 4.612e13, 'a_stein': 3.00e-2,
         'e_up': 8000.0, 'e_low': 6500.0, 'g_up': 41, 'g_low': 39},
        {'nr': 5, 'lev_up': '0|15', 'lev_low': '0|14',
         'lam': 22.500, 'freq': 1.332e13, 'a_stein': 8.00e-3,
         'e_up': 1200.0, 'e_low': 800.0, 'g_up': 31, 'g_low': 29},
    ]


@pytest.fixture
def mll(sample_lines_data):
    """A MoleculeLineList for H2O."""
    return MoleculeLineList(molecule_id='H2O', lines_data=sample_lines_data)


@pytest.fixture
def maker(mll):
    """A LineListMaker initialized from the H2O MoleculeLineList."""
    return LineListMaker(mll)


# ============================================================================
# _ensure_dataframe helper
# ============================================================================

class TestEnsureDataframe:
    """Tests for the _ensure_dataframe helper."""

    def test_from_dataframe(self):
        df = pd.DataFrame({'lam': [10, 15, 20], 'a_stein': [0.01, 0.02, 0.03]})
        result = _ensure_dataframe(df)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        # Should be a copy
        assert result is not df

    def test_from_molecule_line_list(self, mll):
        result = _ensure_dataframe(mll, molecule_id='H2O')
        assert isinstance(result, pd.DataFrame)
        assert 'lam' in result.columns
        assert 'species' in result.columns

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError, match="Expected pd.DataFrame"):
            _ensure_dataframe([1, 2, 3])


# ============================================================================
# Construction
# ============================================================================

class TestLineListMakerInit:
    """Tests for LineListMaker construction."""

    def test_from_mll(self, maker):
        assert len(maker) == 5
        assert maker.species == 'H2O'

    def test_from_dataframe(self):
        df = pd.DataFrame({
            'species': ['CO'] * 3,
            'lam': [4.6, 4.7, 4.8],
            'a_stein': [0.01, 0.02, 0.03],
        })
        m = LineListMaker(df, species='CO')
        assert len(m) == 3
        assert m.species == 'CO'

    def test_repr(self, maker):
        r = repr(maker)
        assert 'H2O' in r
        assert '5' in r

    def test_len(self, maker):
        assert len(maker) == 5

    def test_summary(self, maker):
        s = maker.summary()
        assert 'H2O' in s
        assert 'No filters applied' in s

    def test_no_filters_initially(self, maker):
        assert maker.filters == []


# ============================================================================
# Wavelength filter
# ============================================================================

class TestFilterWavelength:
    """Tests for filter_wavelength."""

    def test_filter_min(self, maker):
        maker.filter_wavelength(min_val=13.0)
        assert all(maker.df['lam'] >= 13.0)
        assert len(maker) < 5

    def test_filter_max(self, maker):
        maker.filter_wavelength(max_val=15.0)
        assert all(maker.df['lam'] <= 15.0)

    def test_filter_range(self, maker):
        maker.filter_wavelength(min_val=10.0, max_val=18.0)
        df = maker.df
        assert all((df['lam'] >= 10.0) & (df['lam'] <= 18.0))
        # Lines at 6.5 and 22.5 should be excluded
        assert 6.5 not in df['lam'].values
        assert 22.5 not in df['lam'].values

    def test_chaining(self, maker):
        result = maker.filter_wavelength(min_val=10.0).filter_wavelength(max_val=18.0)
        assert result is maker  # Returns self
        assert len(maker.filters) == 2


# ============================================================================
# Energy / Einstein-A filters
# ============================================================================

class TestEnergyFilters:
    """Tests for filter_eup, filter_elow, filter_astein."""

    def test_filter_eup_max(self, maker):
        maker.filter_eup(max_val=4000)
        assert all(maker.df['e_up'] <= 4000)

    def test_filter_elow_min(self, maker):
        maker.filter_elow(min_val=2000)
        assert all(maker.df['e_low'] >= 2000)

    def test_filter_astein_min(self, maker):
        maker.filter_astein(min_val=1e-2)
        assert all(maker.df['a_stein'] >= 1e-2)

    def test_filter_gup(self, maker):
        maker.filter_gup(min_val=20)
        assert all(maker.df['g_up'] >= 20)

    def test_filter_glow(self, maker):
        maker.filter_glow(max_val=20)
        assert all(maker.df['g_low'] <= 20)


# ============================================================================
# Quantum filter
# ============================================================================

class TestFilterQuantum:
    """Tests for filter_quantum."""

    def test_exact_match(self, maker):
        maker.filter_quantum(lev_up='0|10')
        assert len(maker) == 1
        assert maker.df['lev_up'].iloc[0] == '0|10'

    def test_substring_match(self, maker):
        maker.filter_quantum(lev_up='0|', contains=True)
        # All lines have '0|' prefix
        assert len(maker) == 5

    def test_no_match(self, maker):
        maker.filter_quantum(lev_up='NONEXISTENT')
        assert len(maker) == 0


# ============================================================================
# Custom filter
# ============================================================================

class TestFilterCustom:
    """Tests for filter_custom."""

    def test_lambda_filter(self, maker):
        maker.filter_custom(
            lambda df: df['a_stein'] > df['a_stein'].median(),
            label='above_median_astein',
        )
        assert len(maker) < 5
        assert all(maker.df['a_stein'] > 0)

    def test_filter_recorded(self, maker):
        maker.filter_custom(lambda df: df['lam'] > 10, label='test')
        assert any(name == 'filter_custom' for name, _ in maker.filters)


# ============================================================================
# filter_range (generic)
# ============================================================================

class TestFilterRange:
    """Tests for the generic filter_range method."""

    def test_range_on_existing_column(self, maker):
        maker.filter_range('e_up', min_val=2000, max_val=5000)
        df = maker.df
        assert all((df['e_up'] >= 2000) & (df['e_up'] <= 5000))

    def test_range_on_missing_column_warns(self, maker):
        with pytest.warns(UserWarning, match="not in DataFrame"):
            maker.filter_range('nonexistent_col', min_val=0)
        # Length unchanged
        assert len(maker) == 5


# ============================================================================
# Sort
# ============================================================================

class TestSort:
    """Tests for the sort method."""

    def test_sort_by_lam(self, maker):
        maker.sort(by='lam')
        lams = maker.df['lam'].tolist()
        assert lams == sorted(lams)

    def test_sort_descending(self, maker):
        maker.sort(by='a_stein', ascending=False)
        vals = maker.df['a_stein'].tolist()
        assert vals == sorted(vals, reverse=True)

    def test_sort_returns_self(self, maker):
        result = maker.sort()
        assert result is maker


# ============================================================================
# Reset / pop_filter / remove_filter
# ============================================================================

class TestFilterManagement:
    """Tests for filter management: reset, pop, remove."""

    def test_reset_restores_all_lines(self, maker):
        maker.filter_wavelength(min_val=13.0)
        assert len(maker) < 5
        maker.reset()
        assert len(maker) == 5
        assert maker.filters == []

    def test_pop_filter(self, maker):
        maker.filter_wavelength(min_val=13.0)
        maker.filter_eup(max_val=4000)
        assert len(maker.filters) == 2
        maker.pop_filter()
        assert len(maker.filters) == 1

    def test_pop_empty_raises(self, maker):
        with pytest.raises(IndexError):
            maker.pop_filter()

    def test_remove_filter_by_index(self, maker):
        maker.filter_wavelength(min_val=10.0)
        maker.filter_eup(max_val=5000)
        maker.remove_filter(0)
        # First filter removed, second should remain
        assert len(maker.filters) == 1
        assert maker.filters[0][0] == 'filter_eup'


# ============================================================================
# Species
# ============================================================================

class TestSpecies:
    """Tests for species property."""

    def test_species_getter(self, maker):
        assert maker.species == 'H2O'

    def test_species_setter(self, maker):
        maker.species = 'CO'
        assert maker.species == 'CO'
        assert all(maker.df['species'] == 'CO')

    def test_filter_species(self):
        df = pd.DataFrame({
            'species': ['H2O', 'CO', 'H2O', 'OH'],
            'lam': [10, 12, 14, 16],
        })
        m = LineListMaker(df)
        m.filter_species('H2O')
        assert len(m) == 2
        assert all(m.df['species'] == 'H2O')


# ============================================================================
# Export — DataFrame
# ============================================================================

class TestExportDataFrame:
    """Tests for to_dataframe."""

    def test_returns_copy(self, maker):
        df = maker.to_dataframe()
        assert isinstance(df, pd.DataFrame)
        # Modifying copy shouldn't affect maker
        df.drop(df.index, inplace=True)
        assert len(maker) == 5

    def test_df_property(self, maker):
        df = maker.df
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 5


# ============================================================================
# Export — CSV
# ============================================================================

class TestExportCSV:
    """Tests for to_csv."""

    def test_basic_csv(self, maker, tmp_path):
        out = tmp_path / "test.csv"
        result_path = maker.to_csv(out)
        assert result_path.exists()

        df = pd.read_csv(result_path)
        assert 'lam' in df.columns
        assert 'species' in df.columns

    def test_extended_csv(self, maker, tmp_path):
        out = tmp_path / "test_ext.csv"
        maker.to_csv(out, extended=True)
        df = pd.read_csv(out)
        assert 'xmin' in df.columns
        assert 'xmax' in df.columns

    def test_extra_columns(self, maker, tmp_path):
        out = tmp_path / "test_extra.csv"
        maker.to_csv(out, extra_columns={'note': 'test_value'})
        df = pd.read_csv(out)
        assert 'note' in df.columns
        assert all(df['note'] == 'test_value')

    def test_csv_sorted_by_lam(self, maker, tmp_path):
        out = tmp_path / "test_sorted.csv"
        maker.to_csv(out)
        df = pd.read_csv(out)
        lams = df['lam'].tolist()
        assert lams == sorted(lams)


# ============================================================================
# Export — to_linelist
# ============================================================================

class TestExportLineList:
    """Tests for to_linelist."""

    def test_roundtrip(self, maker):
        ll = maker.to_linelist()
        assert isinstance(ll, MoleculeLineList)
        assert ll.num_lines == 5
        assert ll.molecule_id == 'H2O'

    def test_filtered_linelist(self, maker):
        maker.filter_wavelength(min_val=10.0, max_val=18.0)
        ll = maker.to_linelist()
        assert ll.num_lines == len(maker)


# ============================================================================
# Copy
# ============================================================================

class TestCopy:
    """Tests for the copy method."""

    def test_independent_copy(self, maker):
        copy = maker.copy()
        copy.filter_wavelength(min_val=20.0)
        # Original unaffected
        assert len(maker) == 5
        assert len(copy) < 5


# ============================================================================
# Append / Merge
# ============================================================================

class TestAppendMerge:
    """Tests for append and merge."""

    def test_append(self, maker):
        extra_df = pd.DataFrame({
            'species': ['CO'],
            'lam': [4.6],
            'a_stein': [0.5],
        })
        maker.append(extra_df, species='CO')
        assert len(maker) == 6

    def test_merge_two_makers(self, mll):
        m1 = LineListMaker(mll)
        m2 = LineListMaker(mll)
        merged = LineListMaker.merge(m1, m2)
        assert len(merged) == 10

    def test_merge_with_species_override(self, mll):
        m1 = LineListMaker(mll)
        m2 = LineListMaker(mll)
        merged = LineListMaker.merge(m1, m2, species_override='Combined')
        assert all(merged.df['species'] == 'Combined')


# ============================================================================
# Chain integration
# ============================================================================

class TestChainIntegration:
    """Integration tests for chained operations."""

    def test_full_chain(self, maker, tmp_path):
        """Typical workflow: filter → sort → export."""
        out = tmp_path / "chain_test.csv"
        (maker
            .filter_wavelength(min_val=10.0, max_val=20.0)
            .filter_eup(max_val=5000)
            .filter_astein(min_val=5e-3)
            .sort(by='lam')
            .to_csv(out))

        df = pd.read_csv(out)
        assert len(df) > 0
        assert all(df['lam'] >= 10.0)
        assert all(df['lam'] <= 20.0)
        assert all(df['e_up'] <= 5000)
        assert all(df['a_stein'] >= 5e-3)

    def test_filter_summary_records_all(self, maker):
        (maker
            .filter_wavelength(min_val=10.0)
            .filter_eup(max_val=5000)
            .filter_astein(min_val=1e-2))
        assert len(maker.filters) == 3
        summary = maker.summary()
        assert 'filter_wavelength' in summary
        assert 'filter_eup' in summary
        assert 'filter_astein' in summary
