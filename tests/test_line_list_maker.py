# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.LineListMaker.

Tests the chainable builder API for filtering and exporting line lists.
"""

import json

import pytest
import numpy as np
import warnings
from pathlib import Path

import pandas as pd

from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker, _ensure_dataframe
from iSLAT.Modules.DataProcessing.LineFilterExpression import (
    Condition,
    ConditionError,
)


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


# ============================================================================
# filter_expression / filter_any / filter_all
# ============================================================================

class TestFilterExpression:
    """Boolean AND/OR filtering through the two-level expression model."""

    def test_user_motivating_example(self, maker):
        """Two wavelength ranges OR an E_up cut - the literal user request."""
        maker.filter_any(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=7.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=13.0),
            Condition('range', field='e_up', op='ge', min_val=5000.0))
        assert sorted(round(v, 3) for v in maker.df['lam']) == [6.5, 12.407]

    def test_dict_form_accepted(self, maker):
        maker.filter_expression({
            'join': 'AND',
            'groups': [{'join': 'OR', 'conditions': [
                {'kind': 'range', 'field': 'lam', 'op': 'between',
                 'min_val': 6.0, 'max_val': 7.0},
                {'kind': 'range', 'field': 'lam', 'op': 'between',
                 'min_val': 12.0, 'max_val': 13.0},
                {'kind': 'range', 'field': 'e_up', 'op': 'ge',
                 'min_val': 5000.0}]}]})
        assert sorted(round(v, 3) for v in maker.df['lam']) == [6.5, 12.407]

    def test_two_groups_and(self, maker):
        """(A OR B) AND e_up <= 5000."""
        maker.filter_expression({'join': 'AND', 'groups': [
            {'join': 'OR', 'conditions': [
                {'kind': 'range', 'field': 'lam', 'op': 'between',
                 'min_val': 6.0, 'max_val': 7.0},
                {'kind': 'range', 'field': 'lam', 'op': 'between',
                 'min_val': 12.0, 'max_val': 13.0}]},
            {'join': 'AND', 'conditions': [
                {'kind': 'range', 'field': 'e_up', 'op': 'le',
                 'max_val': 5000.0}]}]})
        assert sorted(round(v, 3) for v in maker.df['lam']) == [12.407]

    def test_filter_all_is_and(self, maker):
        maker.filter_all(
            Condition('range', field='lam', op='ge', min_val=10.0),
            Condition('range', field='e_up', op='le', max_val=5000.0))
        assert sorted(round(v, 3) for v in maker.df['lam']) == \
            [12.407, 14.95, 17.221, 22.5]

    def test_records_exactly_one_entry(self, maker):
        maker.filter_any(Condition('range', field='lam', op='ge', min_val=10.0))
        assert len(maker.filters) == 1
        assert maker.filters[0][0] == 'filter_expression'

    def test_recorded_kwargs_are_plain_json(self, maker):
        maker.filter_any(Condition('range', field='lam', op='ge', min_val=10.0))
        json.dumps(maker.filters[0][1])  # must not raise

    def test_narrows_current_frame(self, maker):
        """An expression composes as an implicit AND with earlier filters."""
        maker.filter_wavelength(min_val=13.0)
        maker.filter_any(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=7.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=15.0))
        assert sorted(round(v, 3) for v in maker.df['lam']) == [14.95]

    def test_returns_self(self, maker):
        result = maker.filter_any(
            Condition('range', field='lam', op='ge', min_val=1.0))
        assert result is maker

    def test_empty_result_is_legal(self, maker):
        maker.filter_all(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=7.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=13.0))
        assert len(maker) == 0

    def test_raise_records_nothing(self, maker):
        """A bad condition must leave the log and the frame untouched."""
        with pytest.raises(ConditionError):
            maker.filter_any(
                Condition('range', field='nope', op='ge', min_val=1.0))
        assert maker.filters == []
        assert len(maker) == 5

    def test_no_duplicate_rows_in_export(self, maker, tmp_path):
        maker.filter_any(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=15.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=23.0))
        out = maker.to_csv(tmp_path / 'or.csv')
        written = pd.read_csv(out)
        assert len(written) == len(maker)
        assert len(set(written['lam'])) == len(written)


class TestExpressionLogLifecycle:
    """An expression is one ordinary entry in the flat filter log."""

    def test_survives_replay(self, maker):
        maker.filter_wavelength(min_val=5.0)
        maker.filter_any(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=7.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=13.0))
        maker.remove_filter(0)
        assert [n for n, _ in maker.filters] == ['filter_expression']
        assert sorted(round(v, 3) for v in maker.df['lam']) == [6.5, 12.407]

    def test_pop_filter_removes_expression(self, maker):
        maker.filter_any(Condition('range', field='lam', op='ge', min_val=20.0))
        maker.pop_filter()
        assert len(maker) == 5
        assert maker.filters == []

    def test_reset_clears_expression(self, maker):
        maker.filter_any(Condition('range', field='lam', op='ge', min_val=20.0))
        maker.reset()
        assert maker.filters == []
        assert maker.expression is None

    def test_copy_is_independent(self, maker):
        clone = maker.copy()
        clone.filter_any(Condition('range', field='lam', op='ge', min_val=20.0))
        assert len(maker) == 5
        assert maker.filters == []
        assert len(clone) == 1

    def test_expression_property_returns_last(self, maker):
        maker.filter_any(Condition('range', field='lam', op='ge', min_val=10.0))
        expr = maker.expression
        assert expr is not None
        assert expr.groups[0].join == 'OR'
        assert expr.groups[0].conditions[0].field == 'lam'

    def test_summary_renders_expression(self, maker):
        maker.filter_any(
            Condition('range', field='lam', op='between',
                      min_val=6.0, max_val=7.0),
            Condition('range', field='lam', op='between',
                      min_val=12.0, max_val=13.0))
        summary = maker.summary()
        assert 'filter_expression' in summary
        assert 'OR' in summary
        assert "{'kind':" not in summary

    def test_original_df_is_the_unfiltered_snapshot(self, maker):
        maker.filter_wavelength(min_val=20.0)
        assert len(maker) == 1
        assert len(maker.original_df) == 5


class TestFilterSpeciesReplay:
    """filter_species records a kwarg but takes var-positional args."""

    def test_remove_filter_after_filter_species(self):
        df = pd.DataFrame({
            'species': ['H2O', 'H2O', 'CO'],
            'lam': [1.0, 2.0, 3.0],
        })
        maker = LineListMaker(df)
        maker.filter_range('lam', min_val=0)
        maker.filter_species('H2O')
        maker.remove_filter(0)  # must not raise TypeError
        assert [n for n, _ in maker.filters] == ['filter_species']
        assert len(maker) == 2


class TestReplayIsAtomic:
    """A failed replay must never leave an empty log over unfiltered data."""

    def test_failed_replay_preserves_filters_and_rows(self):
        base = pd.DataFrame([
            {'lam': 1.0, 'e_up': 100.0, 'flux': 5.0},
            {'lam': 2.0, 'e_up': 200.0, 'flux': 7.0},
        ])
        maker = LineListMaker(base, species='H2O')
        maker.filter_any(Condition('range', field='flux', op='ge', min_val=6.0))
        maker.filter_range('e_up', min_val=0.0)
        assert len(maker) == 1

        # The baseline loses the column the recorded expression needs.
        maker._df_original = maker._df_original.drop(columns=['flux'])
        with pytest.raises(ConditionError):
            maker.pop_filter()

        # Had this wiped, a caller swallowing the error would export
        # every line believing the list was filtered.
        assert [n for n, _ in maker.filters] == ['filter_expression', 'filter_range']
        assert len(maker) == 1

    def test_failed_remove_filter_preserves_state(self):
        base = pd.DataFrame([
            {'lam': 1.0, 'e_up': 100.0, 'flux': 5.0},
            {'lam': 2.0, 'e_up': 200.0, 'flux': 7.0},
        ])
        maker = LineListMaker(base, species='H2O')
        maker.filter_range('e_up', min_val=0.0)
        maker.filter_any(Condition('range', field='flux', op='ge', min_val=6.0))
        maker._df_original = maker._df_original.drop(columns=['flux'])
        with pytest.raises(ConditionError):
            maker.remove_filter(0)
        assert len(maker.filters) == 2
        assert len(maker) == 1

    def test_successful_removal_still_replays(self):
        base = pd.DataFrame([
            {'lam': 1.0, 'e_up': 100.0, 'flux': 5.0},
            {'lam': 2.0, 'e_up': 200.0, 'flux': 7.0},
        ])
        maker = LineListMaker(base, species='H2O')
        maker.filter_any(Condition('range', field='flux', op='ge', min_val=6.0))
        maker.filter_range('e_up', min_val=0.0)
        maker.remove_filter(1)
        assert [n for n, _ in maker.filters] == ['filter_expression']
        assert len(maker) == 1
