# -*- coding: utf-8 -*-
"""Unit tests for MoleculeLine and LineDataView."""

import pytest
import numpy as np


class TestMoleculeLine:
    """Tests for MoleculeLine and LineDataView."""

    def test_init_basic(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('H2O', sample_line_data)
        assert line.molecule_id == 'H2O'
        assert line.nr == 1
        assert line.lam == 12.407
        assert line.freq == 2.416e13
        assert line.a_stein == 1.05e-2
        assert line.e_up == 4586.4
        assert line.e_low == 3379.0
        assert line.g_up == 21
        assert line.g_low == 19
        assert line.lev_up == '0_0_0|10_2_9'
        assert line.lev_low == '0_0_0|9_1_8'

    def test_init_with_missing_fields(self):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('CO', {})
        assert line.molecule_id == 'CO'
        assert line.nr is None
        assert line.lam is None
        assert line.freq is None
        assert line.a_stein is None
        assert line.e_up is None
        assert line.e_low is None
        assert line.g_up is None
        assert line.g_low is None

    def test_init_partial_fields(self):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('OH', {'lam': 10.5, 'freq': 2.85e13})
        assert line.lam == 10.5
        assert line.freq == 2.85e13
        assert line.nr is None
        assert line.a_stein is None

    def test_get_ndarray(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('H2O', sample_line_data)
        arr = line.get_ndarray()
        assert isinstance(arr, np.ndarray)
        assert len(arr) == 10
        # ndarray is object/string dtype because it mixes strings and numbers
        assert float(arr[0]) == 1       # nr
        assert float(arr[3]) == 12.407  # lam
        assert float(arr[4]) == 2.416e13  # freq

    def test_get_dict(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('H2O', sample_line_data)
        d = line.get_dict()
        assert isinstance(d, dict)
        assert d['nr'] == 1
        assert d['lam'] == 12.407
        assert d['freq'] == 2.416e13
        assert d['a_stein'] == 1.05e-2
        assert d['e_up'] == 4586.4
        assert d['e_low'] == 3379.0
        assert d['g_up'] == 21
        assert d['g_low'] == 19

    def test_get_pandas_table(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        import pandas as pd
        line = MoleculeLine('H2O', sample_line_data)
        df = line.get_pandas_table()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert df['lam'].iloc[0] == 12.407
        assert df['freq'].iloc[0] == 2.416e13

    def test_str_repr(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('H2O', sample_line_data)
        s = str(line)
        assert 'MoleculeLine' in s
        assert 'H2O' in s
        assert '12.407' in s
        assert repr(line) == str(line)

    def test_line_data_view_properties(self, sample_line_data):
        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
        line = MoleculeLine('H2O', sample_line_data)
        view = line.line_data
        assert view.nr == 1
        assert view.lam == 12.407
        assert view.freq == 2.416e13
        assert view.a_stein == 1.05e-2
        assert view.e_up == 4586.4
        assert view.e_low == 3379.0
        assert view.g_up == 21
        assert view.g_low == 19
        assert view.lev_up == '0_0_0|10_2_9'
        assert view.lev_low == '0_0_0|9_1_8'
