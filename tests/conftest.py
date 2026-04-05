# -*- coding: utf-8 -*-
"""
Shared fixtures for iSLAT DataTypes unit tests.
"""

import pytest
import numpy as np

import iSLAT.Constants as c


@pytest.fixture
def sample_line_data():
    """Sample molecular line data dict."""
    return {
        'nr': 1,
        'lev_up': '0_0_0|10_2_9',
        'lev_low': '0_0_0|9_1_8',
        'lam': 12.407,
        'freq': 2.416e13,
        'a_stein': 1.05e-2,
        'e_up': 4586.4,
        'e_low': 3379.0,
        'g_up': 21,
        'g_low': 19,
    }


@pytest.fixture
def sample_lines_data():
    """Sample list of line dicts for MoleculeLineList."""
    return [
        {
            'nr': 1, 'lev_up': '0|10', 'lev_low': '0|9',
            'lam': 12.407, 'freq': 2.416e13, 'a_stein': 1.05e-2,
            'e_up': 4586.4, 'e_low': 3379.0, 'g_up': 21, 'g_low': 19,
        },
        {
            'nr': 2, 'lev_up': '0|8', 'lev_low': '0|7',
            'lam': 14.950, 'freq': 2.005e13, 'a_stein': 2.10e-2,
            'e_up': 3200.0, 'e_low': 2500.0, 'g_up': 17, 'g_low': 15,
        },
        {
            'nr': 3, 'lev_up': '0|6', 'lev_low': '0|5',
            'lam': 17.221, 'freq': 1.741e13, 'a_stein': 5.00e-3,
            'e_up': 2100.0, 'e_low': 1600.0, 'g_up': 13, 'g_low': 11,
        },
        {
            'nr': 4, 'lev_up': '0|20', 'lev_low': '0|19',
            'lam': 6.500, 'freq': 4.612e13, 'a_stein': 3.00e-2,
            'e_up': 8000.0, 'e_low': 6500.0, 'g_up': 41, 'g_low': 39,
        },
        {
            'nr': 5, 'lev_up': '0|15', 'lev_low': '0|14',
            'lam': 22.500, 'freq': 1.332e13, 'a_stein': 8.00e-3,
            'e_up': 1200.0, 'e_low': 800.0, 'g_up': 31, 'g_low': 29,
        },
    ]


@pytest.fixture
def molecule_line_list(sample_lines_data):
    """Create a MoleculeLineList from sample data."""
    from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
    return MoleculeLineList(molecule_id='H2O', lines_data=sample_lines_data)


@pytest.fixture
def molecule_line_list_with_range(sample_lines_data):
    """MoleculeLineList with a wavelength range filter."""
    from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
    return MoleculeLineList(
        molecule_id='H2O',
        lines_data=sample_lines_data,
        wavelength_range=(10.0, 20.0),
    )
