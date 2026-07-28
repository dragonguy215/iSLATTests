# -*- coding: utf-8 -*-
"""Unit tests for the SortMoleculesWindow sort logic."""

from iSLAT.Modules.GUI.Widgets.SortMoleculesWindow import sort_molecule_names


class _Mol:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


def _make_dict():
    return {
        "H2O": _Mol(name="H2O", temp=700.0, radius=2.0, n_mol=1e18),
        "co": _Mol(name="co", temp=300.0, radius=1.0, n_mol=1e20),
        "OH": _Mol(name="OH", temp=500.0, radius=3.0, n_mol=1e19),
    }


def test_sort_numeric_ascending():
    md = _make_dict()
    order = sort_molecule_names(md, "temp", is_numeric=True, descending=False)
    assert order == ["co", "OH", "H2O"]


def test_sort_numeric_descending():
    md = _make_dict()
    order = sort_molecule_names(md, "radius", is_numeric=True, descending=True)
    assert order == ["OH", "H2O", "co"]


def test_sort_name_case_insensitive():
    md = _make_dict()
    order = sort_molecule_names(md, "name", is_numeric=False, descending=False)
    assert order == ["co", "H2O", "OH"]


def test_sort_missing_value_treated_as_zero():
    md = _make_dict()
    md["NoTemp"] = _Mol(name="NoTemp")  # no temp attribute
    order = sort_molecule_names(md, "temp", is_numeric=True, descending=False)
    assert order[0] == "NoTemp"


def test_sort_preserves_all_names():
    md = _make_dict()
    order = sort_molecule_names(md, "n_mol", is_numeric=True, descending=True)
    assert set(order) == set(md.keys())
    assert len(order) == len(md)
