# -*- coding: utf-8 -*-
"""Unit tests for the BulkApplyPropertiesWindow helper logic."""
from iSLAT.Modules.GUI.Widgets.BulkApplyWindow import (
    build_parameter_dict,
    resolve_target_names,
)

class _Mol:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

def test_build_parameter_dict_skips_blanks():
    raw = {"temp": "700", "radius": "", "n_mol": "   "}
    params, invalid = build_parameter_dict(raw)
    assert params == {"temp": 700.0}
    assert invalid == []

def test_build_parameter_dict_parses_floats():
    raw = {"temp": "500.5", "n_mol": "1e18"}
    params, invalid = build_parameter_dict(raw)
    assert params == {"temp": 500.5, "n_mol": 1e18}
    assert invalid == []

def test_build_parameter_dict_reports_invalid():
    raw = {"temp": "abc", "radius": "2.0"}
    params, invalid = build_parameter_dict(raw)
    assert params == {"radius": 2.0}
    assert invalid == ["temp"]

def test_build_parameter_dict_empty():
    params, invalid = build_parameter_dict({"temp": "", "radius": ""})
    assert params == {}
    assert invalid == []

def _make_dict():
    return {
        "H2O": _Mol(name="H2O", is_visible=True),
        "CO": _Mol(name="CO", is_visible=False),
        "OH": _Mol(name="OH", is_visible=True),
    }

class _MolDict(dict):
    def get_visible_molecules(self):
        return {name for name, mol in self.items() if mol.is_visible}

def test_resolve_target_names_all():
    md = _make_dict()
    names = resolve_target_names(md, visible_only=False)
    assert set(names) == {"H2O", "CO", "OH"}

def test_resolve_target_names_visible_only_fallback():
    md = _make_dict()  # plain dict, no get_visible_molecules
    names = resolve_target_names(md, visible_only=True)
    assert set(names) == {"H2O", "OH"}

def test_resolve_target_names_visible_only_uses_method():
    md = _MolDict(_make_dict())
    names = resolve_target_names(md, visible_only=True)
    assert set(names) == {"H2O", "OH"}