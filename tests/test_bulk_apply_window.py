# -*- coding: utf-8 -*-
"""Unit tests for the BulkApplyPropertiesWindow helper logic."""
from iSLAT.Modules.GUI.Widgets.BulkApplyWindow import (
    MODE_ADD,
    MODE_SCALE,
    MODE_SET,
    UNCHANGED_CHOICE,
    apply_bulk_updates,
    apply_mode,
    build_choice_dict,
    build_parameter_dict,
    resolve_molecule_parameters,
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

def test_build_choice_dict_skips_unchanged():
    assert build_choice_dict({"instrumental_profile_key": UNCHANGED_CHOICE}) == {}
    assert build_choice_dict({"instrumental_profile_key": ""}) == {}

def test_build_choice_dict_maps_display_label_to_key():
    from iSLAT.Modules.DataProcessing.InstrumentalProfiles import PROFILE_DISPLAY_NAMES

    label = PROFILE_DISPLAY_NAMES["miri_mrs"]
    assert build_choice_dict({"instrumental_profile_key": label}) == {
        "instrumental_profile_key": "miri_mrs"
    }

def test_build_choice_dict_ignores_unknown_label():
    assert build_choice_dict({"instrumental_profile_key": "Bogus"}) == {}

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

def test_apply_mode_set_add_scale():
    assert apply_mode(100.0, 5.0, MODE_SET) == 5.0
    assert apply_mode(100.0, 5.0, MODE_ADD) == 105.0
    assert apply_mode(100.0, 5.0, MODE_SCALE) == 500.0

def test_apply_mode_skips_relative_on_non_numeric():
    assert apply_mode(None, 5.0, MODE_ADD) is None
    assert apply_mode("abc", 5.0, MODE_SCALE) is None
    assert apply_mode(None, 5.0, MODE_SET) == 5.0

def test_resolve_molecule_parameters_mixed_modes():
    mol = _Mol(temp=500.0, radius=2.0, n_mol=None)
    resolved = resolve_molecule_parameters(
        mol,
        {"temp": 100.0, "radius": 3.0, "n_mol": 2.0},
        {"temp": MODE_ADD, "radius": MODE_SCALE, "n_mol": MODE_ADD},
    )
    assert resolved == {"temp": 600.0, "radius": 6.0}

class _RecordingDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calls = []

    def bulk_update_parameters(self, params, molecules=None):
        self.calls.append((dict(params), list(molecules) if molecules else None))
        for name in molecules or self.keys():
            for attr, value in params.items():
                setattr(self[name], attr, value)

def test_apply_bulk_updates_set_mode_single_call():
    md = _RecordingDict({"H2O": _Mol(temp=500.0), "CO": _Mol(temp=300.0)})
    apply_bulk_updates(md, ["H2O", "CO"], {"temp": 800.0}, {"temp": MODE_SET})
    assert md.calls == [({"temp": 800.0}, ["H2O", "CO"])]

def test_apply_bulk_updates_relative_mode_per_molecule():
    md = _RecordingDict({"H2O": _Mol(temp=500.0), "CO": _Mol(temp=300.0)})
    apply_bulk_updates(md, ["H2O", "CO"], {"temp": 2.0}, {"temp": MODE_SCALE})
    assert md["H2O"].temp == 1000.0
    assert md["CO"].temp == 600.0

def test_apply_bulk_updates_mixes_absolute_relative_and_choices():
    md = _RecordingDict({"H2O": _Mol(temp=500.0, radius=2.0)})
    apply_bulk_updates(
        md,
        ["H2O"],
        {"temp": 100.0, "radius": 4.0},
        {"temp": MODE_ADD, "radius": MODE_SET},
        {"instrumental_profile_key": "miri_mrs"},
    )
    assert md.calls[0] == (
        {"radius": 4.0, "instrumental_profile_key": "miri_mrs"},
        ["H2O"],
    )
    assert md["H2O"].temp == 600.0