# -*- coding: utf-8 -*-
"""Unit tests for MoleculeDict and its helper functions."""

import pytest
import numpy as np

import iSLAT.Constants as c


class TestMoleculeDict:
    """Tests for MoleculeDict."""

    def _make_molecule(self, name='TestH2O', **kwargs):
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        defaults = {
            'name': name,
            'displaylabel': name,
            'filepath': None,
            'color': '#FF0000',
            'is_visible': True,
            'temp': 500.0,
            'radius': 1.0,
            'n_mol': 1e18,
            'distance': 160.0,
            'fwhm': 130.0,
            'initial_molecule_parameters': {
                't_kin': 500.0,
                'scale_exponent': 18.0,
                'scale_number': 1.0,
                'radius_init': 1.0,
            },
        }
        defaults.update(kwargs)
        return Molecule(**defaults)

    def _make_dict(self, **kwargs):
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        return MoleculeDict(**kwargs)

    def test_init_default(self):
        md = self._make_dict()
        assert isinstance(md, dict)
        assert len(md) == 0

    def test_init_global_parameters(self):
        md = self._make_dict(global_distance=200.0, global_stellar_rv=5.0)
        assert md._global_dist == 200.0
        assert md._global_stellar_rv == 5.0

    def test_add_molecule(self):
        md = self._make_dict()
        mol = self._make_molecule('H2O')
        md['H2O'] = mol
        assert 'H2O' in md
        assert len(md) == 1
        assert md['H2O'].name == 'H2O'

    def test_remove_molecule(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O')
        del md['H2O']
        assert 'H2O' not in md

    def test_get_molecule(self):
        md = self._make_dict()
        mol = self._make_molecule('OH')
        md['OH'] = mol
        result = md.get('OH')
        assert result is mol
        assert md.get('nonexistent') is None

    def test_get_visible_molecules(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O', is_visible=True)
        md['CO'] = self._make_molecule('CO', is_visible=False)
        md['OH'] = self._make_molecule('OH', is_visible=True)

        visible = md.get_visible_molecules()
        assert 'H2O' in visible
        assert 'OH' in visible
        assert 'CO' not in visible

    def test_get_visible_molecules_as_objects(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O', is_visible=True)
        md['CO'] = self._make_molecule('CO', is_visible=False)

        visible = md.get_visible_molecules(return_objects=True)
        assert len(visible) == 1
        assert visible[0].name == 'H2O'

    def test_bulk_set_visibility(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O', is_visible=True)
        md['CO'] = self._make_molecule('CO', is_visible=True)
        md['OH'] = self._make_molecule('OH', is_visible=True)

        md.bulk_set_visibility(False)
        for mol in md.values():
            assert mol.is_visible is False

        md.bulk_set_visibility(True, molecule_names=['H2O', 'OH'])
        assert md['H2O'].is_visible is True
        assert md['OH'].is_visible is True
        assert md['CO'].is_visible is False

    def test_apply_stellar_rv_zero(self):
        md = self._make_dict(global_stellar_rv=0.0)
        wave = np.linspace(10, 20, 100)
        result = md.apply_stellar_rv(wave)
        np.testing.assert_array_equal(result, wave)

    def test_apply_stellar_rv_nonzero(self):
        md = self._make_dict(global_stellar_rv=10.0)
        wave = np.array([10.0, 15.0, 20.0])
        result = md.apply_stellar_rv(wave)
        # With positive RV, wavelengths should shift
        assert not np.array_equal(result, wave)
        # The formula: wave - (wave / c_kms * rv)
        expected = wave - (wave / c.SPEED_OF_LIGHT_KMS * 10.0)
        np.testing.assert_allclose(result, expected)

    def test_clear(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O')
        md['CO'] = self._make_molecule('CO')
        md.clear()
        assert len(md) == 0

    def test_iter_and_contains(self):
        md = self._make_dict()
        md['H2O'] = self._make_molecule('H2O')
        md['CO'] = self._make_molecule('CO')
        keys = list(md)
        assert 'H2O' in keys
        assert 'CO' in keys
        assert 'H2O' in md
        assert 'OH' not in md

    def test_values_and_items(self):
        md = self._make_dict()
        mol = self._make_molecule('H2O')
        md['H2O'] = mol
        vals = list(md.values())
        assert len(vals) == 1
        assert vals[0] is mol
        items = list(md.items())
        assert items[0] == ('H2O', mol)

    def test_global_parameter_defaults(self):
        md = self._make_dict()
        assert md._global_parms['dist'] == c.DEFAULT_DISTANCE
        assert md._global_parms['stellar_rv'] == c.DEFAULT_STELLAR_RV
        assert md._global_parms['wavelength_range'] == c.WAVELENGTH_RANGE

    def test_matched_spectral_sampling_off_by_default(self):
        md = self._make_dict()
        assert md._match_spectral_sampling is False

    def test_get_matched_sampling_wavelengths_off(self):
        md = self._make_dict()
        wave = np.linspace(10, 20, 50)
        use_interp, target = md.get_matched_sampling_wavelengths(wave)
        assert use_interp is False
        np.testing.assert_array_equal(target, wave)

    def test_get_matched_sampling_wavelengths_on(self):
        md = self._make_dict(global_stellar_rv=10.0)
        md._match_spectral_sampling = True
        wave = np.linspace(10, 20, 50)
        use_interp, target = md.get_matched_sampling_wavelengths(wave)
        assert use_interp is True
        # target should be RV-corrected wavelengths
        assert not np.array_equal(target, wave)

    # ── global_wavelength_range setter ─────────────────────────────

    def test_global_wavelength_range_no_intermediate_class_callbacks(self):
        """Setting global_wavelength_range must NOT fire per-molecule
        class callbacks (which would trigger intermediate plot refreshes
        before all molecules are consistent)."""
        from iSLAT.Modules.DataTypes.Molecule import Molecule

        md = self._make_dict()
        mol_a = self._make_molecule('MolA')
        mol_b = self._make_molecule('MolB')
        md['MolA'] = mol_a
        md['MolB'] = mol_b

        original_range = (4.5, 28.0)
        md._global_wavelength_range = original_range
        mol_a._wavelength_range = original_range
        mol_b._wavelength_range = original_range

        # Track per-molecule class callbacks
        class_cb_calls = []
        def class_cb(mol_name, param_name, old_val, new_val):
            class_cb_calls.append((mol_name, param_name))

        Molecule.add_molecule_parameter_change_callback(class_cb)
        try:
            md.global_wavelength_range = (5.0, 20.0)

            # No per-molecule class callbacks should have fired
            wl_calls = [c for c in class_cb_calls if c[1] == 'wavelength_range']
            assert len(wl_calls) == 0, (
                f"Expected 0 per-molecule wavelength_range callbacks, got {len(wl_calls)}"
            )

            # But every molecule should have the new range and dirty flags
            for mol in md.values():
                assert mol._wavelength_range == (5.0, 20.0)
                assert mol._dirty_flags['intensity'] is True
                assert mol._dirty_flags['spectrum'] is True
                assert mol._dirty_flags['flux'] is True
        finally:
            Molecule.remove_molecule_parameter_change_callback(class_cb)

    def test_global_wavelength_range_fires_global_callback(self):
        """The global callback should fire exactly once after all molecules
        are updated."""
        md = self._make_dict()
        mol_a = self._make_molecule('MolA')
        md['MolA'] = mol_a
        md._global_wavelength_range = (4.5, 28.0)
        mol_a._wavelength_range = (4.5, 28.0)

        global_cb_calls = []
        md.add_global_parameter_change_callback(
            lambda p, o, n: global_cb_calls.append((p, o, n))
        )
        md.global_wavelength_range = (5.0, 20.0)

        assert len(global_cb_calls) == 1
        assert global_cb_calls[0] == ('wavelength_range', (4.5, 28.0), (5.0, 20.0))


class TestMoleculeDictHelpers:
    """Tests for MoleculeDict helper functions (_ci_get, _safe_float)."""

    def test_ci_get_exact(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _ci_get
        d = {'Temp': 500, 'Name': 'H2O'}
        assert _ci_get(d, 'Temp') == 500
        assert _ci_get(d, 'Name') == 'H2O'

    def test_ci_get_case_insensitive(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _ci_get
        d = {'Temp': 500, 'Name': 'H2O'}
        assert _ci_get(d, 'temp') == 500
        assert _ci_get(d, 'name') == 'H2O'
        assert _ci_get(d, 'TEMP') == 500

    def test_ci_get_missing(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _ci_get
        d = {'Temp': 500}
        assert _ci_get(d, 'missing') is None

    def test_safe_float_basic(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'Temp': '500.5', 'N_mol': 1e18}
        assert _safe_float(d, 'Temp') == 500.5
        assert _safe_float(d, 'N_mol') == 1e18

    def test_safe_float_default(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'Temp': '500'}
        assert _safe_float(d, 'missing', default=42.0) == 42.0

    def test_safe_float_list_keys(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'temperature': 300.0}
        # First key misses, second hits
        assert _safe_float(d, ['Temp', 'temperature']) == 300.0

    def test_safe_float_invalid_value(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'Temp': 'not_a_number'}
        assert _safe_float(d, 'Temp', default=-1.0) == -1.0

    def test_safe_float_none_value(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'Temp': None}
        assert _safe_float(d, 'Temp', default=0.0) == 0.0

    def test_safe_float_case_insensitive(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import _safe_float
        d = {'Temperature': 400.0}
        assert _safe_float(d, 'temperature', case_insensitive=True) == 400.0
