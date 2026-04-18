# -*- coding: utf-8 -*-
"""Unit tests for Molecule."""

import pytest
import numpy as np

import iSLAT.Constants as c


class TestMolecule:
    """Tests for the Molecule class."""

    def _make_molecule(self, **kwargs):
        """Helper to create a Molecule with defaults."""
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        defaults = {
            'name': 'TestH2O',
            'displaylabel': 'H$_2$O',
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

    def test_init_basic(self):
        mol = self._make_molecule()
        assert mol.name == 'TestH2O'
        assert mol.displaylabel == 'H$_2$O'
        assert mol.color == '#FF0000'

    def test_temp_property(self):
        mol = self._make_molecule(temp=500.0)
        assert mol.temp == 500.0
        mol.temp = 600.0
        assert mol.temp == 600.0

    def test_radius_property(self):
        mol = self._make_molecule(radius=1.0)
        assert mol.radius == 1.0
        mol.radius = 2.5
        assert mol.radius == 2.5

    def test_n_mol_property(self):
        mol = self._make_molecule(n_mol=1e18)
        assert mol.n_mol == 1e18
        mol.n_mol = 5e19
        assert mol.n_mol == 5e19

    def test_distance_property(self):
        mol = self._make_molecule(distance=160.0)
        assert mol.distance == 160.0
        mol.distance = 200.0
        assert mol.distance == 200.0

    def test_fwhm_property(self):
        mol = self._make_molecule(fwhm=130.0)
        assert mol.fwhm == 130.0
        mol.fwhm = 150.0
        assert mol.fwhm == 150.0

    def test_broad_property(self):
        mol = self._make_molecule()
        mol.broad = 2.0
        assert mol.broad == 2.0

    def test_rv_shift_property(self):
        mol = self._make_molecule()
        assert mol.rv_shift == 0.0
        mol.rv_shift = 5.0
        assert mol.rv_shift == 5.0

    def test_is_visible_bool(self):
        mol = self._make_molecule(is_visible=True)
        assert mol.is_visible is True
        mol.is_visible = False
        assert mol.is_visible is False

    def test_is_visible_string_conversion(self):
        mol = self._make_molecule(is_visible='true')
        assert mol.is_visible is True
        mol.is_visible = 'false'
        assert mol.is_visible is False
        mol.is_visible = '1'
        assert mol.is_visible is True
        mol.is_visible = '0'
        assert mol.is_visible is False

    def test_wavelength_range_property(self):
        mol = self._make_molecule()
        mol.wavelength_range = (5.0, 25.0)
        assert mol.wavelength_range == (5.0, 25.0)

    def test_intensity_calculation_method_property(self):
        mol = self._make_molecule()
        assert mol.intensity_calculation_method == 'curve_growth'
        mol.intensity_calculation_method = 'radex'
        assert mol.intensity_calculation_method == 'radex'

    def test_str_repr(self):
        mol = self._make_molecule(name='TestMol', temp=300, radius=0.5, n_mol=1e16)
        s = str(mol)
        assert 'Molecule' in s
        assert 'TestMol' in s

    def test_caching_system_initialized(self):
        mol = self._make_molecule()
        assert isinstance(mol._dirty_flags, dict)
        assert 'intensity' in mol._dirty_flags
        assert 'spectrum' in mol._dirty_flags
        assert 'flux' in mol._dirty_flags

    def test_clear_all_caches(self):
        mol = self._make_molecule()
        mol._dirty_flags['intensity'] = False
        mol.clear_all_caches()
        assert mol._dirty_flags['intensity'] is True
        assert mol._dirty_flags['spectrum'] is True
        assert mol._dirty_flags['flux'] is True

    def test_is_cache_valid(self):
        mol = self._make_molecule()
        # All dirty initially
        assert mol.is_cache_valid('intensity') is False
        assert mol.is_cache_valid('spectrum') is False
        assert mol.is_cache_valid('flux') is False
        assert mol.is_cache_valid('full') is False

    def test_parameter_hash_changes_on_update(self):
        mol = self._make_molecule()
        h1 = mol._compute_intensity_hash()
        mol.temp = 700.0
        h2 = mol._compute_intensity_hash()
        assert h1 != h2

    def test_bulk_update_parameters(self):
        mol = self._make_molecule()
        mol.bulk_update_parameters({'temp': 800.0, 'radius': 3.0, 'n_mol': 2e19})
        assert mol.temp == 800.0
        assert mol.radius == 3.0
        assert mol.n_mol == 2e19

    def test_bulk_update_parameters_empty(self):
        mol = self._make_molecule(temp=500.0)
        mol.bulk_update_parameters({})
        assert mol.temp == 500.0  # unchanged

    def test_bulk_update_parameters_skip_notification(self):
        mol = self._make_molecule()
        # Should not raise even with skip_notification
        mol.bulk_update_parameters({'temp': 999.0}, skip_notification=True)
        assert mol.temp == 999.0

    def test_dirty_flags_on_temp_change(self):
        mol = self._make_molecule()
        mol._dirty_flags['intensity'] = False
        mol._dirty_flags['spectrum'] = False
        mol._dirty_flags['flux'] = False
        mol.temp = 700.0
        # temp is intensity-affecting
        assert mol._dirty_flags['intensity'] is True
        assert mol._dirty_flags['spectrum'] is True
        assert mol._dirty_flags['flux'] is True

    def test_dirty_flags_on_radius_change(self):
        mol = self._make_molecule()
        mol._dirty_flags['intensity'] = False
        mol._dirty_flags['spectrum'] = False
        mol._dirty_flags['flux'] = False
        mol.radius = 3.0
        # radius is spectrum-affecting only
        assert mol._dirty_flags['intensity'] is False
        assert mol._dirty_flags['spectrum'] is True
        assert mol._dirty_flags['flux'] is True

    def test_init_from_user_save_data(self):
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        usd = {
            'Molecule Name': 'SavedMol',
            'File Path': None,
            'Molecule Label': 'Saved',
            'Temp': 400.0,
            'Rad': 2.0,
            'N_Mol': 5e17,
            'Color': '#00FF00',
            'Vis': True,
            'Dist': 140.0,
            'FWHM': 120.0,
            'Broad': 1.5,
            'RV Shift': 3.0,
        }
        mol = Molecule(
            user_save_data=usd,
            initial_molecule_parameters={
                't_kin': 400.0,
                'scale_exponent': 17.0,
                'scale_number': 5.0,
                'radius_init': 2.0,
            },
        )
        assert mol.name == 'SavedMol'
        assert mol.temp == 400.0
        assert mol.radius == 2.0
        assert mol.color == '#00FF00'

    def test_molecule_parameter_change_callbacks(self):
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        callback_log = []

        def my_callback(mol_name, param, old, new):
            callback_log.append((mol_name, param, old, new))

        Molecule.add_molecule_parameter_change_callback(my_callback)
        try:
            mol = self._make_molecule(name='CallbackTest')
            mol.temp = 600.0
            # Should have recorded temp change
            assert any(entry[1] == 'temp' for entry in callback_log)
        finally:
            Molecule.remove_molecule_parameter_change_callback(my_callback)

    def test_get_cache_stats(self):
        mol = self._make_molecule()
        stats = mol.get_cache_stats()
        assert 'hits' in stats
        assert 'misses' in stats
        assert 'invalidations' in stats

    def test_no_filepath(self):
        """Molecule with no filepath should not crash."""
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        mol = Molecule(
            name='NoFile', filepath=None,
            initial_molecule_parameters={
                't_kin': 300.0,
                'scale_exponent': 1.0,
                'scale_number': 1.0,
                'radius_init': 1.0,
            },
        )
        assert mol.name == 'NoFile'
        assert mol.filepath is None
