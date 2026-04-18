# -*- coding: utf-8 -*-
"""Integration tests that combine multiple data types, plus helper and edge-case tests."""

import pytest
import numpy as np

import iSLAT.Constants as c


# ============================================================================
# Integration Tests (cross-type interactions)
# ============================================================================

class TestIntegration:
    """Integration tests that combine multiple data types."""

    def test_molecule_line_list_to_intensity(self):
        """MoleculeLineList should feed into Intensity correctly."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        lines_data = [
            {'nr': 1, 'lev_up': '0|4', 'lev_low': '0|3',
             'lam': 14.0, 'freq': c.SPEED_OF_LIGHT_MICRONS / 14.0,
             'a_stein': 0.01, 'e_up': 3000, 'e_low': 2000, 'g_up': 9, 'g_low': 7},
            {'nr': 2, 'lev_up': '0|6', 'lev_low': '0|5',
             'lam': 16.0, 'freq': c.SPEED_OF_LIGHT_MICRONS / 16.0,
             'a_stein': 0.02, 'e_up': 4000, 'e_low': 3000, 'g_up': 13, 'g_low': 11},
        ]
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015

        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        assert inten.molecule is mll
        assert len(inten.intensity) == 2
        assert np.all(inten.intensity >= 0)

    def test_intensity_to_spectrum(self):
        """Intensity should feed into Spectrum via add_intensity."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum

        lines_data = [
            {'nr': 1, 'lev_up': '0|4', 'lev_low': '0|3',
             'lam': 14.0, 'freq': c.SPEED_OF_LIGHT_MICRONS / 14.0,
             'a_stein': 0.01, 'e_up': 3000, 'e_low': 2000, 'g_up': 9, 'g_low': 7},
            {'nr': 2, 'lev_up': '0|6', 'lev_low': '0|5',
             'lam': 16.0, 'freq': c.SPEED_OF_LIGHT_MICRONS / 16.0,
             'a_stein': 0.02, 'e_up': 4000, 'e_low': 3000, 'g_up': 13, 'g_low': 11},
        ]
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015

        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        spec = Spectrum(
            lam_min=10.0, lam_max=20.0, dlambda=0.01,
            R=3000.0, distance=160.0,
        )
        area = np.pi * 1.0 ** 2  # 1 AU radius
        spec.add_intensity(inten, area)

        flux = spec.flux
        assert len(flux) == len(spec.lamgrid)
        # Should have non-zero flux around the line wavelengths
        assert np.max(flux) > 0

    def test_molecule_dict_with_multiple_molecules(self):
        """MoleculeDict should manage multiple Molecule objects."""
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

        md = MoleculeDict()

        for name, temp in [('H2O', 500), ('CO', 300), ('OH', 700)]:
            mol = Molecule(
                name=name, displaylabel=name, filepath=None,
                color='#FF0000', is_visible=True,
                temp=float(temp), radius=1.0, n_mol=1e18,
                distance=160.0, fwhm=130.0,
                initial_molecule_parameters={
                    't_kin': float(temp),
                    'scale_exponent': 18.0,
                    'scale_number': 1.0,
                    'radius_init': 1.0,
                },
            )
            md[name] = mol

        assert len(md) == 3
        visible = md.get_visible_molecules()
        assert len(visible) == 3

        md['CO'].is_visible = False
        visible = md.get_visible_molecules()
        assert len(visible) == 2
        assert 'CO' not in visible


# ============================================================================
# Helper Function Tests (from Molecule.py)
# ============================================================================

class TestSpectresResampling:
    """Tests for the _spectres flux-conserving resampling helper."""

    def test_identity_resampling(self):
        """Resampling onto the same grid should return ~same fluxes."""
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres as _spectres
        wavs = np.linspace(10, 20, 100)
        flux = np.sin(wavs)
        result = _spectres(wavs, wavs, flux)
        np.testing.assert_allclose(result, flux, atol=1e-10)

    def test_resampling_preserves_constant(self):
        """Resampling a constant spectrum should give the same constant."""
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres as _spectres
        old_wavs = np.linspace(10, 20, 200)
        new_wavs = np.linspace(11, 19, 50)
        constant_flux = np.ones_like(old_wavs) * 5.0
        result = _spectres(new_wavs, old_wavs, constant_flux)
        np.testing.assert_allclose(result, 5.0, atol=1e-10)

    def test_out_of_range_fill(self):
        """New wavelengths outside old range should be filled with fill value."""
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres as _spectres
        old_wavs = np.linspace(10, 20, 100)
        new_wavs = np.array([5.0, 15.0, 25.0])
        flux = np.ones_like(old_wavs)
        result = _spectres(new_wavs, old_wavs, flux, fill=-999.0)
        assert result[0] == -999.0  # Below range
        assert result[2] == -999.0  # Above range
        np.testing.assert_allclose(result[1], 1.0, atol=1e-10)


# ============================================================================
# Edge Cases
# ============================================================================

class TestEdgeCases:
    """Edge-case and boundary tests."""

    def test_molecule_line_list_single_line(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        mll = MoleculeLineList(molecule_id='single', lines_data=[{
            'nr': 1, 'lev_up': '0|2', 'lev_low': '0|1',
            'lam': 15.0, 'freq': 2e13, 'a_stein': 0.01,
            'e_up': 2000, 'e_low': 1000, 'g_up': 5, 'g_low': 3,
        }])
        assert mll.num_lines == 1
        assert mll.get_wavelengths()[0] == 15.0

    def test_molecule_line_list_duplicate_wavelengths(self):
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

    def test_molecule_no_filepath(self):
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
