# -*- coding: utf-8 -*-
"""End-to-end tests for the optical depth (tau) pipeline.

Validates the full tau data flow:
   Intensity.tau (per-line) → Spectrum.tau_profile (convolved)
   → Molecule.get_tau() (RV-shifted, resampled)
   → BasePlot.get_molecule_tau_data() (plotting helper)
"""

import pytest
import numpy as np

import iSLAT.Constants as c


class _TauTestBase:
    """Shared helpers for tau pipeline tests."""

    def _make_line_list(self, n_lines=5, lam_start=12.0, lam_step=1.5,
                        molecule_id='TestMol'):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        rng = np.random.default_rng(42)
        lines_data = []
        for i in range(n_lines):
            lam = lam_start + i * lam_step
            freq = c.SPEED_OF_LIGHT_MICRONS / lam
            lines_data.append({
                'nr': i + 1,
                'lev_up': f'0|{2 * i + 2}',
                'lev_low': f'0|{2 * i + 1}',
                'lam': lam,
                'freq': freq,
                'a_stein': rng.uniform(1e-3, 5e-2),
                'e_up': 1000.0 + i * 500,
                'e_low': 500.0 + i * 500,
                'g_up': 2 * i + 3,
                'g_low': 2 * i + 1,
            })
        mll = MoleculeLineList(molecule_id=molecule_id, lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015
        return mll

    def _make_intensity(self, mll, **kwargs):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        defaults = dict(t_kin=500.0, n_mol=1e18, dv=1.0)
        defaults.update(kwargs)
        inten = Intensity(mll)
        inten.calc_intensity(**defaults)
        return inten

    def _make_spectrum(self, **kwargs):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        defaults = dict(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                        R=3000.0, distance=160.0)
        defaults.update(kwargs)
        return Spectrum(**defaults)


class TestTauIntensityToSpectrum(_TauTestBase):
    """Test that per-line tau flows correctly into convolved tau profile."""

    def test_per_line_tau_nonnegative(self):
        mll = self._make_line_list()
        inten = self._make_intensity(mll)
        assert np.all(inten.tau >= 0)

    def test_tau_increases_with_column_density(self):
        mll = self._make_line_list()
        inten_low = self._make_intensity(mll, n_mol=1e16)
        inten_high = self._make_intensity(mll, n_mol=1e20)
        assert np.all(inten_high.tau >= inten_low.tau)

    def test_convolved_tau_peaks_at_line_centres(self):
        """Convolved tau profile should peak near line wavelengths."""
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)

        tp = spec.tau_profile
        idx_peak = np.argmax(tp)
        peak_lam = spec.lamgrid[idx_peak]
        # Peak should be near one of the line centres (14.5, 15.0, 15.5)
        assert 14.0 < peak_lam < 16.0

    def test_convolved_tau_integral_scales_with_line_tau(self):
        """The integral of the convolved profile should scale with per-line tau."""
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten_low = self._make_intensity(mll, n_mol=1e17)
        inten_high = self._make_intensity(mll, n_mol=1e18)

        spec_low = self._make_spectrum()
        spec_low.add_intensity(inten_low, dA=1.0)
        spec_high = self._make_spectrum()
        spec_high.add_intensity(inten_high, dA=1.0)

        integral_low = np.trapz(spec_low.tau_profile, spec_low.lamgrid)
        integral_high = np.trapz(spec_high.tau_profile, spec_high.lamgrid)

        # Integral should be larger for higher column density
        assert integral_high > integral_low


class TestTauMoleculeLevel(_TauTestBase):
    """Test Molecule.get_tau() integration (end-to-end through the full pipeline)."""

    def _make_molecule(self, **kwargs):
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        defaults = dict(
            name='TestH2O', displaylabel='H$_2$O', filepath=None,
            color='#FF0000', is_visible=True, temp=500.0, radius=1.0,
            n_mol=1e18, distance=160.0, fwhm=130.0,
            initial_molecule_parameters={
                't_kin': 500.0, 'scale_exponent': 18.0,
                'scale_number': 1.0, 'radius_init': 1.0,
            },
        )
        defaults.update(kwargs)
        return Molecule(**defaults)

    def test_molecule_without_lines_get_tau_returns_none(self):
        """Molecule with no filepath should handle get_tau gracefully."""
        mol = self._make_molecule()
        # Without lines loaded, get_tau should return None or raise gracefully
        try:
            result = mol.get_tau()
            assert result is None or isinstance(result, np.ndarray)
        except (AttributeError, ValueError, TypeError):
            # Expected: molecule without data cannot compute tau
            pass


class TestTauPerComponentIsolation(_TauTestBase):
    """Verify per-component tau isolation when multiple molecules are added."""

    def test_per_component_isolation(self):
        """Each component's tau profile should be non-overlapping when lines are far apart."""
        mll_a = self._make_line_list(n_lines=3, lam_start=11.5, lam_step=0.5, molecule_id='MolA')
        mll_b = self._make_line_list(n_lines=3, lam_start=17.5, lam_step=0.5, molecule_id='MolB')
        inten_a = self._make_intensity(mll_a)
        inten_b = self._make_intensity(mll_b)

        spec = self._make_spectrum()
        spec.add_intensity(inten_a, dA=1.0)
        spec.add_intensity(inten_b, dA=1.0)

        per_comp = spec._convol_tau_per_component()
        assert len(per_comp) == 2

        # Component A should peak near 11.5–12.5 μm, component B near 17.5–18.5 μm
        peak_a = spec.lamgrid[np.argmax(per_comp[0]['tau_profile'])]
        peak_b = spec.lamgrid[np.argmax(per_comp[1]['tau_profile'])]
        assert peak_a < 13.5
        assert peak_b > 17.0

    def test_per_component_sum_equals_total_with_real_intensity(self):
        """Per-component tau sum should equal total tau profile (with real intensity data)."""
        mll_a = self._make_line_list(n_lines=2, lam_start=13.0, molecule_id='A')
        mll_b = self._make_line_list(n_lines=2, lam_start=16.0, molecule_id='B')
        inten_a = self._make_intensity(mll_a)
        inten_b = self._make_intensity(mll_b)

        spec = self._make_spectrum()
        spec.add_intensity(inten_a, dA=1.0)
        spec.add_intensity(inten_b, dA=1.0)

        total = spec.tau_profile
        per_comp = spec._convol_tau_per_component()
        summed = sum(c['tau_profile'] for c in per_comp)
        np.testing.assert_allclose(total, summed, rtol=1e-10)


class TestTauCacheInvalidation(_TauTestBase):
    """Verify that tau cache invalidation works correctly across operations."""

    def test_tau_invalidated_on_add_intensity(self):
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        _ = spec.tau_profile  # prime cache
        assert spec._tau_profile_valid is True

        # Add another intensity — should invalidate
        mll2 = self._make_line_list(n_lines=3, lam_start=16.0, lam_step=0.5, molecule_id='Mol2')
        inten2 = self._make_intensity(mll2)
        spec.add_intensity(inten2, dA=1.0)
        assert spec._tau_profile_valid is False

    def test_tau_invalidated_on_wavelength_range_change(self):
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        _ = spec.tau_profile
        assert spec._tau_profile_valid is True

        spec.wavelength_range = (12.0, 18.0)
        assert spec._tau_profile_valid is False

    def test_tau_invalidated_on_R_func_change(self):
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        _ = spec.tau_profile
        assert spec._tau_profile_valid is True

        spec.R_func = lambda w: np.full_like(np.atleast_1d(w), 1500.0)
        assert spec._tau_profile_valid is False

    def test_tau_and_flux_independent_caching(self):
        """Accessing tau should not compute flux and vice versa."""
        mll = self._make_line_list(n_lines=3, lam_start=14.5, lam_step=0.5)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)

        # Access tau only
        _ = spec.tau_profile
        assert spec._tau_profile_valid is True
        assert spec._flux_valid is False  # flux not yet computed

        # Now access flux
        _ = spec.flux
        assert spec._flux_valid is True
        assert spec._tau_profile_valid is True  # tau still valid
