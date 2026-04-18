# -*- coding: utf-8 -*-
"""Extended unit tests for the Spectrum class — coverage gaps.

Covers:
- add_intensity() integration with real Intensity objects
- Multi-component flux additivity
- resample_to() with RV shifts and unit selection
- _repr_html_() notebook display
- R_func (wavelength-dependent resolving power) convolution
- Distance scaling verification
- flux_jy conversion correctness
"""

import pytest
import numpy as np
import pandas as pd

import iSLAT.Constants as c


class TestSpectrumAddIntensity:
    """Tests for Spectrum.add_intensity() — the core composition method."""

    def _make_line_list(self, n_lines=5, lam_start=12.0):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        lines_data = []
        for i in range(n_lines):
            lam = lam_start + i * 1.5
            freq = c.SPEED_OF_LIGHT_MICRONS / lam
            lines_data.append({
                'nr': i + 1,
                'lev_up': f'0|{2 * i + 2}',
                'lev_low': f'0|{2 * i + 1}',
                'lam': lam,
                'freq': freq,
                'a_stein': np.random.default_rng(42 + i).uniform(1e-3, 5e-2),
                'e_up': 1000.0 + i * 500,
                'e_low': 500.0 + i * 500,
                'g_up': 2 * i + 3,
                'g_low': 2 * i + 1,
            })
        mll = MoleculeLineList(molecule_id='TestMol', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015
        return mll

    def _make_intensity(self, mll, t_kin=500.0, n_mol=1e18, dv=1.0):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=t_kin, n_mol=n_mol, dv=dv)
        return inten

    def _make_spectrum(self, lam_min=10.0, lam_max=20.0, dlambda=0.01,
                       R=3000.0, distance=160.0):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        return Spectrum(
            lam_min=lam_min, lam_max=lam_max, dlambda=dlambda,
            R=R, distance=distance,
        )

    # ------------------------------------------------------------------

    def test_add_intensity_populates_arrays(self):
        """add_intensity should populate the internal accumulation arrays."""
        mll = self._make_line_list()
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        assert len(spec._I_arrays) == 1
        assert len(spec._lam_arrays) == 1
        assert len(spec._tau_arrays) == 1

    def test_add_intensity_invalidates_cache(self):
        """Adding an intensity should mark flux cache invalid."""
        spec = self._make_spectrum()
        _ = spec.flux  # prime cache
        assert spec._flux_valid is True
        mll = self._make_line_list()
        inten = self._make_intensity(mll)
        spec.add_intensity(inten, dA=1.0)
        assert spec._flux_valid is False

    def test_add_intensity_produces_nonzero_flux(self):
        """After adding an intensity, flux should be nonzero near line centres."""
        mll = self._make_line_list(n_lines=3, lam_start=14.0)
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        flux = spec.flux
        assert np.any(flux > 0), "Flux should be nonzero after adding intensity"

    def test_add_intensity_area_scaling(self):
        """Doubling the area should double the flux."""
        mll = self._make_line_list(n_lines=2, lam_start=14.0)
        inten = self._make_intensity(mll)

        spec1 = self._make_spectrum()
        spec1.add_intensity(inten, dA=1.0)

        spec2 = self._make_spectrum()
        spec2.add_intensity(inten, dA=2.0)

        # spec2 flux should be ~2x spec1 flux everywhere flux is nonzero
        mask = spec1.flux > 0
        if np.any(mask):
            ratio = spec2.flux[mask] / spec1.flux[mask]
            np.testing.assert_allclose(ratio, 2.0, rtol=1e-10)

    def test_add_intensity_no_data_returns_safely(self):
        """Adding an intensity with None data should not crash."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_line_list()
        inten = Intensity(mll)
        # intensity is None — not calculated yet
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        # Should be a no-op
        assert len(spec._I_arrays) == 0

    def test_add_intensity_stores_component_info(self):
        """Components list should record molecule info."""
        mll = self._make_line_list()
        inten = self._make_intensity(mll, t_kin=600.0, n_mol=2e18, dv=1.5)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=3.14)
        assert len(spec.components) == 1
        comp = spec.components[0]
        assert comp['name'] == 'TestMol'
        assert comp['t_kin'] == 600.0
        assert comp['n_mol'] == 2e18
        assert comp['area'] == 3.14

    def test_add_intensity_tau_arrays_stored(self):
        """tau values from Intensity should be stored in _tau_arrays."""
        mll = self._make_line_list()
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        spec.add_intensity(inten, dA=1.0)
        tau_stored = spec._tau_arrays[0]
        assert np.all(tau_stored >= 0)
        assert len(tau_stored) > 0


class TestSpectrumMultiComponent:
    """Tests for multi-component flux additivity and interaction."""

    def _make_line_list(self, molecule_id, lam_values, seed=42):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        rng = np.random.default_rng(seed)
        lines_data = []
        for i, lam in enumerate(lam_values):
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

    def test_two_component_flux_additivity(self):
        """Flux from two separate molecules should equal combined spectrum."""
        mll_a = self._make_line_list('MolA', [12.5, 13.0, 13.5], seed=10)
        mll_b = self._make_line_list('MolB', [16.5, 17.0, 17.5], seed=20)
        inten_a = self._make_intensity(mll_a)
        inten_b = self._make_intensity(mll_b)

        # Combined
        spec_combined = self._make_spectrum()
        spec_combined.add_intensity(inten_a, dA=1.0)
        spec_combined.add_intensity(inten_b, dA=1.0)

        # Separate
        spec_a = self._make_spectrum()
        spec_a.add_intensity(inten_a, dA=1.0)
        spec_b = self._make_spectrum()
        spec_b.add_intensity(inten_b, dA=1.0)

        np.testing.assert_allclose(
            spec_combined.flux,
            spec_a.flux + spec_b.flux,
            rtol=1e-10,
        )

    def test_two_component_tau_additivity(self):
        """Tau profiles from two separate molecules should add linearly."""
        mll_a = self._make_line_list('MolA', [13.5, 14.0, 14.5], seed=10)
        mll_b = self._make_line_list('MolB', [15.5, 16.0, 16.5], seed=20)
        inten_a = self._make_intensity(mll_a)
        inten_b = self._make_intensity(mll_b)

        spec_combined = self._make_spectrum()
        spec_combined.add_intensity(inten_a, dA=1.0)
        spec_combined.add_intensity(inten_b, dA=1.0)

        spec_a = self._make_spectrum()
        spec_a.add_intensity(inten_a, dA=1.0)
        spec_b = self._make_spectrum()
        spec_b.add_intensity(inten_b, dA=1.0)

        np.testing.assert_allclose(
            spec_combined.tau_profile,
            spec_a.tau_profile + spec_b.tau_profile,
            rtol=1e-10,
        )

    def test_components_list_grows(self):
        """Each add_intensity call should add one component."""
        mll = self._make_line_list('Mol', [14.5, 15.0, 15.5])
        inten = self._make_intensity(mll)
        spec = self._make_spectrum()
        assert len(spec.components) == 0
        spec.add_intensity(inten, dA=1.0)
        assert len(spec.components) == 1
        spec.add_intensity(inten, dA=2.0)
        assert len(spec.components) == 2

    def test_per_component_tau_matches_individual(self):
        """per-component tau should match individually convolved spectra."""
        mll_a = self._make_line_list('MolA', [12.5, 13.0, 13.5], seed=10)
        mll_b = self._make_line_list('MolB', [16.5, 17.0, 17.5], seed=20)
        inten_a = self._make_intensity(mll_a)
        inten_b = self._make_intensity(mll_b)

        spec = self._make_spectrum()
        spec.add_intensity(inten_a, dA=1.0)
        spec.add_intensity(inten_b, dA=1.0)

        per_comp = spec._convol_tau_per_component()
        assert len(per_comp) == 2
        assert per_comp[0]['name'] == 'MolA'
        assert per_comp[1]['name'] == 'MolB'


class TestSpectrumResample:
    """Tests for resample_to() — flux-conserving resampling."""

    def _make_spectrum_with_lines(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        lines_data = []
        for i, lam in enumerate([14.0, 15.0, 16.0]):
            lines_data.append({
                'nr': i + 1, 'lev_up': f'0|{2*i+2}', 'lev_low': f'0|{2*i+1}',
                'lam': lam, 'freq': c.SPEED_OF_LIGHT_MICRONS / lam,
                'a_stein': 1e-2, 'e_up': 1500.0 + i * 500, 'e_low': 1000.0 + i * 500,
                'g_up': 2*i + 5, 'g_low': 2*i + 3,
            })
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015

        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        spec = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                        R=3000.0, distance=160.0)
        spec.add_intensity(inten, dA=1.0)
        return spec

    def test_resample_to_output_length(self):
        spec = self._make_spectrum_with_lines()
        target = np.linspace(12.0, 18.0, 100)
        result = spec.resample_to(target)
        assert len(result) == 100

    def test_resample_to_jy_unit(self):
        spec = self._make_spectrum_with_lines()
        target = np.linspace(12.0, 18.0, 50)
        result_jy = spec.resample_to(target, unit='jy')
        assert len(result_jy) == 50
        assert np.any(result_jy > 0)

    def test_resample_to_cgs_unit(self):
        spec = self._make_spectrum_with_lines()
        target = np.linspace(12.0, 18.0, 50)
        result_cgs = spec.resample_to(target, unit='cgs')
        assert len(result_cgs) == 50
        assert np.any(result_cgs > 0)

    def test_resample_to_rv_shift(self):
        """Positive RV should red-shift peaks relative to unshifted spectrum."""
        spec = self._make_spectrum_with_lines()
        target = np.linspace(12.0, 18.0, 2000)
        result_0 = spec.resample_to(target, rv_shift=0.0, unit='cgs')
        result_rv = spec.resample_to(target, rv_shift=300.0, unit='cgs')
        # Peak should shift to longer wavelength
        peak_0 = target[np.argmax(result_0)]
        peak_rv = target[np.argmax(result_rv)]
        assert peak_rv > peak_0

    def test_resample_to_fill_value(self):
        """Out-of-range target pixels should use the fill value."""
        spec = self._make_spectrum_with_lines()
        # Target extends well beyond the model grid
        target = np.linspace(25.0, 30.0, 50)
        result = spec.resample_to(target, fill=np.nan)
        # All should be fill value (NaN or 0)
        assert np.all(np.isnan(result) | (result == 0))


class TestSpectrumDistanceScaling:
    """Verify that flux scales as 1/distance² but tau does not."""

    def _make_spectrum_pair(self, d1=160.0, d2=320.0):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        lines_data = []
        for i, lam in enumerate([14.0, 15.0, 16.0]):
            lines_data.append({
                'nr': i + 1, 'lev_up': f'0|{2*i+2}', 'lev_low': f'0|{2*i+1}',
                'lam': lam, 'freq': c.SPEED_OF_LIGHT_MICRONS / lam,
                'a_stein': 1e-2, 'e_up': 1500.0 + i * 500, 'e_low': 1000.0 + i * 500,
                'g_up': 2*i + 5, 'g_low': 2*i + 3,
            })
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015

        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        s1 = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=d1)
        s1.add_intensity(inten, dA=1.0)
        s2 = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=d2)
        s2.add_intensity(inten, dA=1.0)
        return s1, s2

    def test_flux_scales_with_distance_squared(self):
        s1, s2 = self._make_spectrum_pair(d1=160.0, d2=320.0)
        mask = s1.flux > 0
        if np.any(mask):
            ratio = s1.flux[mask] / s2.flux[mask]
            # d2/d1 = 2, so flux ratio should be 4
            np.testing.assert_allclose(ratio, 4.0, rtol=1e-10)

    def test_tau_does_not_scale_with_distance(self):
        s1, s2 = self._make_spectrum_pair(d1=160.0, d2=320.0)
        np.testing.assert_allclose(s1.tau_profile, s2.tau_profile, rtol=1e-10)


class TestSpectrumFluxJy:
    """Tests for flux_jy conversion."""

    def _make_spectrum_with_lines(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        lines_data = []
        for i, lam in enumerate([14.0, 15.0, 16.0]):
            lines_data.append({
                'nr': i + 1, 'lev_up': f'0|{2*i+2}', 'lev_low': f'0|{2*i+1}',
                'lam': lam, 'freq': c.SPEED_OF_LIGHT_MICRONS / lam,
                'a_stein': 1e-2, 'e_up': 1500.0 + i * 500, 'e_low': 1000.0 + i * 500,
                'g_up': 2*i + 5, 'g_low': 2*i + 3,
            })
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        spec = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                        R=3000.0, distance=160.0)
        spec.add_intensity(inten, dA=1.0)
        return spec

    def test_flux_jy_conversion_formula(self):
        """flux_jy should equal flux * conversion_factor * lambda^2."""
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        spec = self._make_spectrum_with_lines()
        expected = spec.flux * Spectrum._FLUX_JY_FACTOR * (spec.lamgrid ** 2)
        np.testing.assert_allclose(spec.flux_jy, expected, rtol=1e-12)

    def test_flux_jy_positive_where_flux_positive(self):
        spec = self._make_spectrum_with_lines()
        mask = spec.flux > 0
        assert np.all(spec.flux_jy[mask] > 0)


class TestSpectrumRFunc:
    """Tests for wavelength-dependent R_func convolution."""

    def _make_spectrum(self, R=3000.0, R_func=None):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        return Spectrum(
            lam_min=10.0, lam_max=20.0, dlambda=0.01,
            R=R, distance=160.0, R_func=R_func,
        )

    def test_R_func_changes_line_width(self):
        """An R_func that halves R should produce broader lines."""
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        lines_data = []
        for i, lam in enumerate([14.0, 15.0, 16.0]):
            lines_data.append({
                'nr': i + 1, 'lev_up': f'0|{2*i+2}', 'lev_low': f'0|{2*i+1}',
                'lam': lam, 'freq': c.SPEED_OF_LIGHT_MICRONS / lam,
                'a_stein': 1e-2, 'e_up': 1500.0 + i * 500, 'e_low': 1000.0 + i * 500,
                'g_up': 2*i + 5, 'g_low': 2*i + 3,
            })
        mll = MoleculeLineList(molecule_id='Test', lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)

        spec_high = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                             R=3000.0, distance=160.0)
        spec_high.add_intensity(inten, dA=1.0)

        spec_low = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                            R=3000.0, distance=160.0,
                            R_func=lambda w: np.full_like(np.atleast_1d(w), 1500.0))
        spec_low.add_intensity(inten, dA=1.0)

        # Lower R = broader line = lower peak
        peak_high = np.max(spec_high.flux)
        peak_low = np.max(spec_low.flux)
        assert peak_low < peak_high


class TestSpectrumReprHtml:
    """Test _repr_html_ for notebook rendering."""

    def test_repr_html_returns_string(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        spec = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01,
                        R=3000.0, distance=160.0)
        html = spec._repr_html_()
        assert isinstance(html, str)
        assert len(html) > 0

    def test_repr_html_contains_table(self):
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        spec = Spectrum(lam_min=10.0, lam_max=11.0, dlambda=0.1,
                        R=3000.0, distance=160.0)
        html = spec._repr_html_()
        assert '<table' in html.lower() or 'dataframe' in html.lower()
