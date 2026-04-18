# -*- coding: utf-8 -*-
"""Extended integration tests that exercise multi-module workflows.

These tests go beyond single-class unit tests to validate the interactions
between MoleculeLineList, Intensity, Spectrum, Molecule, MoleculeDict,
and the Plotting subsystem working together.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import iSLAT.Constants as c


# ======================================================================
# Helpers
# ======================================================================

def _make_line_list(molecule_id='H2O', n_lines=5, lam_start=12.0):
    from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
    lines_data = []
    for i in range(n_lines):
        lam = lam_start + i * 1.0
        lines_data.append({
            'nr': i + 1,
            'lev_up': f'0|{2 * i + 2}',
            'lev_low': f'0|{2 * i + 1}',
            'lam': lam,
            'freq': c.SPEED_OF_LIGHT_MICRONS / lam,
            'a_stein': 0.01 + i * 0.005,
            'e_up': 1000 + i * 500,
            'e_low': 500 + i * 500,
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


def _make_intensity(mll, t_kin=500.0, n_mol=1e18):
    from iSLAT.Modules.DataTypes.Intensity import Intensity
    inten = Intensity(mll)
    inten.calc_intensity(t_kin=t_kin, n_mol=n_mol, dv=1.0)
    return inten


def _make_spectrum(inten, lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=160.0):
    from iSLAT.Modules.DataTypes.Spectrum import Spectrum
    spec = Spectrum(
        lam_min=lam_min, lam_max=lam_max, dlambda=dlambda,
        R=R, distance=distance,
    )
    area = np.pi * 1.0 ** 2
    spec.add_intensity(inten, area)
    return spec


def _make_molecule(name='H2O', temp=500.0, n_lines=5, lam_start=12.0):
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    mll = _make_line_list(molecule_id=name, n_lines=n_lines, lam_start=lam_start)
    mol = Molecule(
        name=name,
        displaylabel=name,
        filepath=None,
        color='#FF0000',
        is_visible=True,
        temp=float(temp),
        radius=1.0,
        n_mol=1e18,
        distance=160.0,
        fwhm=130.0,
        initial_molecule_parameters={
            't_kin': float(temp),
            'scale_exponent': 18.0,
            'scale_number': 1.0,
            'radius_init': 1.0,
        },
    )
    mol.lines = mll
    inten = _make_intensity(mll, t_kin=temp)
    mol.intensity = inten
    spec = _make_spectrum(inten)
    mol.spectrum = spec
    return mol


# ======================================================================
# Full pipeline: LineList -> Intensity -> Spectrum -> Plot
# ======================================================================

class TestFullPipeline:
    """End-to-end pipeline from line data to rendered plot."""

    def test_linelist_to_spectrum_plot(self):
        """LineList -> Intensity -> Spectrum -> SpectrumPanel renders."""
        from iSLAT.Modules.Plotting import SpectrumPanel

        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)

        wave = spec.lamgrid
        flux = spec.flux
        assert len(wave) == len(flux)
        assert np.max(flux) > 0

        # Render into a SpectrumPanel
        panel = SpectrumPanel(
            wave, flux, float(wave[0]), float(wave[-1]),
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert len(panel.ax.lines) > 0
        plt.close("all")

    def test_molecule_to_population_diagram(self):
        """Molecule -> PopulationDiagramPlot renders without error."""
        from iSLAT.Modules.Plotting import PopulationDiagramPlot

        mol = _make_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        # Synthetic test data may be too small for valid population diagram
        # but the plot should still be generated
        assert pdp.fig is not None
        assert pdp.ax is not None
        pdp.close()

    def test_molecule_to_full_spectrum_plot(self):
        """Molecule -> FullSpectrumPlot with mol_cache renders."""
        from iSLAT.Modules.Plotting import FullSpectrumPlot

        mol = _make_molecule()
        spec = mol.spectrum
        wave = spec.lamgrid
        flux = spec.flux

        # Build mol_cache the same way the GUI would
        mol_wave = spec.lamgrid
        mol_flux = spec.flux
        cache = [(mol_wave, mol_flux, mol.color, 'H2O', 'h2o')]

        fsp = FullSpectrumPlot(
            wave, flux, n_panels=2, mol_cache=cache,
        )
        fsp.generate_plot()
        assert fsp.fig is not None
        assert len(fsp.panels) == 2
        fsp.close()

    def test_linelist_to_tau_to_optical_depth_panel(self):
        """LineList -> Intensity -> Spectrum -> tau -> OpticalDepthPanel."""
        from iSLAT.Modules.Plotting import OpticalDepthPanel

        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)

        tau = spec.tau_profile
        assert tau is not None and len(tau) > 0

        # Build tau cache as list of (lam, tau, color, label, name) tuples
        tau_cache = [
            (spec.lamgrid, tau, '#FF0000', 'H2O', 'H2O'),
        ]
        panel = OpticalDepthPanel(
            spec.lamgrid, spec.flux,
            float(spec.lamgrid[0]), float(spec.lamgrid[-1]),
            mol_tau_cache=tau_cache,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert panel.ax is not None
        plt.close("all")


# ======================================================================
# Multi-molecule interactions
# ======================================================================

class TestMultiMoleculeWorkflow:
    """Tests combining multiple molecules in a MoleculeDict."""

    def test_molecule_dict_summed_flux(self):
        """MoleculeDict.get_summed_flux from multiple molecules."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

        md = MoleculeDict()
        mol_h2o = _make_molecule(name='H2O', temp=500.0, lam_start=12.0)
        mol_co = _make_molecule(name='CO', temp=300.0, lam_start=13.0)
        md['H2O'] = mol_h2o
        md['CO'] = mol_co

        wave_obs = np.linspace(10.0, 20.0, 500)
        s_wave, s_flux = md.get_summed_flux(wave_obs, visible_only=True)
        assert len(s_wave) == len(s_flux)
        assert np.max(s_flux) > 0

    def test_visibility_toggle_affects_sum(self):
        """Toggling molecule visibility changes the summed flux."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

        md = MoleculeDict()
        mol = _make_molecule(name='H2O', temp=500.0)
        md['H2O'] = mol

        wave_obs = np.linspace(10.0, 20.0, 500)
        _, flux_visible = md.get_summed_flux(wave_obs, visible_only=True)
        peak_visible = np.max(flux_visible) if len(flux_visible) > 0 else 0.0

        mol.is_visible = False
        _, flux_hidden = md.get_summed_flux(wave_obs, visible_only=True)
        peak_hidden = np.max(flux_hidden) if len(flux_hidden) > 0 else 0.0

        # Hidden molecule should produce zero or lower peak flux
        assert peak_visible >= peak_hidden

    def test_multi_molecule_to_spectrum_panel(self):
        """Multiple molecules -> mol_cache -> SpectrumPanel."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        from iSLAT.Modules.Plotting import SpectrumPanel

        md = MoleculeDict()
        colors = ['#FF0000', '#00FF00', '#0000FF']
        names = ['H2O', 'CO', 'OH']
        for i, (name, color) in enumerate(zip(names, colors)):
            mol = _make_molecule(name=name, temp=300 + i * 200, lam_start=12.0 + i)
            mol.color = color
            md[name] = mol

        wave_obs = np.linspace(10.0, 20.0, 500)
        s_wave, s_flux = md.get_summed_flux(wave_obs, visible_only=True)

        cache = []
        for name, mol in md.items():
            if mol.is_visible and mol.spectrum is not None:
                cache.append((
                    mol.spectrum.lamgrid, mol.spectrum.flux,
                    mol.color, name, name.lower(),
                ))

        panel = SpectrumPanel(
            s_wave, s_flux, 10.0, 20.0,
            mol_cache=cache,
            summed_wave=s_wave,
            summed_flux=s_flux,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        mol_lines = [l for l in panel.ax.lines if hasattr(l, "_molecule_name")]
        assert len(mol_lines) == 3
        plt.close("all")


# ======================================================================
# Spectrum + Tau cross-validation
# ======================================================================

class TestSpectrumTauIntegration:
    """Spectrum flux and tau arrays should be self-consistent."""

    def test_flux_and_tau_same_grid(self):
        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)
        assert len(spec.flux) == len(spec.tau_profile)
        assert len(spec.flux) == len(spec.lamgrid)

    def test_tau_positive_where_flux_positive(self):
        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)
        # Where flux is significantly above zero, tau should be > 0
        significant = spec.flux > np.max(spec.flux) * 0.01
        if np.any(significant):
            assert np.all(spec.tau_profile[significant] >= 0)

    def test_per_component_tau_matches_total(self):
        """Adding two intensities to one Spectrum: total tau >= each component."""
        from iSLAT.Modules.DataTypes.Spectrum import Spectrum
        mll1 = _make_line_list(molecule_id='H2O', n_lines=5, lam_start=12.0)
        mll2 = _make_line_list(molecule_id='CO', n_lines=5, lam_start=13.0)
        inten1 = _make_intensity(mll1, t_kin=500.0)
        inten2 = _make_intensity(mll2, t_kin=300.0)

        # Build separate spectra for each component
        spec1 = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=160.0)
        spec2 = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=160.0)
        spec_both = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=160.0)
        area = np.pi * 1.0 ** 2
        spec1.add_intensity(inten1, area)
        spec2.add_intensity(inten2, area)
        spec_both.add_intensity(inten1, area)
        spec_both.add_intensity(inten2, area)

        tau_both = spec_both.tau_profile
        tau_1 = spec1.tau_profile
        tau_2 = spec2.tau_profile

        # Combined tau should be >= each individual (tau is additive)
        assert np.all(tau_both >= tau_1 - 1e-30)
        assert np.all(tau_both >= tau_2 - 1e-30)


# ======================================================================
# Distance scaling pipeline
# ======================================================================

class TestDistanceScalingPipeline:
    """Verify distance scaling propagates correctly through the pipeline."""

    def test_flux_scales_with_distance_squared(self):
        mll = _make_line_list()
        inten = _make_intensity(mll)

        spec_near = _make_spectrum(inten, distance=100.0)
        spec_far = _make_spectrum(inten, distance=200.0)

        peak_near = np.max(spec_near.flux)
        peak_far = np.max(spec_far.flux)

        # Flux ~ 1/d^2, so doubling distance should quarter the flux
        ratio = peak_near / peak_far
        np.testing.assert_allclose(ratio, 4.0, rtol=0.01)

    def test_tau_independent_of_distance(self):
        """Optical depth should NOT change with distance."""
        mll = _make_line_list()
        inten = _make_intensity(mll)

        spec_near = _make_spectrum(inten, distance=100.0)
        spec_far = _make_spectrum(inten, distance=200.0)

        np.testing.assert_allclose(spec_near.tau_profile, spec_far.tau_profile, atol=1e-30)


# ======================================================================
# Cache invalidation cascade
# ======================================================================

class TestCacheInvalidationCascade:
    """Verify that changing parameters invalidates caches properly."""

    def test_intensity_recalc_changes_spectrum(self):
        mll = _make_line_list()
        inten = _make_intensity(mll, t_kin=500.0)
        intens_500 = inten.intensity.copy()

        # Recalculate at different temperature
        inten.calc_intensity(t_kin=800.0, n_mol=1e18, dv=1.0)
        intens_800 = inten.intensity.copy()

        # Different temperature should yield different intensities
        assert not np.allclose(intens_500, intens_800)

    def test_molecule_temp_change_invalidates(self):
        """Changing molecule temperature produces different peak flux."""
        mol = _make_molecule(temp=500.0)
        spec1 = mol.spectrum
        peak1 = np.max(spec1.flux)

        # Build a new molecule at different temp
        mol2 = _make_molecule(temp=800.0)
        peak2 = np.max(mol2.spectrum.flux)

        assert peak1 != peak2

    def test_molecule_dict_cache_refresh(self):
        """MoleculeDict summed flux updates when molecule changes."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

        md = MoleculeDict()
        mol = _make_molecule(name='H2O', temp=500.0)
        md['H2O'] = mol

        wave_obs = np.linspace(10.0, 20.0, 500)
        _, flux1 = md.get_summed_flux(wave_obs, visible_only=True)
        peak1 = np.max(flux1)

        # Replace with hotter molecule => higher peak flux
        mol2 = _make_molecule(name='H2O', temp=900.0)
        md['H2O'] = mol2
        _, flux2 = md.get_summed_flux(wave_obs, visible_only=True)
        peak2 = np.max(flux2)

        assert peak2 != peak1


# ======================================================================
# Intensity batch operations + Spectrum
# ======================================================================

class TestIntensityBatchToSpectrum:
    """Batch intensity methods feeding into the spectrum pipeline."""

    def test_batch_fast_matches_individual(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        mll = _make_line_list(n_lines=5)
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        individual_result = inten.intensity.copy()

        if hasattr(inten, 'calc_intensity_batch_fast'):
            # batch_fast requires arrays + pre-computed context
            t_arr = np.array([500.0])
            n_arr = np.array([1e18])
            dv_arr = np.array([1.0])
            if hasattr(inten, 'prepare_batch_context'):
                ctx = inten.prepare_batch_context()
                batch_result = inten.calc_intensity_batch_fast(
                    t_arr, n_arr, dv_arr, ctx,
                )
                # batch_result shape: (n_walkers, n_lines) => [0] for first walker
                np.testing.assert_allclose(
                    individual_result, batch_result[0], rtol=1e-6,
                )
            else:
                pytest.skip("prepare_batch_context not available")


# ======================================================================
# Line detection + fitting round trip
# ======================================================================

class TestLineDetectionWorkflow:
    """Line detection on synthetic data produces expected results."""

    def test_detect_emission_peaks(self):
        from iSLAT.Modules.DataProcessing.LineAnalyzer import LineAnalyzer

        # Synthetic spectrum with two clear emission features
        wave = np.linspace(10.0, 20.0, 2000)
        flux = np.zeros_like(wave)
        flux += 0.01  # continuum
        # Two Gaussian emission lines
        flux += 0.5 * np.exp(-((wave - 14.0) ** 2) / (2 * 0.02 ** 2))
        flux += 0.3 * np.exp(-((wave - 17.0) ** 2) / (2 * 0.02 ** 2))

        la = LineAnalyzer()
        detected = la.find_single_lines(wave, flux, specsep=0.01, line_threshold=0.1)
        # Should detect at least the two emission lines
        assert len(detected) >= 2

    def test_fitting_engine_gaussian(self):
        """FittingEngine should fit a simple Gaussian."""
        from iSLAT.Modules.DataProcessing.FittingEngine import FittingEngine

        wave = np.linspace(13.5, 14.5, 200)
        center = 14.0
        sigma = 0.03
        amplitude = 0.5
        flux = amplitude * np.exp(-((wave - center) ** 2) / (2 * sigma ** 2))

        fe = FittingEngine()
        result, fitted_wave, fitted_flux = fe.fit_gaussian_line(wave, flux)
        assert result is not None
        # Fitted center should be near the true center
        fit_center = result.params['center'].value
        assert abs(fit_center - center) < 0.01


# ======================================================================
# Stacked plot + residual integration
# ======================================================================

class TestStackedPlotIntegration:
    """Integration between FullSpectrumPlot and ResidualSpectrumPlot."""

    def test_residual_from_model_spectrum(self):
        """Build model flux from molecules, then render residuals."""
        from iSLAT.Modules.Plotting import ResidualSpectrumPlot

        mol = _make_molecule()
        obs_wave = mol.spectrum.lamgrid
        obs_flux = mol.spectrum.flux + np.random.default_rng(42).normal(0, 0.001, len(mol.spectrum.flux))
        model_flux = mol.spectrum.flux

        rsp = ResidualSpectrumPlot(
            obs_wave, obs_flux,
            error_data=np.full_like(obs_flux, 0.001),
            n_panels=2,
            model_flux=model_flux,
        )
        rsp.generate_plot()
        assert rsp.fig is not None
        assert len(rsp.panels) == 2
        rsp.close()

    def test_full_spectrum_and_residual_composite(self):
        """FullSpectrumPlot + ResidualSpectrumPlot -> CompositeStackedPanel."""
        from iSLAT.Modules.Plotting import (
            FullSpectrumPlot,
            ResidualSpectrumPlot,
            CompositeStackedPanel,
        )

        mol = _make_molecule()
        obs_wave = mol.spectrum.lamgrid
        obs_flux = mol.spectrum.flux + np.random.default_rng(42).normal(0, 0.001, len(mol.spectrum.flux))
        model_flux = mol.spectrum.flux

        fsp = FullSpectrumPlot(obs_wave, obs_flux, n_panels=2)
        rsp = ResidualSpectrumPlot(
            obs_wave, obs_flux, n_panels=2, model_flux=model_flux,
        )

        comp = fsp + rsp
        assert isinstance(comp, CompositeStackedPanel)
        comp.generate_plot()
        assert comp.fig is not None
        comp.close()


# ======================================================================
# MainPlotGrid with live molecule data
# ======================================================================

class TestMainPlotGridLive:
    """MainPlotGrid with real molecule/intensity data."""

    def test_full_three_panel_workflow(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        from iSLAT.Modules.Plotting import MainPlotGrid

        mol = _make_molecule()
        md = MoleculeDict()
        md['H2O'] = mol

        # Simulate observed spectrum from model
        obs_wave = mol.spectrum.lamgrid
        obs_flux = mol.spectrum.flux + np.random.default_rng(42).normal(0, 0.001, len(mol.spectrum.flux))

        mpg = MainPlotGrid(
            obs_wave, obs_flux,
            molecules=md,
            active_molecule=mol,
            spectrum_range=(11.0, 19.0),
            inspection_range=(13.0, 15.0),
        )
        mpg.generate_plot()

        # All three panels should have content
        assert len(mpg.ax_spectrum.lines) > 0
        assert len(mpg.ax_inspection.lines) > 0
        # Population diagram should have been rendered (may show
        # "No valid data" for synthetic test data)
        assert mpg.ax_popdiagram is not None
        mpg.close()


# ======================================================================
# Chi2Spectrum integration
# ======================================================================

class TestChi2Integration:
    """Chi2Spectrum cross-module integration."""

    def test_chi2_from_model_spectrum(self):
        from iSLAT.Modules.DataProcessing.chi2spectrum import Chi2Spectrum, FluxMeasurement

        mol = _make_molecule()
        obs_wave = mol.spectrum.lamgrid
        model_flux = mol.spectrum.flux

        chi2 = Chi2Spectrum()
        # Add a measurement in a window where the model has flux
        peak_idx = np.argmax(model_flux)
        peak_lam = obs_wave[peak_idx]
        fm = FluxMeasurement(
            lam_min=peak_lam - 0.05,
            lam_max=peak_lam + 0.05,
            flux=float(model_flux[peak_idx]),
            flux_error=float(model_flux[peak_idx]) * 0.1,
        )
        chi2.add_measurement(fm)
        assert len(chi2.measurements) == 1


# ======================================================================
# Spectral resampling in pipeline context
# ======================================================================

class TestResamplingInPipeline:
    """Spectral resampling preserves features through the pipeline."""

    def test_resample_preserves_peak(self):
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres

        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)

        # Resample to coarser grid
        new_wave = np.linspace(11.0, 19.0, 200)
        resampled = spectres(new_wave, spec.lamgrid, spec.flux)

        # Peak should be roughly preserved
        original_peak = np.max(spec.flux)
        resampled_peak = np.max(resampled)
        assert resampled_peak > 0
        # Coarse resampling loses peak height; just verify within an order of magnitude
        assert resampled_peak / original_peak > 0.1

    def test_identity_resample_within_pipeline(self):
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres

        mll = _make_line_list()
        inten = _make_intensity(mll)
        spec = _make_spectrum(inten)

        # Resample to same grid should be identity
        resampled = spectres(spec.lamgrid, spec.lamgrid, spec.flux)
        np.testing.assert_allclose(resampled, spec.flux, atol=1e-10)


# ======================================================================
# Optical depth full-stack rendering
# ======================================================================

class TestOpticalDepthFullStack:
    """Full stack: build molecule -> compute tau -> render OpticalDepthSpectrumPlot."""

    def test_od_spectrum_plot(self):
        from iSLAT.Modules.Plotting import OpticalDepthSpectrumPlot
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict

        mol = _make_molecule()
        md = MoleculeDict()
        md['H2O'] = mol
        spec = mol.spectrum
        odsp = OpticalDepthSpectrumPlot(
            spec.lamgrid, spec.flux,
            n_panels=2,
            molecules=md,
        )
        odsp.generate_plot()
        assert odsp.fig is not None
        assert len(odsp.panels) == 2
        odsp.close()

    def test_od_with_molecule_object(self):
        from iSLAT.Modules.Plotting import OpticalDepthPanel

        mol = _make_molecule()
        panel = OpticalDepthPanel(
            mol.spectrum.lamgrid,
            mol.spectrum.flux,
            float(mol.spectrum.lamgrid[0]),
            float(mol.spectrum.lamgrid[-1]),
            molecule=mol,
            ax=plt.subplots()[1],
        )
        panel.generate_plot()
        assert panel.ax is not None
        plt.close("all")


# ======================================================================
# MoleculeLineList + Intensity build_table round-trip
# ======================================================================

class TestBuildTableRoundTrip:
    """Intensity.build_table feeds correctly into PopulationDiagramPlot._build_table_from_intensity."""

    def test_build_table_has_required_columns(self):
        inten = _make_intensity(_make_line_list())
        if hasattr(inten, 'build_table'):
            df = inten.build_table(full_range=True)
            assert df is not None
            for col in ['lam', 'intens', 'a_stein', 'g_up', 'e_up']:
                assert col in df.columns

    def test_population_diagram_from_build_table(self):
        from iSLAT.Modules.Plotting import PopulationDiagramPlot

        inten = _make_intensity(_make_line_list())
        pdp = PopulationDiagramPlot(
            intensity=inten, name="TestMol",
            color="blue", radius=1.0, distance=160.0,
        )
        pdp.generate_plot()
        # With synthetic data the intensity values may be too small
        # to produce valid population diagram points, so just check
        # that the plot was generated without error
        assert pdp.fig is not None
        assert pdp.ax is not None
        pdp.close()


# ======================================================================
# Molecule parameter hashing + caching correctness
# ======================================================================

class TestParameterHashCaching:
    """Molecule parameter hash should change when parameters change."""

    def test_different_temps_different_hashes(self):
        mol1 = _make_molecule(temp=500.0)
        mol2 = _make_molecule(temp=800.0)
        if hasattr(mol1, 'get_parameter_hash') and hasattr(mol2, 'get_parameter_hash'):
            assert mol1.get_parameter_hash() != mol2.get_parameter_hash()

    def test_same_params_same_hash(self):
        mol1 = _make_molecule(temp=500.0)
        mol2 = _make_molecule(temp=500.0)
        if hasattr(mol1, 'get_parameter_hash') and hasattr(mol2, 'get_parameter_hash'):
            assert mol1.get_parameter_hash() == mol2.get_parameter_hash()
