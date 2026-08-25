# -*- coding: utf-8 -*-
"""Regression tests for pixel-sampling bugs.

Bug 1: Changing ``model_pixel_res`` (e.g. via the GUI "Pixel res." box) must
rebuild the molecule's internal Spectrum grid.  Previously the spectrum cache
hash did not include ``model_pixel_res``, so requesting a *finer* resolution
left the model stuck on the old coarse grid (rendered as a histogram-style
staircase after flux-conserving upsampling in ``get_summed_flux``).

Bug 2: "Export Models" (``generate_csv`` / ``generate_all_csv``) ignored the
"Match Pix. Sampling" GUI state and always exported on the model's native
grid.  Exports must now follow ``MoleculeDict.match_spectral_sampling`` by
default, with an explicit ``match_pixel_sampling`` argument to override.
"""

import csv
import os

import numpy as np
import pytest

import iSLAT.Constants as c_mod


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_molecule_with_lines(wavelength_range=(4.9, 5.5), model_pixel_res=0.02,
                              name='PixTestMol'):
    """Create a real Molecule with an injected line list (no file I/O)."""
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

    lines_data = [
        {'nr': i + 1,
         'lev_up': f'0|{2 * i + 2}', 'lev_low': f'0|{2 * i + 1}',
         'lam': lam,
         'freq': c_mod.SPEED_OF_LIGHT_MICRONS / lam,
         'a_stein': 0.02,
         'e_up': 3000.0 - i * 200,
         'e_low': 2000.0 - i * 200,
         'g_up': 2 * i + 3,
         'g_low': 2 * i + 1}
        for i, lam in enumerate([5.0, 5.15, 5.3, 5.45])
    ]
    mll = MoleculeLineList(
        molecule_id=name,
        lines_data=lines_data,
        wavelength_range=wavelength_range,
    )
    mll.partition_function = mll._partition_type(
        t=np.array([100, 300, 500, 1000, 2000], dtype=np.float64),
        q=np.array([10, 150, 500, 2000, 8000], dtype=np.float64),
    )
    mll._molar_mass = 18.015

    mol = Molecule(
        name=name, displaylabel=name,
        filepath=None, color='#0000FF', is_visible=True,
        temp=500.0, radius=1.0, n_mol=1e18, distance=140.0,
        fwhm=130.0,
        wavelength_range=wavelength_range,
        model_pixel_res=model_pixel_res,
        initial_molecule_parameters={
            't_kin': 500.0,
            'scale_exponent': 18.0,
            'scale_number': 1.0,
            'radius_init': 1.0,
        },
    )
    mol.lines = mll
    mol._wavelength_range = wavelength_range
    return mol


def _grid_spacing(lam):
    """Median grid spacing of a wavelength array."""
    return float(np.median(np.diff(lam)))


# Spectrum builds its grid with np.linspace(lam_min, lam_max, int(1 + range/dlambda)),
# so the realized spacing matches the requested dlambda only within a few percent.
_GRID_RTOL = 0.05


# ===========================================================================
# Bug 1 — model_pixel_res must rebuild the spectrum grid
# ===========================================================================

class TestModelPixelResInvalidation:
    """model_pixel_res is spectrum-affecting: the Spectrum is constructed
    with dlambda=model_pixel_res, so changing it must invalidate the
    spectrum cache (dirty flag AND hash)."""

    def test_model_pixel_res_in_spectrum_affecting_params(self):
        from iSLAT.Modules.DataTypes.Molecule import Molecule
        assert 'model_pixel_res' in Molecule.SPECTRUM_AFFECTING_PARAMS

    def test_spectrum_hash_changes_with_model_pixel_res(self):
        mol = _make_molecule_with_lines(model_pixel_res=0.02)
        hash_before = mol._compute_spectrum_hash()
        mol._model_pixel_res = 0.001
        hash_after = mol._compute_spectrum_hash()
        assert hash_before != hash_after, (
            "Spectrum hash must include model_pixel_res, otherwise the "
            "cached spectrum (built with dlambda=model_pixel_res) is never rebuilt"
        )

    def test_setter_marks_spectrum_dirty(self):
        mol = _make_molecule_with_lines(model_pixel_res=0.02)
        mol.get_flux(return_wavelengths=True)  # populate caches
        assert not mol._dirty_flags['spectrum']
        mol.model_pixel_res = 0.001
        assert mol._dirty_flags['spectrum'], (
            "Changing model_pixel_res must mark the spectrum cache dirty"
        )
        assert mol._dirty_flags['flux']

    def test_finer_pixel_res_refines_output_grid(self):
        """The reported bug: increasing pixel sampling (smaller pixel res)
        left the model stuck on the old coarse grid."""
        mol = _make_molecule_with_lines(model_pixel_res=0.02)
        lam_coarse, _ = mol.get_flux(return_wavelengths=True)
        assert _grid_spacing(lam_coarse) == pytest.approx(0.02, rel=_GRID_RTOL)

        mol.model_pixel_res = 0.002  # 10x finer sampling
        lam_fine, flux_fine = mol.get_flux(return_wavelengths=True)

        assert _grid_spacing(lam_fine) == pytest.approx(0.002, rel=_GRID_RTOL), (
            f"Grid spacing stuck at {_grid_spacing(lam_fine):.4g} µm after "
            "requesting 0.002 µm — spectrum was not rebuilt"
        )
        assert len(lam_fine) > len(lam_coarse)
        assert len(flux_fine) == len(lam_fine)

    def test_coarser_pixel_res_coarsens_output_grid(self):
        mol = _make_molecule_with_lines(model_pixel_res=0.002)
        lam_fine, _ = mol.get_flux(return_wavelengths=True)

        mol.model_pixel_res = 0.02
        lam_coarse, _ = mol.get_flux(return_wavelengths=True)

        assert _grid_spacing(lam_coarse) == pytest.approx(0.02, rel=_GRID_RTOL)
        assert len(lam_coarse) < len(lam_fine)

    def test_get_tau_grid_follows_pixel_res(self):
        mol = _make_molecule_with_lines(model_pixel_res=0.02)
        lam_before, _ = mol.get_tau(return_wavelengths=True)
        assert _grid_spacing(lam_before) == pytest.approx(0.02, rel=_GRID_RTOL)

        mol.model_pixel_res = 0.002
        lam_after, _ = mol.get_tau(return_wavelengths=True)
        assert _grid_spacing(lam_after) == pytest.approx(0.002, rel=_GRID_RTOL)


class TestGlobalModelPixelResPropagation:
    """The GUI "Pixel res." box sets MoleculeDict.global_model_pixel_res,
    which must propagate to molecules and refine the summed-flux grid."""

    def _make_dict(self, pixel_res=0.02):
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        md = MoleculeDict(
            global_wavelength_range=(4.9, 5.5),
            global_model_pixel_res=pixel_res,
        )
        mol = _make_molecule_with_lines(model_pixel_res=pixel_res)
        md[mol.name] = mol
        return md

    def test_global_setter_updates_molecule(self):
        md = self._make_dict(pixel_res=0.02)
        md.global_model_pixel_res = 0.002
        assert md['PixTestMol'].model_pixel_res == pytest.approx(0.002)

    def test_summed_flux_grid_refines_after_global_change(self):
        """End-to-end regression of the histogram bug: after increasing the
        pixel sampling globally, the summed-flux canonical grid must be
        genuinely fine-sampled (not a staircase upsampled from a stale
        coarse molecule grid)."""
        md = self._make_dict(pixel_res=0.02)
        wave_obs = np.linspace(4.9, 5.5, 300)

        # Prime caches at coarse resolution
        md.get_summed_flux(wave_obs, visible_only=False)

        md.global_model_pixel_res = 0.002
        lam, flux = md.get_summed_flux(wave_obs, visible_only=False)

        assert _grid_spacing(lam) == pytest.approx(0.002, rel=1e-3)

        # A staircase produced by upsampling a stale 0.02 µm grid would have
        # long runs (~10 pixels) of identical flux values. On a genuinely
        # recomputed fine grid, consecutive samples across a line profile differ.
        line_region = (lam > 4.99) & (lam < 5.01) & (flux > 0)
        if np.count_nonzero(line_region) > 5:
            diffs = np.diff(flux[line_region])
            frac_repeated = np.mean(diffs == 0.0)
            assert frac_repeated < 0.5, (
                "Summed flux looks like a staircase (histogram-style) — "
                "molecule spectrum was not rebuilt at the new pixel resolution"
            )


# ===========================================================================
# Bug 2 — Export Models must honor "Match Pix. Sampling"
# ===========================================================================

class _FakeDataField:
    def insert_text(self, *a, **k):
        pass


def _read_csv_waves(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    return header, np.array([float(r[0]) for r in data])


@pytest.fixture
def export_setup(tmp_path):
    """A MoleculeDict with one real molecule plus an observed grid whose
    sampling differs from the model's native grid."""
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    md = MoleculeDict(
        global_wavelength_range=(4.9, 5.5),
        global_model_pixel_res=0.02,
    )
    mol = _make_molecule_with_lines(model_pixel_res=0.02)
    md[mol.name] = mol
    # Uneven, non-native observed sampling (like real instrument data)
    wave_obs = np.linspace(4.95, 5.45, 173)
    return md, wave_obs, str(tmp_path)


class TestGenerateCsvMatchPixelSampling:

    def test_default_follows_dict_state_off(self, export_setup):
        """match_pixel_sampling=None + state off → native model grid."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = False

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert len(waves) != len(wave_obs)
        assert _grid_spacing(waves) == pytest.approx(0.02, rel=_GRID_RTOL)

    def test_default_follows_dict_state_on(self, export_setup):
        """The reported bug: with 'Match Pix. Sampling' enabled in the GUI,
        the export used the model pixel resolution instead of the data grid."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = True

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert len(waves) == len(wave_obs), (
            "Exported grid must match the observed pixel sampling when the "
            "GUI 'Match Pix. Sampling' state is enabled"
        )
        np.testing.assert_allclose(waves, wave_obs)

    def test_explicit_true_overrides_dict_state(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = False

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, match_pixel_sampling=True)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        np.testing.assert_allclose(waves, wave_obs)

    def test_explicit_false_overrides_dict_state(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = True

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, match_pixel_sampling=False)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert len(waves) != len(wave_obs)
        assert _grid_spacing(waves) == pytest.approx(0.02, rel=_GRID_RTOL)

    def test_matched_export_applies_stellar_rv(self, export_setup):
        """With a nonzero stellar RV, matched exports resample onto the
        rest-frame grid (mirroring get_matched_sampling_wavelengths)."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.global_stellar_rv = 30.0
        md.match_spectral_sampling = True

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        expected = md.apply_stellar_rv(wave_obs)
        np.testing.assert_allclose(waves, expected)

    def test_tau_column_length_matches_matched_grid(self, export_setup):
        """Tau export must use the same grid as flux (previously it always
        used the native grid, which would mismatch a matched flux grid)."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = True

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, save_tau=True)

        header, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert header == ['wave', 'flux', 'tau']
        assert len(waves) == len(wave_obs)

    def test_sum_export_matched(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = False  # override should win

        generate_csv(md, 'SUM', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, match_pixel_sampling=True)

        _, waves = _read_csv_waves(os.path.join(outdir, 'SUM_spec_output.csv'))
        assert len(waves) == len(wave_obs)
        np.testing.assert_allclose(waves, wave_obs)

    def test_sum_export_restores_dict_state(self, export_setup):
        """Temporary override during SUM export must not leak state."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = False

        generate_csv(md, 'SUM', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, match_pixel_sampling=True)

        assert md.match_spectral_sampling is False

    def test_all_export_matched(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.match_spectral_sampling = True

        generate_csv(md, 'ALL', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs)

        for fname in ('PixTestMol_spec_output.csv', 'SUM_spec_output.csv'):
            _, waves = _read_csv_waves(os.path.join(outdir, fname))
            assert len(waves) == len(wave_obs), f"{fname} not on matched grid"

    def test_export_uses_current_pixel_res_after_change(self, export_setup):
        """Combined regression: change pixel res, export unmatched — the
        exported grid must reflect the *new* resolution."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md.get_summed_flux(wave_obs, visible_only=False)  # prime caches

        md.global_model_pixel_res = 0.002
        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, match_pixel_sampling=False)

        _, waves = _read_csv_waves(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert _grid_spacing(waves) == pytest.approx(0.002, rel=_GRID_RTOL)


class TestExportFileNameOption:
    """The Export Models window lets the user choose the output folder and
    file name; generate_csv accepts them via output_dir / file_name."""

    def test_custom_file_name_used(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='my_model.csv')

        assert os.path.exists(os.path.join(outdir, 'my_model.csv'))
        assert not os.path.exists(os.path.join(outdir, 'PixTestMol_spec_output.csv'))

    def test_csv_extension_appended(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='my_model')

        assert os.path.exists(os.path.join(outdir, 'my_model.csv'))

    def test_line_params_follow_custom_stem(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup
        md['PixTestMol'].get_flux(return_wavelengths=True)  # build intensity table

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='my_model.csv')

        assert os.path.exists(os.path.join(outdir, 'my_model_line_params.csv'))

    def test_none_file_name_uses_default(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name=None)

        assert os.path.exists(os.path.join(outdir, 'PixTestMol_spec_output.csv'))

    def test_sum_custom_file_name(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'SUM', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='combined')

        assert os.path.exists(os.path.join(outdir, 'combined.csv'))

    def test_output_dir_is_created(self, export_setup, tmp_path):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, _ = export_setup
        nested = str(tmp_path / 'exports' / 'run1')

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=nested,
                     wave_data=wave_obs, file_name='m.csv')

        assert os.path.exists(os.path.join(nested, 'm.csv'))

    def test_directory_component_in_name_is_stripped(self, export_setup, tmp_path):
        """A path typed into the file-name box must not escape output_dir."""
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'PixTestMol', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='../escaped.csv')

        assert os.path.exists(os.path.join(outdir, 'escaped.csv'))
        assert not os.path.exists(os.path.join(os.path.dirname(outdir), 'escaped.csv'))

    def test_all_ignores_file_name(self, export_setup):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import generate_csv
        md, wave_obs, outdir = export_setup

        generate_csv(md, 'ALL', _FakeDataField(), output_dir=outdir,
                     wave_data=wave_obs, file_name='ignored.csv')

        assert os.path.exists(os.path.join(outdir, 'PixTestMol_spec_output.csv'))
        assert not os.path.exists(os.path.join(outdir, 'ignored.csv'))


class TestResolveOutputPath:
    """Unit tests for the output-path helper."""

    def test_default_when_empty(self, tmp_path):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_output_path
        for value in (None, '', '   '):
            assert _resolve_output_path(str(tmp_path), value, 'default.csv') == \
                os.path.join(str(tmp_path), 'default.csv')

    def test_extension_enforced(self, tmp_path):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_output_path
        assert _resolve_output_path(str(tmp_path), 'name', 'default.csv').endswith('name.csv')
        assert _resolve_output_path(str(tmp_path), 'name.CSV', 'default.csv').endswith('name.CSV')

    def test_basename_only(self, tmp_path):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_output_path
        result = _resolve_output_path(str(tmp_path), '/etc/evil.csv', 'default.csv')
        assert result == os.path.join(str(tmp_path), 'evil.csv')


class TestResolveMatchPixelSampling:
    """Unit tests for the resolver helper."""

    def test_none_reads_state(self):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_match_pixel_sampling

        class _Stub:
            match_spectral_sampling = True

        assert _resolve_match_pixel_sampling(_Stub(), None) is True
        _Stub.match_spectral_sampling = False
        assert _resolve_match_pixel_sampling(_Stub(), None) is False

    def test_explicit_bool_wins(self):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_match_pixel_sampling

        class _Stub:
            match_spectral_sampling = False

        assert _resolve_match_pixel_sampling(_Stub(), True) is True
        _Stub.match_spectral_sampling = True
        assert _resolve_match_pixel_sampling(_Stub(), False) is False

    def test_object_without_state_defaults_false(self):
        from iSLAT.Modules.FileHandling.iSLATFileHandling import _resolve_match_pixel_sampling
        assert _resolve_match_pixel_sampling(object(), None) is False
