# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.subtract_models_from_data.

Covers:
- The import fix (_get_spectral_utils, not _spectres directly)
- Core spectres resampling correctness
- subtract_models_from_data happy path (produces correct subtracted flux)
- subtract_models_from_data with no visible molecules
- subtract_models_from_data with missing wave/flux data
- subtract_models_from_data output file is written
- subtract_models_from_data with multiple visible molecules
"""
import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_wave(start=10.0, stop=20.0, n=200):
    return np.linspace(start, stop, n)


def _make_flux(wave, amp=1.0):
    """Simple Gaussian bump as a fake 'molecule' flux."""
    center = (wave[0] + wave[-1]) / 2
    return amp * np.exp(-0.5 * ((wave - center) / 1.0) ** 2)


# ---------------------------------------------------------------------------
# Import / lazy-loader tests
# ---------------------------------------------------------------------------

class TestSpectralUtilsImport:
    """The _spectres name no longer exists at module level; _get_spectral_utils must be used."""

    def test_spectres_not_importable_as_module_attribute(self):
        from iSLAT.Modules.DataTypes import Molecule as _mol_module
        assert not hasattr(_mol_module, '_spectres'), (
            "_spectres should NOT be a top-level name in Molecule.py; "
            "use _get_spectral_utils() instead."
        )

    def test_get_spectral_utils_returns_tuple(self):
        from iSLAT.Modules.DataTypes.Molecule import _get_spectral_utils
        result = _get_spectral_utils()
        assert isinstance(result, tuple) and len(result) == 2

    def test_get_spectral_utils_returns_callables(self):
        from iSLAT.Modules.DataTypes.Molecule import _get_spectral_utils
        make_bins, spectres = _get_spectral_utils()
        assert callable(make_bins)
        assert callable(spectres)

    def test_spectres_via_lazy_loader_runs(self):
        from iSLAT.Modules.DataTypes.Molecule import _get_spectral_utils
        _, spectres = _get_spectral_utils()
        old_w = np.linspace(10, 20, 100)
        old_f = np.ones(100)
        new_w = np.linspace(11, 19, 50)
        result = spectres(new_w, old_w, old_f, fill=0.0)
        assert result.shape == new_w.shape
        np.testing.assert_allclose(result, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# spectres correctness tests
# ---------------------------------------------------------------------------

class TestSpectresResampling:
    """Direct tests of the spectres resampling function."""

    @pytest.fixture(autouse=True)
    def _get_spectres(self):
        from iSLAT.Modules.DataProcessing.spectral_utils import spectres
        self.spectres = spectres

    def test_identity_resampling_flat_spectrum(self):
        """Resampling a flat spectrum onto the same grid gives the same values."""
        w = np.linspace(10, 20, 100)
        f = np.ones(100)
        out = self.spectres(w, w, f, fill=0.0)
        np.testing.assert_allclose(out, f, atol=1e-10)

    def test_coarser_grid_preserves_flux_density(self):
        """Coarser output grid should give ~1.0 for a flat input."""
        old_w = np.linspace(10, 20, 200)
        old_f = np.ones(200)
        new_w = np.linspace(10.1, 19.9, 50)
        out = self.spectres(new_w, old_w, old_f, fill=0.0)
        np.testing.assert_allclose(out, 1.0, atol=1e-6)

    def test_finer_grid_preserves_flux_density(self):
        """Finer output grid should still give ~1.0 for a flat input."""
        old_w = np.linspace(10, 20, 50)
        old_f = np.ones(50)
        new_w = np.linspace(10.2, 19.8, 200)
        out = self.spectres(new_w, old_w, old_f, fill=0.0)
        np.testing.assert_allclose(out, 1.0, atol=1e-6)

    def test_fill_value_used_outside_range(self):
        """Bins outside the source range receive the fill value."""
        old_w = np.linspace(12, 18, 60)
        old_f = np.ones(60)
        new_w = np.linspace(10, 20, 100)
        out = self.spectres(new_w, old_w, old_f, fill=-99.0)
        # Points outside [12, 18] must be fill
        outside = (new_w < old_w[0]) | (new_w > old_w[-1])
        assert np.all(out[outside] == -99.0), "Fill value not applied outside source range"

    def test_output_shape_matches_new_wavs(self):
        old_w = np.linspace(10, 20, 100)
        old_f = np.ones(100)
        for n in (10, 50, 150, 300):
            new_w = np.linspace(11, 19, n)
            out = self.spectres(new_w, old_w, old_f)
            assert out.shape == (n,)

    def test_zero_flux_returns_zero(self):
        w = np.linspace(10, 20, 100)
        new_w = np.linspace(11, 19, 50)
        out = self.spectres(new_w, w, np.zeros(100), fill=0.0)
        np.testing.assert_allclose(out, 0.0, atol=1e-15)


# ---------------------------------------------------------------------------
# subtract_models_from_data integration tests (mocked iSLAT object)
# ---------------------------------------------------------------------------

def _build_islat_mock(wave, flux, mol_fluxes, visible_names=None, tmpdir=None):
    """
    Build a minimal mock of the iSLAT class that exercises
    subtract_models_from_data.

    Parameters
    ----------
    wave : np.ndarray   — data wavelength grid
    flux : np.ndarray   — data flux
    mol_fluxes : dict   — {mol_name: flux_array} on the *same* grid as wave
    visible_names : list|None — which molecules are "visible" (default: all)
    tmpdir : str|None   — directory for saving output (default: system tmp)
    """
    from iSLAT.Modules.DataTypes.Molecule import _get_spectral_utils

    if visible_names is None:
        visible_names = list(mol_fluxes.keys())

    # Build molecule mocks
    molecules_dict = {}
    for name, mflux in mol_fluxes.items():
        mol = MagicMock()
        mol.name = name
        mol.get_flux.return_value = (wave.copy(), mflux.copy())
        mol.is_visible = name in visible_names
        molecules_dict[name] = mol

    mol_dict_mock = MagicMock()
    mol_dict_mock.__contains__ = lambda self, k: k in molecules_dict
    mol_dict_mock.__getitem__ = lambda self, k: molecules_dict[k]
    mol_dict_mock.get_visible_molecules.return_value = visible_names
    mol_dict_mock.keys.return_value = list(mol_fluxes.keys())
    # Make bool(mol_dict_mock) True
    mol_dict_mock.__bool__ = lambda self: True

    islat = MagicMock()
    islat.wave_data = wave
    islat.flux_data = flux
    islat.err_data = np.ones_like(flux) * 0.01
    islat.continuum_data = None
    islat.molecules_dict = mol_dict_mock
    islat.GUI = None

    # tmpdir for output file
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    fake_spectrum = os.path.join(tmpdir, "spectrum.csv")
    # Create the fake file so Path(...).parent exists
    Path(fake_spectrum).write_text("wave,flux\n")
    islat.loaded_spectrum_file = fake_spectrum
    islat.load_spectrum = MagicMock()

    # Bind the real method to the mock
    from iSLAT.iSLATClass import iSLAT as iSLATClass
    islat.subtract_models_from_data = lambda visible_only=True: (
        iSLATClass.subtract_models_from_data(islat, visible_only=visible_only)
    )

    return islat, tmpdir


class TestSubtractModelsFromData:

    def test_import_fix_no_error(self):
        """Calling subtract_models_from_data must not raise ImportError."""
        wave = _make_wave(10, 20, 200)
        flux = np.ones(200)
        mol_flux = _make_flux(wave, amp=0.3)
        islat, _ = _build_islat_mock(wave, flux, {"H2O": mol_flux})
        # Should not raise
        result = islat.subtract_models_from_data(visible_only=True)
        assert result is not None

    def test_subtracted_flux_correct_single_molecule(self):
        """Output flux = data flux − model flux (resampled, but same grid here)."""
        wave = _make_wave(10, 20, 200)
        data_flux = 2.0 * np.ones(200)
        mol_flux = _make_flux(wave, amp=1.0)

        islat, tmpdir = _build_islat_mock(wave, data_flux, {"H2O": mol_flux})
        result = islat.subtract_models_from_data(visible_only=True)
        assert result is not None

        out_df = __import__('pandas').read_csv(result)
        expected = data_flux - mol_flux
        np.testing.assert_allclose(out_df['flux'].values, expected, atol=1e-6)

    def test_subtracted_flux_correct_two_molecules(self):
        """Subtraction of two molecules sums their fluxes first."""
        wave = _make_wave(10, 20, 200)
        data_flux = 3.0 * np.ones(200)
        mol1 = _make_flux(wave, amp=1.0)
        mol2 = _make_flux(wave, amp=0.5)

        islat, tmpdir = _build_islat_mock(
            wave, data_flux, {"H2O": mol1, "CO2": mol2}
        )
        result = islat.subtract_models_from_data(visible_only=True)
        assert result is not None

        out_df = __import__('pandas').read_csv(result)
        expected = data_flux - mol1 - mol2
        np.testing.assert_allclose(out_df['flux'].values, expected, atol=1e-6)

    def test_visible_only_true_skips_hidden_molecules(self):
        """With visible_only=True, hidden molecules must not be subtracted."""
        wave = _make_wave(10, 20, 200)
        data_flux = 2.0 * np.ones(200)
        mol_vis = _make_flux(wave, amp=0.5)
        mol_hid = _make_flux(wave, amp=1.0)

        islat, tmpdir = _build_islat_mock(
            wave, data_flux,
            {"H2O": mol_vis, "CO": mol_hid},
            visible_names=["H2O"],
        )
        result = islat.subtract_models_from_data(visible_only=True)
        assert result is not None

        out_df = __import__('pandas').read_csv(result)
        # Only mol_vis should have been subtracted
        expected = data_flux - mol_vis
        np.testing.assert_allclose(out_df['flux'].values, expected, atol=1e-6)

    def test_visible_only_false_includes_all(self):
        """With visible_only=False, all molecules (even hidden) are subtracted."""
        wave = _make_wave(10, 20, 200)
        data_flux = 3.0 * np.ones(200)
        mol1 = _make_flux(wave, amp=0.5)
        mol2 = _make_flux(wave, amp=0.8)

        islat, tmpdir = _build_islat_mock(
            wave, data_flux,
            {"H2O": mol1, "CO": mol2},
            visible_names=["H2O"],   # CO hidden
        )
        # Override: visible_only=False should include CO too
        # We need mol_dict_mock.keys() to return both
        result = islat.subtract_models_from_data(visible_only=False)
        assert result is not None

        out_df = __import__('pandas').read_csv(result)
        expected = data_flux - mol1 - mol2
        np.testing.assert_allclose(out_df['flux'].values, expected, atol=1e-6)

    def test_output_file_has_correct_suffix(self):
        """The output filename must contain '_iSLATsub'."""
        wave = _make_wave(10, 20, 100)
        islat, _ = _build_islat_mock(wave, np.ones(100), {"H2O": np.zeros(100)})
        result = islat.subtract_models_from_data()
        assert result is not None
        assert "_iSLATsub" in Path(result).name

    def test_output_file_contains_wave_and_flux_columns(self):
        """CSV output must have at least 'wave' and 'flux' columns."""
        import pandas as pd
        wave = _make_wave(10, 20, 100)
        islat, _ = _build_islat_mock(wave, np.ones(100), {"H2O": np.zeros(100)})
        result = islat.subtract_models_from_data()
        df = pd.read_csv(result)
        assert 'wave' in df.columns
        assert 'flux' in df.columns

    def test_output_file_includes_err_column_when_available(self):
        """If err_data is set, the output CSV must include an 'err' column."""
        import pandas as pd
        wave = _make_wave(10, 20, 100)
        islat, _ = _build_islat_mock(wave, np.ones(100), {"H2O": np.zeros(100)})
        islat.err_data = np.full(100, 0.05)
        result = islat.subtract_models_from_data()
        df = pd.read_csv(result)
        assert 'err' in df.columns

    def test_returns_none_when_no_wave_data(self):
        """Returns None (and does not crash) when wave_data is missing."""
        islat = MagicMock()
        del islat.wave_data  # attribute missing
        islat.flux_data = np.ones(100)
        islat.molecules_dict = MagicMock()
        islat.molecules_dict.__bool__ = lambda s: True
        from iSLAT.iSLATClass import iSLAT as iSLATClass
        result = iSLATClass.subtract_models_from_data(islat)
        assert result is None

    def test_returns_none_when_no_molecules(self):
        """Returns None when molecules_dict is empty."""
        wave = _make_wave(10, 20, 100)
        islat = MagicMock()
        islat.wave_data = wave
        islat.flux_data = np.ones(100)
        islat.molecules_dict = MagicMock()
        islat.molecules_dict.__bool__ = lambda s: False
        from iSLAT.iSLATClass import iSLAT as iSLATClass
        result = iSLATClass.subtract_models_from_data(islat)
        assert result is None

    def test_returns_none_when_no_visible_molecules(self):
        """Returns None when all molecules are hidden (visible_only=True)."""
        wave = _make_wave(10, 20, 100)
        islat, _ = _build_islat_mock(
            wave, np.ones(100), {"H2O": np.zeros(100)},
            visible_names=[],  # nothing visible
        )
        result = islat.subtract_models_from_data(visible_only=True)
        assert result is None

    def test_load_spectrum_called_after_subtraction(self):
        """subtract_models_from_data must call load_spectrum on the output path."""
        wave = _make_wave(10, 20, 100)
        islat, _ = _build_islat_mock(wave, np.ones(100), {"H2O": np.zeros(100)})
        result = islat.subtract_models_from_data()
        assert result is not None
        islat.load_spectrum.assert_called_once()
        call_arg = islat.load_spectrum.call_args[0][0]
        assert "_iSLATsub" in call_arg

    def test_subtraction_does_not_modify_original_flux(self):
        """The original islat.flux_data must not be mutated."""
        wave = _make_wave(10, 20, 200)
        original_flux = np.ones(200)
        islat, _ = _build_islat_mock(
            wave, original_flux.copy(), {"H2O": _make_flux(wave, amp=0.5)}
        )
        islat.flux_data = original_flux.copy()
        _ = islat.subtract_models_from_data()
        np.testing.assert_array_equal(islat.flux_data, original_flux)
