# -*- coding: utf-8 -*-
"""Tests for PopulationDiagramPlot axis-change, color-by and active-line
scatter behaviour.

Covers:
- render_active_lines plots scatter at the correct position for the
  default axes (eu / rd_yax).
- After set_axes(), the scatter is re-rendered at coordinates that match
  the new axis properties.
- After color_by(), the scatter is re-rendered and present on the axes.
- clear_active_lines() removes scatter artists and clears the cache so
  generate_plot() does not resurrect them.
- Multiple set_axes() / color_by() calls keep the scatter consistent.
- render_active_lines returns None when called without a molecule.
- Scatter coordinates fall within the axes data limits after axis change.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import iSLAT.Constants as c
from iSLAT.Modules.Plotting import PopulationDiagramPlot


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_line_list(n=5, lam_start=12.0):
    from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
    lines = []
    for i in range(n):
        lam = lam_start + i * 1.5
        lines.append({
            "nr": i + 1,
            "lev_up": f"0|{2*i+2}",
            "lev_low": f"0|{2*i+1}",
            "lam": lam,
            "freq": c.SPEED_OF_LIGHT_MICRONS / lam,
            "a_stein": 0.01 + i * 0.005,
            "e_up": 1000 + i * 600,
            "e_low": 500 + i * 600,
            "g_up": 2 * i + 3,
            "g_low": 2 * i + 1,
        })
    mll = MoleculeLineList(molecule_id="TestMol", lines_data=lines)
    mll.partition_function = mll._partition_type(
        t=np.array([100, 300, 500, 1000, 2000], dtype=float),
        q=np.array([10, 150, 500, 2000, 8000], dtype=float),
    )
    mll._molar_mass = 18.015
    return mll


def _make_intensity(n=5, lam_start=12.0, t_kin=500.0, n_mol=1e18):
    from iSLAT.Modules.DataTypes.Intensity import Intensity
    mll = _make_line_list(n=n, lam_start=lam_start)
    inten = Intensity(mll)
    inten.calc_intensity(t_kin=t_kin, n_mol=n_mol, dv=1.0)
    return inten


def _make_molecule(name="H2O", temp=500.0, n=5, lam_start=12.0):
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.Spectrum import Spectrum
    mll = _make_line_list(n=n, lam_start=lam_start)
    mol = Molecule(
        name=name,
        displaylabel=name,
        filepath=None,
        color="#3399FF",
        is_visible=True,
        temp=float(temp),
        radius=1.0,
        n_mol=1e18,
        distance=160.0,
        fwhm=130.0,
        initial_molecule_parameters={
            "t_kin": float(temp),
            "scale_exponent": 18.0,
            "scale_number": 1.0,
            "radius_init": 1.0,
        },
    )
    mol.lines = mll
    mol.intensity = _make_intensity(n=n, lam_start=lam_start, t_kin=temp)
    spec = Spectrum(lam_min=10.0, lam_max=30.0, dlambda=0.01, R=3000.0, distance=160.0)
    area = np.pi * 1.0 ** 2
    spec.add_intensity(mol.intensity, area)
    mol.spectrum = spec
    return mol


def _line_data_from_molecule(mol):
    """Return (MoleculeLine, intensity, tau) triples from molecule intensity."""
    inten = mol.intensity
    table = inten.build_table()
    if table is None or table.empty:
        return []

    # Ensure MoleculeLine objects have been created lazily
    mll = mol.lines  # MoleculeLineList
    mll._ensure_lines_created()
    molecule_lines = mll.lines or []

    result = []
    for line in molecule_lines:
        if line.lam is None:
            continue
        mask = table["lam"].apply(lambda v: abs(v - line.lam) < 0.001)
        row = table[mask]
        if row.empty:
            continue
        intensity_val = float(row.iloc[0].get("intens", 0.0))
        tau_val = float(row.iloc[0].get("tau", 0.0)) if "tau" in row.columns else 0.0
        result.append((line, intensity_val, tau_val))
    return result


def _scatter_offsets(pdp):
    """Return all scatter offsets currently drawn on pdp's axes."""
    offsets = []
    for coll in pdp.ax.collections:
        try:
            off = coll.get_offsets()
            if len(off):
                offsets.append(np.asarray(off))
        except Exception:
            pass
    return offsets


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mol():
    return _make_molecule()


@pytest.fixture
def line_data(mol):
    return _line_data_from_molecule(mol)


@pytest.fixture
def pdp_with_lines(mol, line_data):
    """A PopulationDiagramPlot that has had render_active_lines called."""
    pdp = PopulationDiagramPlot(molecule=mol)
    pdp.generate_plot()
    active_lines = []
    pdp.render_active_lines(line_data, active_lines, molecule=mol, color="green")
    yield pdp, active_lines
    pdp.close()


# ===========================================================================
# 1. render_active_lines — basic contract
# ===========================================================================

class TestRenderActiveLines:

    def test_returns_scatter_collection(self, mol, line_data):
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        sc = pdp.render_active_lines(line_data, active_lines, molecule=mol, color="green")
        assert sc is not None, "Expected a PathCollection to be returned"
        pdp.close()

    def test_active_lines_populated(self, mol, line_data):
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)
        assert len(active_lines) > 0, "active_lines should have been populated"
        pdp.close()

    def test_scatter_present_on_axes(self, mol, line_data):
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol, color="green")
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), "No scatter points found on axes"
        pdp.close()

    def test_returns_none_without_molecule(self, line_data):
        pdp = PopulationDiagramPlot()
        pdp.generate_plot()
        active_lines = []
        result = pdp.render_active_lines(line_data, active_lines, molecule=None)
        assert result is None
        pdp.close()

    def test_returns_none_with_empty_line_data(self, mol):
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        result = pdp.render_active_lines([], active_lines, molecule=mol)
        assert result is None
        pdp.close()

    def test_cache_populated_after_render(self, mol, line_data):
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)
        assert pdp._active_lines_cache is not None
        assert pdp._active_lines_cache["molecule"] is mol
        pdp.close()

    def test_scatter_x_coords_match_e_up(self, mol, line_data):
        """Default axes are eu (x) / rd_yax (y) — x values should be e_up values."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        expected_e_ups = sorted(
            line.e_up for line, _, _ in line_data
            if line.e_up is not None
        )
        all_offsets = np.concatenate(_scatter_offsets(pdp))
        actual_x = sorted(all_offsets[:, 0].tolist())

        # The scatter x values should be a subset of e_up values
        for xv in actual_x:
            assert any(abs(xv - eu) < 1.0 for eu in expected_e_ups), (
                f"Scatter x={xv} not found in e_up values {expected_e_ups}"
            )
        pdp.close()


# ===========================================================================
# 2. clear_active_lines
# ===========================================================================

class TestClearActiveLines:

    def test_clears_active_lines_list(self, pdp_with_lines):
        pdp, active_lines = pdp_with_lines
        assert len(active_lines) > 0, "Pre-condition: active_lines must be populated"
        pdp.clear_active_lines(active_lines)
        assert active_lines == []

    def test_clears_cache(self, pdp_with_lines):
        pdp, active_lines = pdp_with_lines
        assert pdp._active_lines_cache is not None, "Pre-condition: cache must be set"
        pdp.clear_active_lines(active_lines)
        assert pdp._active_lines_cache is None

    def test_generate_plot_does_not_restore_after_clear(self, mol, line_data):
        """After clear_active_lines, generate_plot must NOT re-render active-line scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()

        # Record scatter offset count from the base plot only (no active lines)
        n_base = sum(len(c.get_offsets()) for c in pdp.ax.collections)

        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        # After rendering active lines there should be more scatter
        n_before = sum(len(c.get_offsets()) for c in pdp.ax.collections)
        assert n_before > n_base

        pdp.clear_active_lines(active_lines)
        pdp.generate_plot()  # should NOT re-add active-line scatter

        n_after = sum(len(c.get_offsets()) for c in pdp.ax.collections)
        assert n_after == n_base, (
            f"Expected {n_base} scatter offsets after clear+regenerate (base only), got {n_after}"
        )
        pdp.close()


# ===========================================================================
# 3. set_axes — scatter persists and updates coordinates
# ===========================================================================

class TestSetAxesScatterConsistency:

    def test_scatter_persists_after_set_axes(self, mol, line_data):
        """set_axes must re-render the active-line scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        # Change to wavelength (x) vs intensity (y)
        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)

        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), (
            "Active-line scatter disappeared after set_axes()"
        )
        pdp.close()

    def test_scatter_x_coords_change_with_x_prop(self, mol, line_data):
        """After changing x_prop to 'wavelength', scatter x values are wavelengths."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="wavelength", y_prop="rd_yax", regenerate=True)

        expected_lams = sorted(line.lam for line, _, _ in line_data if line.lam is not None)
        all_offsets = np.concatenate(_scatter_offsets(pdp))
        actual_x = sorted(all_offsets[:, 0].tolist())

        for xv in actual_x:
            assert any(abs(xv - lam) < 0.1 for lam in expected_lams), (
                f"Scatter x={xv:.3f} not in expected wavelengths {expected_lams}"
            )
        pdp.close()

    def test_scatter_x_coords_change_with_e_low(self, mol, line_data):
        """After x_prop='e_low', scatter x values are lower-level energies."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="e_low", y_prop="rd_yax", regenerate=True)

        expected_e_lows = sorted(line.e_low for line, _, _ in line_data if line.e_low is not None)
        all_offsets = np.concatenate(_scatter_offsets(pdp))
        actual_x = sorted(all_offsets[:, 0].tolist())

        for xv in actual_x:
            assert any(abs(xv - el) < 1.0 for el in expected_e_lows), (
                f"Scatter x={xv} not in expected e_low values {expected_e_lows}"
            )
        pdp.close()

    def test_scatter_y_coords_change_with_y_prop(self, mol, line_data):
        """After y_prop='intens', scatter y values are raw intensities."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="eu", y_prop="intens", regenerate=True)

        expected_intensities = [intensity for _, intensity, _ in line_data if intensity > 0]
        all_offsets = np.concatenate(_scatter_offsets(pdp))
        actual_y = sorted(all_offsets[:, 1].tolist())

        for yv in actual_y:
            assert any(abs(yv - iv) / max(iv, 1e-30) < 0.01 for iv in expected_intensities), (
                f"Scatter y={yv} not within 1% of any intensity value"
            )
        pdp.close()

    def test_multiple_axis_changes_scatter_consistent(self, mol, line_data):
        """Multiple sequential set_axes() calls all preserve scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        for x_prop, y_prop in [
            ("eu", "rd_yax"),
            ("wavelength", "intens"),
            ("e_low", "rd_yax"),
            ("eu", "rd_yax"),
        ]:
            pdp.set_axes(x_prop=x_prop, y_prop=y_prop, regenerate=True)
            offsets = _scatter_offsets(pdp)
            assert any(len(o) > 0 for o in offsets), (
                f"Scatter gone after set_axes(x_prop={x_prop!r}, y_prop={y_prop!r})"
            )
        pdp.close()

    def test_log_scale_axes_scatter_present(self, mol, line_data):
        """Log scale axes: scatter should still be plotted."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="eu", y_prop="intens", x_log=True, y_log=True, regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), "Scatter gone after log scale axis change"
        pdp.close()

    def test_reset_axes_scatter_present(self, mol, line_data):
        """Resetting axes back to defaults also restores scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        pdp.set_axes()  # reset defaults
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), "Scatter gone after reset to default axes"
        pdp.close()

    def test_active_lines_list_refreshed_after_set_axes(self, mol, line_data):
        """active_lines list should still hold entries after set_axes."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)
        n_before = len(active_lines)
        assert n_before > 0

        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        assert len(active_lines) == n_before, (
            f"active_lines count changed from {n_before} to {len(active_lines)}"
        )
        pdp.close()


# ===========================================================================
# 4. color_by — scatter persists after color remapping
# ===========================================================================

class TestColorByScatterConsistency:

    def test_scatter_persists_after_color_by(self, mol, line_data):
        """color_by() must not remove the active-line scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.color_by("e_up", regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), (
            "Active-line scatter disappeared after color_by()"
        )
        pdp.close()

    def test_scatter_persists_after_categorical_color_by(self, mol, line_data):
        """Categorical color_by (e.g., molecule) must also keep scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.color_by("molecule", regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), (
            "Active-line scatter disappeared after categorical color_by('molecule')"
        )
        pdp.close()

    def test_scatter_persists_after_clear_color_mapping(self, mol, line_data):
        """clear_color_mapping() must also preserve the active-line scatter."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.color_by("e_up", regenerate=True)
        pdp.clear_color_mapping(regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), (
            "Active-line scatter disappeared after clear_color_mapping()"
        )
        pdp.close()

    def test_color_by_then_set_axes_scatter_present(self, mol, line_data):
        """color_by() followed by set_axes() must keep scatter in new coords."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.color_by("a_stein", regenerate=True)
        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), (
            "Scatter gone after color_by() + set_axes()"
        )
        pdp.close()

    def test_scatter_x_coords_unchanged_by_color_by(self, mol, line_data):
        """color_by() changes colour but must not shift scatter positions."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        # Record scatter positions before color_by
        before_offsets = np.concatenate(_scatter_offsets(pdp))

        pdp.color_by("e_up", regenerate=True)
        after_offsets = np.concatenate(_scatter_offsets(pdp))

        # x positions should be unchanged (still eu)
        np.testing.assert_allclose(
            sorted(before_offsets[:, 0]),
            sorted(after_offsets[:, 0]),
            rtol=1e-6,
            err_msg="Scatter x positions changed after color_by()",
        )
        pdp.close()

    def test_multiple_color_by_calls(self, mol, line_data):
        """Repeated color_by() calls must all keep scatter present."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        for prop in ["e_up", "a_stein", "g_up", "wavelength"]:
            pdp.color_by(prop, regenerate=True)
            offsets = _scatter_offsets(pdp)
            assert any(len(o) > 0 for o in offsets), (
                f"Scatter gone after color_by({prop!r})"
            )
        pdp.close()


# ===========================================================================
# 5. Combined axis + color interaction
# ===========================================================================

class TestCombinedAxisColorInteraction:

    def test_set_axes_before_render_then_render(self, mol, line_data):
        """If axes are set BEFORE render_active_lines, coords should match."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=False)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        expected_lams = sorted(line.lam for line, _, _ in line_data if line.lam)
        all_offsets = np.concatenate(_scatter_offsets(pdp))
        actual_x = sorted(all_offsets[:, 0].tolist())

        for xv in actual_x:
            assert any(abs(xv - lam) < 0.1 for lam in expected_lams), (
                f"Scatter x={xv:.3f} does not match any wavelength {expected_lams}"
            )
        pdp.close()

    def test_scatter_info_rd_yax_always_stored(self, mol, line_data):
        """rd_yax must always be stored in active_lines info dict even after axis change."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)
        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)

        for entry in active_lines:
            info = entry[3]
            assert "rd_yax" in info, "rd_yax missing from active_lines info after axis change"
            assert isinstance(info["rd_yax"], float)
        pdp.close()

    def test_cache_molecule_reference_preserved(self, mol, line_data):
        """Cache must retain the original molecule reference through re-renders."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        pdp.color_by("e_up", regenerate=True)
        pdp.set_axes(x_prop="eu", y_prop="rd_yax", regenerate=True)

        assert pdp._active_lines_cache["molecule"] is mol
        pdp.close()

    def test_no_active_lines_no_extra_scatter_on_axis_change(self, mol):
        """With no active lines, set_axes must not add spurious scatter points."""
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        # No render_active_lines called

        n_before = sum(len(c.get_offsets()) for c in pdp.ax.collections)
        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        n_after = sum(len(c.get_offsets()) for c in pdp.ax.collections)

        assert n_after == n_before, (
            "Axis change without active lines changed scatter count"
        )
        pdp.close()

    def test_external_axes_scatter_on_correct_ax(self, mol, line_data):
        """With external axes, scatter must appear on that axes object."""
        fig, ax = plt.subplots()
        pdp = PopulationDiagramPlot(molecule=mol, ax=ax)
        pdp.generate_plot()
        active_lines = []
        pdp.render_active_lines(line_data, active_lines, molecule=mol)

        # All scatter must live on ax
        for coll in ax.collections:
            assert coll.axes is ax

        pdp.set_axes(x_prop="wavelength", y_prop="intens", regenerate=True)
        offsets = _scatter_offsets(pdp)
        assert any(len(o) > 0 for o in offsets), "Scatter missing after set_axes on external ax"
        plt.close(fig)
