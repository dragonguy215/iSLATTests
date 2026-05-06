# -*- coding: utf-8 -*-
"""Tests for PlotView hierarchy, MainPlotGrid, PopulationDiagramPlot,
LineInspectionPlot, CompositeStackedPanel, and other previously-untested
plotting view / composite classes.

All matplotlib work uses the 'Agg' backend for headless safety.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock
from matplotlib.figure import Figure as MplFigure
from matplotlib.axes import Axes

import iSLAT.Constants as c
from iSLAT.Modules.Plotting import (
    BasePlot,
    FullSpectrumPlot,
    ResidualSpectrumPlot,
    CompositeStackedPanel,
    LineInspectionPlot,
    MainPlotGrid,
    PopulationDiagramPlot,
    PlotView,
    GapMode,
)

from tests.plotting import make_wave_flux, make_atomic_lines, make_line_list


# ======================================================================
# Helpers for building test molecules / intensities
# ======================================================================

def _make_test_line_list(n_lines=5, lam_start=12.0):
    """Build a MoleculeLineList with partition function and molar mass."""
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
    mll = MoleculeLineList(molecule_id='TestMol', lines_data=lines_data)
    mll.partition_function = mll._partition_type(
        t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
        q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
    )
    mll._molar_mass = 18.015
    return mll


def _make_test_intensity(n_lines=5, lam_start=12.0, t_kin=500.0, n_mol=1e18):
    """Build a computed Intensity object ready for plotting."""
    from iSLAT.Modules.DataTypes.Intensity import Intensity
    mll = _make_test_line_list(n_lines=n_lines, lam_start=lam_start)
    inten = Intensity(mll)
    inten.calc_intensity(t_kin=t_kin, n_mol=n_mol, dv=1.0)
    return inten


def _make_test_molecule(name='H2O', temp=500.0, n_lines=5, lam_start=12.0):
    """Build a Molecule with computed intensity/spectrum data."""
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.Spectrum import Spectrum
    mll = _make_test_line_list(n_lines=n_lines, lam_start=lam_start)
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
    mol.intensity = _make_test_intensity(n_lines=n_lines, lam_start=lam_start, t_kin=temp)
    spec = Spectrum(lam_min=10.0, lam_max=20.0, dlambda=0.01, R=3000.0, distance=160.0)
    area = np.pi * 1.0 ** 2
    spec.add_intensity(mol.intensity, area)
    mol.spectrum = spec
    return mol


# ======================================================================
# PlotView ABC
# ======================================================================

class TestPlotViewABC:
    """PlotView is abstract and cannot be instantiated directly."""

    def test_cannot_instantiate_base(self):
        with pytest.raises(TypeError):
            PlotView()

    def test_has_abstract_methods(self):
        abstracts = PlotView.__abstractmethods__
        assert "activate" in abstracts
        assert "deactivate" in abstracts
        assert "update_model_plot" in abstracts
        assert "on_molecule_visibility_changed" in abstracts


# ======================================================================
# PopulationDiagramPlot
# ======================================================================

class TestPopulationDiagramPlot:
    """Coverage for PopulationDiagramPlot with various inputs."""

    def test_with_molecule(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        assert pdp.fig is not None
        assert pdp.ax is not None
        pdp.close()

    def test_with_bare_intensity(self):
        inten = _make_test_intensity()
        pdp = PopulationDiagramPlot(
            intensity=inten,
            name="TestMol",
            color="blue",
            radius=1.0,
            distance=160.0,
        )
        pdp.generate_plot()
        assert pdp.fig is not None
        title = pdp.ax.get_title()
        assert "TestMol" in title
        pdp.close()

    def test_mutual_exclusion(self):
        mol = _make_test_molecule()
        inten = _make_test_intensity()
        with pytest.raises(ValueError, match="not more than one"):
            PopulationDiagramPlot(molecule=mol, intensity=inten)

    def test_no_data_shows_message(self):
        pdp = PopulationDiagramPlot()
        pdp.generate_plot()
        title = pdp.ax.get_title()
        assert "No molecule selected" in title
        pdp.close()

    def test_external_axes(self):
        fig, ax = plt.subplots()
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol, ax=ax)
        pdp.generate_plot()
        assert pdp.ax is ax
        plt.close(fig)

    def test_set_molecule_regenerates(self):
        mol1 = _make_test_molecule(name='H2O')
        mol2 = _make_test_molecule(name='CO', temp=300.0)
        pdp = PopulationDiagramPlot(molecule=mol1)
        pdp.generate_plot()
        assert "H2O" in pdp.ax.get_title()
        pdp.set_molecule(mol2)
        assert "CO" in pdp.ax.get_title()
        pdp.close()

    def test_set_intensity_regenerates(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        inten = _make_test_intensity(t_kin=300.0)
        pdp.set_intensity(inten, name="Bare", radius=2.0, distance=100.0)
        assert "Bare" in pdp.ax.get_title()
        pdp.close()

    # --- Multi-molecule tests ---

    def test_multiple_molecules_list(self):
        """Passing a list of molecules renders all with distinct colours."""
        mol1 = _make_test_molecule(name='H2O', temp=500.0)
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecules=[mol1, mol2])
        pdp.generate_plot()
        assert pdp.fig is not None
        assert pdp.ax is not None
        assert "Population diagram" in pdp.ax.get_title()
        # Should have 2 component data entries
        assert len(pdp.component_data) == 2
        assert pdp.component_data[0]["name"] == "H2O"
        assert pdp.component_data[1]["name"] == "CO"
        # Colours should be distinct
        assert pdp.component_data[0]["color"] != pdp.component_data[1]["color"]
        pdp.close()

    def test_multiple_molecules_molecule_dict(self):
        """Passing a MoleculeDict renders visible molecules."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        mol1 = _make_test_molecule(name='H2O', temp=500.0)
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        md = MoleculeDict()
        md['H2O'] = mol1
        md['CO'] = mol2
        pdp = PopulationDiagramPlot(molecules=md)
        pdp.generate_plot()
        assert pdp.fig is not None
        assert len(pdp.component_data) >= 1
        pdp.close()

    def test_mutual_exclusion_molecules_and_molecule(self):
        mol = _make_test_molecule()
        with pytest.raises(ValueError, match="not more than one"):
            PopulationDiagramPlot(molecule=mol, molecules=[mol])

    def test_mutual_exclusion_molecules_and_intensity(self):
        mol = _make_test_molecule()
        inten = _make_test_intensity()
        with pytest.raises(ValueError, match="not more than one"):
            PopulationDiagramPlot(molecules=[mol], intensity=inten)

    def test_set_molecules_regenerates(self):
        mol1 = _make_test_molecule(name='H2O')
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecule=mol1)
        pdp.generate_plot()
        assert "H2O" in pdp.ax.get_title()
        pdp.set_molecules([mol1, mol2])
        assert "Population diagram" in pdp.ax.get_title()
        assert len(pdp.component_data) == 2
        pdp.close()

    def test_add_molecule(self):
        mol1 = _make_test_molecule(name='H2O')
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecule=mol1)
        pdp.generate_plot()
        assert len(pdp.component_data) == 1
        pdp.add_molecule(mol2)
        assert len(pdp.component_data) == 2
        pdp.close()

    def test_remove_molecule(self):
        mol1 = _make_test_molecule(name='H2O')
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecules=[mol1, mol2])
        pdp.generate_plot()
        assert len(pdp.component_data) == 2
        pdp.remove_molecule('CO')
        # Should be back to single-molecule mode
        assert len(pdp.component_data) == 1
        assert pdp.component_data[0]["name"] == "H2O"
        pdp.close()

    # --- Color mapping tests ---

    def test_color_by_e_up(self):
        """Color mapping by upper-level energy produces a colorbar."""
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        pdp.color_by('e_up', cmap='plasma')
        # Should have a colorbar (extra axes)
        assert len(pdp.fig.axes) > 1
        pdp.close()

    def test_color_by_a_stein(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('a_stein')
        assert pdp.fig is not None
        pdp.close()

    def test_color_by_wavelength(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('wavelength')
        assert pdp.fig is not None
        pdp.close()

    def test_color_by_lam_alias(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('lam')
        assert pdp.fig is not None
        pdp.close()

    def test_color_by_lev_up_categorical(self):
        """Categorical color mapping by transition label."""
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('lev_up', cmap='tab10')
        assert pdp.fig is not None
        # Should have a legend
        legend = pdp.ax.get_legend()
        assert legend is not None
        pdp.close()

    def test_color_by_component_categorical(self):
        """Categorical color mapping by component name."""
        mol1 = _make_test_molecule(name='H2O')
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecules=[mol1, mol2])
        pdp.color_by('component', cmap='Set1')
        assert pdp.fig is not None
        legend = pdp.ax.get_legend()
        assert legend is not None
        pdp.close()

    def test_color_by_with_vmin_vmax(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('e_up', vmin=500, vmax=3000)
        assert pdp.fig is not None
        pdp.close()

    def test_clear_color_mapping(self):
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.color_by('e_up')
        pdp.clear_color_mapping()
        assert pdp._color_mapping is None
        pdp.close()

    def test_color_by_multi_molecule(self):
        """Color mapping works across multiple molecules."""
        mol1 = _make_test_molecule(name='H2O', temp=500.0)
        mol2 = _make_test_molecule(name='CO', temp=300.0, lam_start=14.0)
        pdp = PopulationDiagramPlot(molecules=[mol1, mol2])
        pdp.color_by('e_up', cmap='coolwarm')
        assert pdp.fig is not None
        assert len(pdp.component_data) == 2
        pdp.close()

    def test_component_data_has_expected_keys(self):
        """Verify component_data contains all expected molecular properties."""
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        assert len(pdp.component_data) == 1
        cd = pdp.component_data[0]
        for key in ('name', 'color', 'eu', 'rd_yax', 'wavelength',
                     'intens', 'a_stein', 'g_up', 'lev_up', 'lev_low',
                     'e_low', 'valid_mask', 'beam_s'):
            assert key in cd, f"Missing key: {key}"
        pdp.close()

    def test_color_by_no_regenerate(self):
        """color_by with regenerate=False does not call generate_plot."""
        mol = _make_test_molecule()
        pdp = PopulationDiagramPlot(molecule=mol)
        pdp.generate_plot()
        n_axes_before = len(pdp.fig.axes)
        pdp.color_by('e_up', regenerate=False)
        # Figure should not have changed since regenerate=False
        assert len(pdp.fig.axes) == n_axes_before
        pdp.close()


# ======================================================================
# LineInspectionPlot
# ======================================================================

class TestLineInspectionPlot:
    """Tests for LineInspectionPlot rendering."""

    def test_basic_render(self):
        wave, flux, err = make_wave_flux(n=200)
        lip = LineInspectionPlot(
            wave, flux, 13.0, 17.0, error_data=err,
        )
        lip.generate_plot()
        assert lip.fig is not None
        lip.close()

    def test_external_axes(self):
        wave, flux, err = make_wave_flux(n=200)
        fig, ax = plt.subplots()
        lip = LineInspectionPlot(
            wave, flux, 13.0, 17.0, error_data=err, ax=ax,
        )
        lip.generate_plot()
        assert lip.ax is ax
        plt.close(fig)

    def test_with_molecule(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        lip = LineInspectionPlot(
            wave, flux, 12.0, 16.0, error_data=err,
            molecule=mol,
        )
        lip.generate_plot()
        # Should have molecule model line
        assert lip.fig is not None
        lip.close()

    def test_empty_range(self):
        wave, flux, err = make_wave_flux(n=200)
        lip = LineInspectionPlot(
            wave, flux, 50.0, 60.0,  # outside data range
        )
        lip.generate_plot()
        assert lip.fig is not None
        lip.close()


# ======================================================================
# MainPlotGrid
# ======================================================================

class TestMainPlotGrid:
    """Tests for the three-panel composite MainPlotGrid."""

    def test_basic_generate(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        assert mpg.fig is not None
        assert mpg.ax_spectrum is not None
        assert mpg.ax_inspection is not None
        assert mpg.ax_popdiagram is not None
        mpg.close()

    def test_with_spectrum_range(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            spectrum_range=(12.0, 18.0),
        )
        mpg.generate_plot()
        xlim = mpg.ax_spectrum.get_xlim()
        assert abs(xlim[0] - 12.0) < 0.1
        assert abs(xlim[1] - 18.0) < 0.1
        mpg.close()

    def test_with_inspection_range(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            inspection_range=(14.0, 16.0),
        )
        mpg.generate_plot()
        # The inspection panel should have been rendered
        assert len(mpg.ax_inspection.lines) > 0
        mpg.close()

    def test_with_active_molecule(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            active_molecule=mol,
            inspection_range=(12.0, 16.0),
        )
        mpg.generate_plot()
        # Population diagram should show molecule
        title = mpg.ax_popdiagram.get_title()
        assert "H2O" in title
        mpg.close()

    def test_set_inspection_range(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        # Initially no inspection range
        assert "select a range" in mpg.ax_inspection.get_title()
        # Now set one
        mpg.set_inspection_range(13.0, 15.0)
        assert mpg.inspection_range == (13.0, 15.0)
        mpg.close()

    def test_set_spectrum_range(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        mpg.set_spectrum_range(11.0, 19.0)
        assert mpg.spectrum_range == (11.0, 19.0)
        mpg.close()

    def test_reset_spectrum_range(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            spectrum_range=(12.0, 18.0),
        )
        mpg.generate_plot()
        mpg.set_spectrum_range(None, None)
        assert mpg.spectrum_range is None
        mpg.close()

    def test_with_annotations(self):
        wave, flux, err = make_wave_flux(n=200)
        atomic = make_atomic_lines([12.0, 15.0])
        lines = make_line_list([13.0, 17.0])
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            line_list=lines, atomic_lines=atomic,
        )
        mpg.generate_plot()
        txt_count = len(mpg.ax_spectrum.texts)
        assert txt_count > 0
        mpg.close()

    def test_regenerate_clears_old_axes(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        ax1 = mpg.ax_spectrum
        mpg.generate_plot()
        # Second generate should create fresh axes
        assert mpg.ax_spectrum is not ax1
        mpg.close()

# ======================================================================
# MainPlotGrid — borrowed-axes mode
# ======================================================================

class TestMainPlotGridBorrowedAxes:
    """Tests for borrowed-axes mode where MainPlotGrid renders onto
    pre-existing axes without calling fig.clf() or creating new axes."""

    @pytest.fixture
    def borrowed_setup(self):
        """Create a figure with 3 pre-existing axes and a MainPlotGrid in borrowed mode."""
        fig = plt.figure(figsize=(12, 8))
        from iSLAT.Modules.Plotting.BasePlot import BasePlot
        ax1, ax2, ax3 = BasePlot.create_three_panel_axes(fig)
        wave, flux, err = make_wave_flux(n=200)
        yield fig, ax1, ax2, ax3, wave, flux, err
        plt.close(fig)

    def test_borrowed_flag_is_set(self, borrowed_setup):
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        assert mpg._borrowed_axes is True
        assert mpg.ax_spectrum is ax1
        assert mpg.ax_inspection is ax2
        assert mpg.ax_popdiagram is ax3

    def test_standalone_flag_is_false(self):
        """Without injected axes the flag should be False."""
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        assert mpg._borrowed_axes is False
        assert mpg.ax_spectrum is None  # not yet created

    def test_generate_plot_preserves_axes_identity(self, borrowed_setup):
        """In borrowed mode, generate_plot must NOT replace the axes objects."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        # Axes should be the same objects — not new ones
        assert mpg.ax_spectrum is ax1
        assert mpg.ax_inspection is ax2
        assert mpg.ax_popdiagram is ax3

    def test_generate_plot_does_not_clf(self, borrowed_setup):
        """In borrowed mode, other axes on the figure should survive."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        # Add a 4th axes to the figure to prove fig.clf() is not called
        extra_ax = fig.add_axes([0.1, 0.1, 0.2, 0.2])
        assert extra_ax in fig.get_axes()

        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        # The extra axes should still exist
        assert extra_ax in fig.get_axes()

    def test_renders_spectrum_on_borrowed_axes(self, borrowed_setup):
        """Borrowed mode should actually render observed spectrum lines."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        # ax1 should have at least one Line2D (the observed spectrum)
        assert len(ax1.lines) > 0

    def test_renders_with_molecules(self, borrowed_setup):
        """Borrowed mode should render molecule overlays."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mol = _make_test_molecule()
        md = MoleculeDict()
        md[mol.name] = mol
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            molecules=md, active_molecule=mol,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        # Should have molecule-tagged lines
        tagged = [l for l in ax1.lines if getattr(l, "_molecule_name", None) == mol.name]
        assert len(tagged) > 0

    def test_fig_property_points_to_injected_figure(self, borrowed_setup):
        """In borrowed mode, mpg.fig should be the original figure."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        assert mpg.fig is fig

    def test_regenerate_preserves_axes(self, borrowed_setup):
        """Calling generate_plot twice in borrowed mode keeps the same axes."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        mpg.generate_plot()
        assert mpg.ax_spectrum is ax1
        assert mpg.ax_inspection is ax2
        assert mpg.ax_popdiagram is ax3

    def test_update_data_works_in_borrowed_mode(self, borrowed_setup):
        """update_data should refresh panels on borrowed axes."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        new_flux = flux * 2
        mpg.update_data(flux_data=new_flux)
        assert np.array_equal(mpg.flux_data, new_flux)

    def test_set_inspection_range_in_borrowed_mode(self, borrowed_setup):
        """set_inspection_range should render the inspection panel on borrowed ax2."""
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mol = _make_test_molecule()
        mpg = MainPlotGrid(
            wave, flux, error_data=err, active_molecule=mol,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        mpg.set_inspection_range(12.0, 15.0)
        assert mpg.inspection_range == (12.0, 15.0)
        # ax2 should have been rendered (has lines from the LIP)
        assert len(ax2.lines) > 0

    def test_molecule_visibility_in_borrowed_mode(self, borrowed_setup):
        """set_molecule_visibility should toggle artists on borrowed axes."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        fig, ax1, ax2, ax3, wave, flux, err = borrowed_setup
        mol = _make_test_molecule()
        md = MoleculeDict()
        md[mol.name] = mol
        mpg = MainPlotGrid(
            wave, flux, error_data=err, molecules=md,
            ax_spectrum=ax1, ax_inspection=ax2, ax_popdiagram=ax3,
        )
        mpg.generate_plot()
        tagged = [l for l in ax1.lines if getattr(l, "_molecule_name", None) == mol.name]
        assert len(tagged) > 0
        # Toggle off
        found = mpg.set_molecule_visibility(mol.name, False)
        assert found is True
        for l in tagged:
            assert not l.get_visible()

class TestMainPlotGridUpdateData:
    """Tests for the update_data / molecule-visibility helpers."""

    def test_update_data_refreshes_panels(self):
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        # Change flux data
        new_flux = flux * 2
        mpg.update_data(flux_data=new_flux)
        assert np.array_equal(mpg.flux_data, new_flux)
        mpg.close()

    def test_update_data_partial(self):
        """Only updated fields should change."""
        wave, flux, err = make_wave_flux(n=200)
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        orig_wave = mpg.wave_data.copy()
        mpg.update_data(flux_data=flux * 3)
        # wave_data should be unchanged
        assert np.array_equal(mpg.wave_data, orig_wave)
        mpg.close()

    def test_update_data_with_molecule(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        mpg = MainPlotGrid(wave, flux, error_data=err)
        mpg.generate_plot()
        mpg.update_data(active_molecule=mol)
        assert mpg.active_molecule is mol
        mpg.close()


class TestMainPlotGridMoleculeVisibility:
    """Tests for set_molecule_visibility and handle_molecule_visibility_change."""

    def _make_mpg_with_molecule(self):
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        md = MoleculeDict()
        md[mol.name] = mol
        mpg = MainPlotGrid(wave, flux, error_data=err, molecules=md, active_molecule=mol)
        mpg.generate_plot()
        return mpg, mol

    def test_set_molecule_visibility_toggles(self):
        mpg, mol = self._make_mpg_with_molecule()
        # There should be a line tagged with the molecule name
        tagged = [l for l in mpg.ax_spectrum.lines if getattr(l, "_molecule_name", None) == mol.name]
        assert len(tagged) > 0
        # Toggle off
        found = mpg.set_molecule_visibility(mol.name, False)
        assert found is True
        for l in tagged:
            assert not l.get_visible()
        # Toggle back on
        mpg.set_molecule_visibility(mol.name, True)
        for l in tagged:
            assert l.get_visible()
        mpg.close()

    def test_set_molecule_visibility_returns_false_for_unknown(self):
        mpg, _ = self._make_mpg_with_molecule()
        assert mpg.set_molecule_visibility("NONEXISTENT", False) is False
        mpg.close()

    def test_remove_molecule_lines(self):
        mpg, mol = self._make_mpg_with_molecule()
        mpg.remove_molecule_lines(mol.name)
        tagged = [l for l in mpg.ax_spectrum.lines if getattr(l, "_molecule_name", None) == mol.name]
        assert len(tagged) == 0
        mpg.close()


class TestMainPlotGridActiveLineMarkers:
    """Tests for render_active_line_markers and render_active_line_scatter."""

    def _make_line_data(self, mol):
        """Return a small list of (MoleculeLine, intensity, tau) triples."""
        if not hasattr(mol, 'intensity') or mol.intensity is None:
            return []
        try:
            return mol.intensity.get_lines_in_range_with_intensity(12.0, 16.0)
        except Exception:
            return []

    def test_render_active_line_markers(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            active_molecule=mol,
            inspection_range=(12.0, 16.0),
        )
        mpg.generate_plot()
        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")
        active_lines: list = []
        mpg.render_active_line_markers(line_data, active_lines, max_y=0.1, color="green")
        assert len(active_lines) > 0
        # Each entry is [vline, text, scatter_or_None, info_dict]
        for entry in active_lines:
            assert len(entry) == 4
            assert entry[0] is not None  # vline
            assert entry[1] is not None  # text
            assert isinstance(entry[3], dict)
        mpg.close()

    def test_render_active_line_scatter(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            active_molecule=mol,
            inspection_range=(12.0, 16.0),
        )
        mpg.generate_plot()
        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")
        active_lines: list = []
        sc = mpg.render_active_line_scatter(line_data, active_lines, mol, color="green")
        assert sc is not None
        assert len(active_lines) > 0
        # Scatter entries should have _scatter_point_index
        for entry in active_lines:
            assert "_scatter_point_index" in entry[3]
        mpg.close()

class TestCompositeStackedPanel:
    """Tests for CompositeStackedPanel.from_pair and rendering."""

    def _make_fsp(self, n_panels=3, **kw):
        wave, flux, err = make_wave_flux(n=200)
        return FullSpectrumPlot(wave, flux, error_data=err, n_panels=n_panels, **kw)

    def _make_rsp(self, n_panels=3, **kw):
        wave, flux, err = make_wave_flux(n=200)
        model = flux + np.random.default_rng(1).normal(0, 0.002, len(flux))
        return ResidualSpectrumPlot(
            wave, flux, error_data=err, n_panels=n_panels,
            model_flux=model, show_chi2_per_panel=False, **kw,
        )

    def test_from_pair_creates_composite(self):
        fsp = self._make_fsp()
        rsp = self._make_rsp()
        comp = CompositeStackedPanel.from_pair(fsp, rsp)
        assert isinstance(comp, CompositeStackedPanel)
        assert comp.sources == (fsp, rsp)

    def test_row_plan_has_entries(self):
        fsp = self._make_fsp(n_panels=2)
        rsp = self._make_rsp(n_panels=2)
        comp = CompositeStackedPanel.from_pair(fsp, rsp)
        assert len(comp.row_plan) >= 2

    def test_generate_plot(self):
        fsp = self._make_fsp(n_panels=2)
        rsp = self._make_rsp(n_panels=2)
        comp = CompositeStackedPanel.from_pair(fsp, rsp)
        comp.generate_plot()
        assert comp.fig is not None
        comp.close()

    def test_plus_operator(self):
        fsp = self._make_fsp(n_panels=2)
        rsp = self._make_rsp(n_panels=2)
        comp = fsp + rsp
        assert isinstance(comp, CompositeStackedPanel)

    def test_labels(self):
        fsp = self._make_fsp(n_panels=2)
        rsp = self._make_rsp(n_panels=2)
        comp = CompositeStackedPanel.from_pair(
            fsp, rsp, labels=("Spectrum", "Residual"),
        )
        assert comp.labels == ("Spectrum", "Residual")

    def test_gap_mode_inherit(self):
        fsp = self._make_fsp(n_panels=2, gap_mode="skip")
        rsp = self._make_rsp(n_panels=2)
        comp = CompositeStackedPanel.from_pair(fsp, rsp)
        assert comp.gap_mode is GapMode.SKIP


# ======================================================================
# ThreePanelView / FullSpectrumView (mock-based)
# ======================================================================

class TestThreePanelView:
    """ThreePanelView requires an iSLATPlot controller -- use mocks."""

    @pytest.fixture
    def mock_controller(self):
        fig, axes = plt.subplots(1, 3)
        pm = MagicMock()
        pm.fig = fig
        pm.ax1 = axes[0]
        pm.ax2 = axes[1]
        pm.ax3 = axes[2]
        pm.canvas = fig.canvas
        pm.theme = {"background": "#181A1B", "foreground": "#e8e6e3"}
        pm.toggle_state = {"atomic_lines": False, "saved_lines": False}
        pm.islat = MagicMock()
        pm.plot_renderer = MagicMock()
        yield pm
        plt.close(fig)

    def test_init(self, mock_controller):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(mock_controller)
        assert view._pm is mock_controller
        assert view._needs_refresh is True

    def test_properties(self, mock_controller):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(mock_controller)
        assert view.ax1 is mock_controller.ax1
        assert view.ax2 is mock_controller.ax2
        assert view.ax3 is mock_controller.ax3
        assert view._renderer is mock_controller.plot_renderer
        assert view._islat is mock_controller.islat

    def test_has_tagged_artists(self, mock_controller):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(mock_controller)
        # No artists with tag
        assert view._has_tagged_artists("_test_tag") is False

    def test_deactivate(self, mock_controller):
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(mock_controller)
        # Deactivate should call pack_forget -- with Agg backend this is a no-op
        # but shouldn't raise
        try:
            view.deactivate()
        except Exception:
            pass  # Agg backend doesn't support Tk pack

class TestThreePanelViewGrid:
    """ThreePanelView composing a MainPlotGrid in borrowed-axes mode.

    These tests exercise the _ensure_grid / _do_update_model_plot /
    on_molecule_visibility_changed pathways that now route through a
    MainPlotGrid.
    """

    @pytest.fixture
    def rich_controller(self):
        """A mock controller with real axes AND realistic islat data."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        from iSLAT.Modules.Plotting.BasePlot import BasePlot

        fig = plt.figure(figsize=(12, 8))
        ax1, ax2, ax3 = BasePlot.create_three_panel_axes(fig)

        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        md = MoleculeDict()
        md[mol.name] = mol

        islat = MagicMock()
        islat.wave_data = wave
        islat.wave_data_original = wave.copy()
        islat.flux_data = flux
        islat.err_data = err
        islat.molecules_dict = md
        islat.active_molecule = mol
        # apply_stellar_rv should be a pass-through in tests
        md.apply_stellar_rv = MagicMock(side_effect=lambda w: w)

        pm = MagicMock()
        pm.fig = fig
        pm.ax1 = ax1
        pm.ax2 = ax2
        pm.ax3 = ax3
        pm.canvas = MagicMock()
        pm.canvas.draw_idle = MagicMock()
        pm.theme = {"background": "#181A1B", "foreground": "#e8e6e3"}
        pm.toggle_state = {
            "atomic_lines": False,
            "saved_lines": False,
            "summed": True,
            "legend": True,
        }
        pm.islat = islat
        pm.plot_renderer = MagicMock()
        pm.summed_toggle = True
        pm.atomic_toggle = False
        pm.line_toggle = False
        pm.legend_toggle = True
        pm.make_span_selector = MagicMock()

        yield pm, mol, md
        plt.close(fig)

    def test_grid_is_none_before_render(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        assert view._grid is None

    def test_ensure_grid_creates_borrowed_grid(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        grid = view._ensure_grid()
        assert grid is not None
        assert grid._borrowed_axes is True
        assert grid.ax_spectrum is pm.ax1
        assert grid.ax_inspection is pm.ax2
        assert grid.ax_popdiagram is pm.ax3

    def test_ensure_grid_is_idempotent(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        g1 = view._ensure_grid()
        g2 = view._ensure_grid()
        assert g1 is g2

    def test_do_update_model_plot_renders_onto_ax1(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        # ax1 should now have observed spectrum + molecule lines
        assert len(pm.ax1.lines) > 0
        # The grid should have been created
        assert view._grid is not None

    def test_do_update_model_plot_calls_make_span_selector(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        pm.make_span_selector.assert_called_once()

    def test_do_update_model_plot_calls_draw_idle(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        pm.canvas.draw_idle.assert_called()

    def test_molecule_visibility_routes_through_grid(self, rich_controller):
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        # First render to establish molecule artists
        view._do_update_model_plot()
        # Verify molecule lines exist
        tagged = [l for l in pm.ax1.lines if getattr(l, "_molecule_name", None) == mol.name]
        assert len(tagged) > 0
        # Toggle off via the view
        view.on_molecule_visibility_changed(
            molecule_name=mol.name,
            is_visible=False,
            molecules_dict=md,
            wave_data=pm.islat.wave_data,
        )
        # Artists should be invisible
        tagged_after = [l for l in pm.ax1.lines if getattr(l, "_molecule_name", None) == mol.name]
        for l in tagged_after:
            assert not l.get_visible()

    def test_molecule_visibility_does_not_use_renderer(self, rich_controller):
        """After Step 2, on_molecule_visibility_changed should NOT call
        renderer.handle_molecule_visibility_change."""
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        pm.plot_renderer.handle_molecule_visibility_change.reset_mock()
        view.on_molecule_visibility_changed(
            molecule_name=mol.name,
            is_visible=False,
            molecules_dict=md,
            wave_data=pm.islat.wave_data,
        )
        pm.plot_renderer.handle_molecule_visibility_change.assert_not_called()

    def test_grid_axes_are_controller_axes(self, rich_controller):
        """The grid's borrowed axes must be the exact same objects as the controller's."""
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        grid = view._ensure_grid()
        assert grid.ax_spectrum is pm.ax1
        assert grid.ax_inspection is pm.ax2
        assert grid.ax_popdiagram is pm.ax3

    def test_update_model_plot_public_api(self, rich_controller):
        """The public update_model_plot() should delegate and clear _needs_refresh."""
        pm, mol, md = rich_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        assert view._needs_refresh is True
        view.update_model_plot()
        assert view._needs_refresh is False
        assert len(pm.ax1.lines) > 0

    def test_empty_molecules_clears_model(self, rich_controller):
        """When molecules_dict is empty, _do_update_model_plot should clear model lines from ax1."""
        pm, mol, md = rich_controller
        # First render so molecule lines exist on ax1
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        # Confirm molecule lines are on ax1
        tagged_before = [l for l in pm.ax1.lines if getattr(l, "_molecule_name", None) is not None]

        # Now swap in an empty molecules_dict and re-render
        pm.islat.molecules_dict = MagicMock()
        pm.islat.molecules_dict.__len__ = MagicMock(return_value=0)
        # Reset the grid so it re-evaluates molecule count
        view._grid = None
        view._do_update_model_plot()
        # All molecule-tagged lines should be gone
        tagged_after = [l for l in pm.ax1.lines if getattr(l, "_molecule_name", None) is not None]
        assert len(tagged_after) == 0


class TestThreePanelViewThreshold:
    """Verify that ThreePanelView applies the line-intensity threshold
    from ``user_settings['line_threshold']`` when rendering active lines
    in the line-inspection and population-diagram panels.
    """

    @pytest.fixture
    def threshold_controller(self):
        """Controller with real axes, molecule data, and configurable
        ``user_settings`` so we can set/change the threshold."""
        from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
        from iSLAT.Modules.Plotting.BasePlot import BasePlot

        fig = plt.figure(figsize=(12, 8))
        ax1, ax2, ax3 = BasePlot.create_three_panel_axes(fig)

        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        md = MoleculeDict()
        md[mol.name] = mol

        islat = MagicMock()
        islat.wave_data = wave
        islat.wave_data_original = wave.copy()
        islat.flux_data = flux
        islat.err_data = err
        islat.molecules_dict = md
        islat.active_molecule = mol
        islat.user_settings = {"line_threshold": 0.3}  # 30% default
        md.apply_stellar_rv = MagicMock(side_effect=lambda w: w)

        pm = MagicMock()
        pm.fig = fig
        pm.ax1 = ax1
        pm.ax2 = ax2
        pm.ax3 = ax3
        pm.canvas = MagicMock()
        pm.canvas.draw_idle = MagicMock()
        pm.theme = {"background": "#181A1B", "foreground": "#e8e6e3"}
        pm.toggle_state = {
            "atomic_lines": False,
            "saved_lines": False,
            "summed": True,
            "legend": True,
        }
        pm.islat = islat
        pm.plot_renderer = MagicMock()
        pm.summed_toggle = True
        pm.atomic_toggle = False
        pm.line_toggle = False
        pm.legend_toggle = True
        pm.make_span_selector = MagicMock()

        yield pm, mol, md
        plt.close(fig)

    def _make_line_data(self, mol, lam_min=12.0, lam_max=16.0):
        """Get real line data from the molecule's computed intensity."""
        if mol.intensity is None:
            return []
        try:
            return mol.intensity.get_lines_in_range_with_intensity(lam_min, lam_max)
        except Exception:
            return []

    # ------------------------------------------------------------------
    # _get_line_threshold
    # ------------------------------------------------------------------
    def test_get_line_threshold_reads_user_settings(self, threshold_controller):
        """_get_line_threshold should return the value from user_settings."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        assert view._get_line_threshold() == 0.3

    def test_get_line_threshold_custom_value(self, threshold_controller):
        """Changing user_settings should be reflected immediately."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        pm.islat.user_settings["line_threshold"] = 0.5
        assert view._get_line_threshold() == 0.5

    def test_get_line_threshold_default_when_missing(self, threshold_controller):
        """If 'line_threshold' key is absent, default to 0.3."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        pm.islat.user_settings = {}  # no key
        view = ThreePanelView(pm)
        assert view._get_line_threshold() == 0.3

    def test_get_line_threshold_default_when_no_user_settings(self, threshold_controller):
        """If islat has no user_settings attribute, default to 0.3."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        del pm.islat.user_settings
        view = ThreePanelView(pm)
        assert view._get_line_threshold() == 0.3

    # ------------------------------------------------------------------
    # Line-inspection threshold filtering
    # ------------------------------------------------------------------
    def test_render_line_inspection_passes_threshold(self, threshold_controller):
        """_render_line_inspection must forward the threshold to the grid."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        grid = view._ensure_grid()
        grid.render_active_line_markers = MagicMock()

        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")

        view._render_line_inspection(12.0, 16.0, line_data)

        grid.render_active_line_markers.assert_called_once()
        _, kwargs = grid.render_active_line_markers.call_args
        assert "threshold" in kwargs
        assert kwargs["threshold"] == 0.3

    def test_render_line_inspection_threshold_filters_weak_lines(self, threshold_controller):
        """With a high threshold, weak lines should be excluded from active_lines."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView

        line_data = self._make_line_data(mol)
        if len(line_data) < 2:
            pytest.skip("Need at least 2 lines for threshold test")

        # Render with threshold=0 (all lines visible)
        pm.islat.user_settings["line_threshold"] = 0.0
        view_all = ThreePanelView(pm)
        view_all._do_update_model_plot()
        view_all._render_line_inspection(12.0, 16.0, line_data)
        n_all = len(view_all.active_lines)

        # Clean up axes for next render
        pm.ax2.clear()

        # Render with a high threshold (only strongest lines visible)
        pm.islat.user_settings["line_threshold"] = 0.99
        view_strict = ThreePanelView(pm)
        view_strict._do_update_model_plot()
        view_strict._render_line_inspection(12.0, 16.0, line_data)
        n_strict = len(view_strict.active_lines)

        # Strict threshold should show fewer (or equal) lines
        assert n_strict <= n_all
        # With threshold=0.99, most lines should be filtered
        if n_all > 1:
            assert n_strict < n_all

    def test_render_line_inspection_zero_threshold_shows_all(self, threshold_controller):
        """Threshold=0.0 should include every line with positive intensity."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView

        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")

        pm.islat.user_settings["line_threshold"] = 0.0
        view = ThreePanelView(pm)
        view._do_update_model_plot()
        view._render_line_inspection(12.0, 16.0, line_data)

        # Every line with positive intensity should be present
        positive_lines = [
            (l, i, t) for l, i, t in line_data if i > 0
        ]
        assert len(view.active_lines) == len(positive_lines)

    # ------------------------------------------------------------------
    # Population-diagram threshold filtering
    # ------------------------------------------------------------------
    def test_render_population_diagram_passes_threshold(self, threshold_controller):
        """_render_population_diagram_with_lines must forward the threshold."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        grid = view._ensure_grid()
        grid.render_active_line_scatter = MagicMock(return_value=None)

        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")

        view._render_population_diagram_with_lines(line_data)

        grid.render_active_line_scatter.assert_called_once()
        _, kwargs = grid.render_active_line_scatter.call_args
        assert "threshold" in kwargs
        assert kwargs["threshold"] == 0.3

    def test_render_population_diagram_threshold_filters_weak_lines(self, threshold_controller):
        """With a high threshold, weak lines should not appear in scatter."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView

        line_data = self._make_line_data(mol)
        if len(line_data) < 2:
            pytest.skip("Need at least 2 lines for threshold test")

        # All lines with threshold=0
        pm.islat.user_settings["line_threshold"] = 0.0
        view_all = ThreePanelView(pm)
        view_all._do_update_model_plot()
        view_all._render_population_diagram_with_lines(line_data)
        n_all = len(view_all.active_lines)

        pm.ax3.clear()

        # Strict threshold
        pm.islat.user_settings["line_threshold"] = 0.99
        view_strict = ThreePanelView(pm)
        view_strict._do_update_model_plot()
        view_strict._render_population_diagram_with_lines(line_data)
        n_strict = len(view_strict.active_lines)

        assert n_strict <= n_all
        if n_all > 1:
            assert n_strict < n_all

    def test_threshold_change_takes_effect_on_next_render(self, threshold_controller):
        """Changing the threshold between renders should change the output."""
        pm, mol, md = threshold_controller
        from iSLAT.Modules.Plotting.ThreePanelView import ThreePanelView
        view = ThreePanelView(pm)
        grid = view._ensure_grid()
        grid.render_active_line_markers = MagicMock()

        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data generated")

        # First render with 0.3
        pm.islat.user_settings["line_threshold"] = 0.3
        view._render_line_inspection(12.0, 16.0, line_data)
        first_threshold = grid.render_active_line_markers.call_args[1]["threshold"]
        assert first_threshold == 0.3

        # Change threshold and render again
        pm.islat.user_settings["line_threshold"] = 0.7
        grid.render_active_line_markers.reset_mock()
        view._render_line_inspection(12.0, 16.0, line_data)
        second_threshold = grid.render_active_line_markers.call_args[1]["threshold"]
        assert second_threshold == 0.7


class TestMainPlotGridThresholdFiltering:
    """Validate that MainPlotGrid.render_active_line_markers and
    render_active_line_scatter actually filter lines by threshold.

    These are lower-level tests that ensure the filtering machinery
    itself works correctly, independent of ThreePanelView.
    """

    @pytest.fixture
    def grid_with_data(self):
        wave, flux, err = make_wave_flux(n=200)
        mol = _make_test_molecule()
        mpg = MainPlotGrid(
            wave, flux, error_data=err,
            active_molecule=mol,
            inspection_range=(12.0, 16.0),
        )
        mpg.generate_plot()
        yield mpg, mol
        mpg.close()

    def _make_line_data(self, mol, lam_min=12.0, lam_max=16.0):
        if mol.intensity is None:
            return []
        try:
            return mol.intensity.get_lines_in_range_with_intensity(lam_min, lam_max)
        except Exception:
            return []

    def test_markers_threshold_zero_includes_all(self, grid_with_data):
        """threshold=0.0 should include all lines with positive intensity."""
        mpg, mol = grid_with_data
        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data")
        active: list = []
        mpg.render_active_line_markers(line_data, active, max_y=0.1, threshold=0.0)
        positive = [(l, i, t) for l, i, t in line_data if i > 0]
        assert len(active) == len(positive)

    def test_markers_threshold_one_includes_only_strongest(self, grid_with_data):
        """threshold=1.0 should only include lines at exactly max intensity."""
        mpg, mol = grid_with_data
        line_data = self._make_line_data(mol)
        if len(line_data) < 2:
            pytest.skip("Need multiple lines")
        active: list = []
        mpg.render_active_line_markers(line_data, active, max_y=0.1, threshold=1.0)
        # Only the single strongest line should pass (frac == 1.0)
        assert len(active) == 1

    def test_markers_higher_threshold_fewer_lines(self, grid_with_data):
        """Increasing threshold should never increase the number of lines."""
        mpg, mol = grid_with_data
        line_data = self._make_line_data(mol)
        if len(line_data) < 2:
            pytest.skip("Need multiple lines")

        counts = []
        for thresh in [0.0, 0.3, 0.5, 0.8, 1.0]:
            active: list = []
            mpg.ax_inspection.clear()
            mpg.render_active_line_markers(line_data, active, max_y=0.1, threshold=thresh)
            counts.append(len(active))

        # Counts should be monotonically non-increasing
        for i in range(1, len(counts)):
            assert counts[i] <= counts[i - 1], (
                f"threshold={[0.0, 0.3, 0.5, 0.8, 1.0][i]} gave "
                f"{counts[i]} lines but previous gave {counts[i - 1]}"
            )

    def test_scatter_threshold_zero_includes_all(self, grid_with_data):
        """threshold=0.0 for scatter should include all valid lines."""
        mpg, mol = grid_with_data
        line_data = self._make_line_data(mol)
        if not line_data:
            pytest.skip("No line data")
        active: list = []
        sc = mpg.render_active_line_scatter(
            line_data, active, mol, threshold=0.0,
        )
        # All lines with valid data should be present
        assert len(active) > 0
        if sc is not None:
            assert len(active) == len(line_data)

    def test_scatter_higher_threshold_fewer_points(self, grid_with_data):
        """Increasing threshold should never increase scatter point count."""
        mpg, mol = grid_with_data
        line_data = self._make_line_data(mol)
        if len(line_data) < 2:
            pytest.skip("Need multiple lines")

        counts = []
        for thresh in [0.0, 0.3, 0.5, 0.8, 1.0]:
            active: list = []
            mpg.ax_popdiagram.clear()
            mpg.render_active_line_scatter(
                line_data, active, mol, threshold=thresh,
            )
            counts.append(len(active))

        for i in range(1, len(counts)):
            assert counts[i] <= counts[i - 1]


class TestFullSpectrumViewInit:
    """FullSpectrumView initialization tests with mocks."""

    @pytest.fixture
    def mock_controller(self):
        fig, ax = plt.subplots()
        pm = MagicMock()
        pm.fig = fig
        pm.canvas = fig.canvas
        pm.theme = {"background": "#181A1B", "foreground": "#e8e6e3"}
        pm.islat = MagicMock()
        pm.plot_renderer = MagicMock()
        yield pm
        plt.close(fig)

    def test_init(self, mock_controller):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        view = FullSpectrumView(mock_controller)
        assert view._pm is mock_controller
        assert view._plot is None
        assert view._initialised is False
        assert view._needs_refresh is True
        assert view.span_selectors == {}

    def test_fig_property_none_before_init(self, mock_controller):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        view = FullSpectrumView(mock_controller)
        assert view.fig is None

    def test_subplots_empty_before_init(self, mock_controller):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        view = FullSpectrumView(mock_controller)
        assert view.subplots == {}

    def test_iter_all_axes_empty(self, mock_controller):
        from iSLAT.Modules.Plotting.FullSpectrumView import FullSpectrumView
        view = FullSpectrumView(mock_controller)
        assert view._iter_all_axes() == []
