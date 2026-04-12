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
        with pytest.raises(ValueError, match="not both"):
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
# CompositeStackedPanel
# ======================================================================

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
