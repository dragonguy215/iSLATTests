"""
iSLAT Plotting Module
=====================

Provides standalone, class-based plot objects that work **both** inside the
iSLAT GUI and as regular matplotlib figures in scripts / Jupyter notebooks.

Quick-start (notebook)::

    from iSLAT.Modules.Plotting import (
        LineInspectionPlot,
        PopulationDiagramPlot,
        FullSpectrumPlot,
        MainPlotGrid,
        FitLinesPlotGrid,
    )

Class hierarchy::

    BasePlot (ABC)
    ├── SpectralPanel (ABC)        — single panel of spectral data in a range
    │   ├── LineInspectionPlot     — zoomed wavelength region
    │   ├── SpectrumPanel          — concrete panel: spectrum + molecules
    │   └── ResidualPanel          — concrete panel: residual + chi^2
    ├── StackedSpectralPanel (ABC) — composer for vertically stacked panels
    │   ├── FullSpectrumPlot       — multi-panel full spectrum overview
    │   │   └── ResidualSpectrumPlot — ...with per-panel residuals & chi^2
    │   └── CompositeStackedPanel  — merged cells from 2 source plots
    ├── PopulationDiagramPlot      — Boltzmann / rotation diagram
    ├── MainPlotGrid               — 3-panel composite (spectrum + inspection + pop-diagram)
    └── FitLinesPlotGrid           — grid of individual line-fit results

    GapMode (Enum)             — CONNECT (default) | SKIP (break lines at gaps)
    XScaling (Enum)            — WAVELENGTH (default) | DATA_DENSITY (uniform point density)

    PlotView (ABC)             — switchable view interface (GUI only)
    ├── ThreePanelView         — standard 3-panel GUI layout
    └── FullSpectrumView       — multi-panel full spectrum GUI layout
"""

from .BasePlot import BasePlot, DEFAULT_THEME, _detect_system_theme
from .BasePlot import BasePlot as _BP
load_theme = _BP.load_theme  # Convenience alias at package level
from .SpectralPanel import SpectralPanel, GapMode, XScaling
from .StackedSpectralPanel import StackedSpectralPanel
from .CompositeStackedPanel import CompositeStackedPanel
from .SpectrumPanel import SpectrumPanel
from .ResidualPanel import ResidualPanel
from .LineInspectionPlot import LineInspectionPlot
from .PopulationDiagramPlot import PopulationDiagramPlot
from .FullSpectrumPlot import FullSpectrumPlot
from .ResidualSpectrumPlot import ResidualSpectrumPlot
from .MainPlotGrid import MainPlotGrid
from .FitLinesPlotGrid import FitLinesPlotGrid
from .PlotView import PlotView
from .ThreePanelView import ThreePanelView
from .FullSpectrumView import FullSpectrumView

__all__ = [
    "BasePlot",
    "DEFAULT_THEME",
    "load_theme",
    "GapMode",
    "XScaling",
    "SpectralPanel",
    "StackedSpectralPanel",
    "CompositeStackedPanel",
    "SpectrumPanel",
    "ResidualPanel",
    "LineInspectionPlot",
    "PopulationDiagramPlot",
    "FullSpectrumPlot",
    "MainPlotGrid",
    "FitLinesPlotGrid",
    "PlotView",
    "ThreePanelView",
    "FullSpectrumView",
    "ResidualSpectrumPlot",
]