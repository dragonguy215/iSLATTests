import numpy as np
#import pandas as pd

from typing import Dict, List, Optional, Tuple, Callable, Any, Union, TYPE_CHECKING
from matplotlib.ticker import MaxNLocator

from .BasePlot import BasePlot
from .CompositePlot import CompositePlot
from .SpectrumPanel import SpectrumPanel

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from lmfit.model import ModelResult

class FitLinesPlotGrid(CompositePlot):
    def __init__(self, #files : List[str], 
                 #molecules_dict: 'MoleculeDict', 
                 fit_data: Dict[str, Any] = None,
                 rows: int = None, cols: int = 10, 
                 #figsize: Tuple[int, int] = (15, 15),
                 **kwargs
                 ):
        #self.files = files
        #self.molecules_dict = molecules_dict
        self.fit_csv_dict: Dict[str, Any]
        self.fit_data_tuple_list: List[Tuple[ModelResult, np.ndarray, np.ndarray]]
        self.fit_csv_dict, self.fit_data_tuple_list = fit_data
        #self.gauss_fit, self.fitted_wave, self.fitted_flux = fit_data
        
        # Set rows and cols by the length of fit data if not specified
        if rows is None or cols is None:
            n_plots = len(self.fit_csv_dict)
            if rows is None and cols is None:
                cols = int(np.ceil(np.sqrt(n_plots)))
                rows = int(np.ceil(n_plots / cols))
            elif rows is None:
                rows = int(np.ceil(n_plots / cols))
            elif cols is None:
                cols = int(np.ceil(n_plots / rows))
        self.rows = rows
        self.cols = cols

        self.spectrum_name = kwargs.get('spectrum_name', 'Spectrum')
        self.figsize = kwargs.get('figsize', (2.5 * self.cols, 2 * self.rows))

        # Initialize CompositePlot with computed figsize
        super().__init__(figsize=self.figsize, **kwargs)

        self.plt_extra_range = kwargs.get('plt_extra_range', 0.015)  # extra range to plot for each line
        self.wave_data = kwargs.get('wave_data', None)
        self.flux_data = kwargs.get('flux_data', None)
        self.err_data = kwargs.get('err_data', None)
        self.fit_line_uncertainty = kwargs.get('fit_line_uncertainty', 3.0)
        # axs is set by _build_layout when generate_plot is called
        self.axs: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # CompositePlot abstract-method implementations
    # ------------------------------------------------------------------

    def _build_layout(self) -> np.ndarray:
        """Create the 2-D axes grid for the fit-line subplots.

        Returns
        -------
        np.ndarray
            2-D array of shape ``(rows, cols)`` containing the
            :class:`~matplotlib.axes.Axes` objects.
        """
        axs = self.fig.subplots(self.rows, self.cols)
        # Normalise to a 2-D array regardless of grid shape.
        if self.rows == 1 and self.cols == 1:
            axs = np.array([[axs]])
        elif self.rows == 1:
            axs = axs.reshape(1, -1)
        elif self.cols == 1:
            axs = axs.reshape(-1, 1)
        self.axs = axs
        return axs

    def _create_panels(self, layout: np.ndarray) -> None:
        """Create one :class:`SpectrumPanel` per fit entry and register it.

        Each panel covers the cropped wavelength range
        ``[xmin - extra, xmax + extra]`` for its fit.  The Gaussian
        fit overlay and axis decoration are handled in
        :meth:`_post_render` after all panels have been rendered.
        """
        gauss_fits, fitted_waves, fitted_fluxes = self.fit_data_tuple_list
        for idx, _ in enumerate(zip(gauss_fits, fitted_waves, fitted_fluxes)):
            if idx >= self.rows * self.cols:
                break
            row = idx // self.cols
            col = idx % self.cols
            ax = layout[row, col]

            xmin = self.fit_csv_dict[idx]['xmin']
            xmax = self.fit_csv_dict[idx]['xmax']

            if self.wave_data is None or self.flux_data is None:
                continue

            fit_mask = (
                (self.wave_data >= xmin - self.plt_extra_range)
                & (self.wave_data <= xmax + self.plt_extra_range)
            )
            err_slice = (
                self.err_data[fit_mask]
                if self.err_data is not None
                else None
            )
            panel = SpectrumPanel(
                wave_data=self.wave_data[fit_mask],
                flux_data=self.flux_data[fit_mask],
                xmin=xmin - self.plt_extra_range,
                xmax=xmax + self.plt_extra_range,
                error_data=err_slice,
                ax=ax,
                theme=self.theme,
            )
            self._register_panel(f"fit_{idx}", panel)

    def _post_render(self, **kwargs) -> None:
        """Overlay Gaussian fits; apply titles, tick settings, and hide unused subplots."""
        gauss_fits, fitted_waves, fitted_fluxes = self.fit_data_tuple_list
        fit_pairs = list(zip(gauss_fits, fitted_waves, fitted_fluxes))
        n_plots = len(fit_pairs)

        for idx, (gauss_fit, fitted_wave, fitted_flux) in enumerate(fit_pairs):
            if idx >= self.rows * self.cols:
                break
            panel = self.get_panel(f"fit_{idx}")
            if panel is None:
                continue
            ax = panel.ax  # set by SpectrumPanel._resolve_axes() during generate_plot

            xmin = self.fit_csv_dict[idx]['xmin']
            xmax = self.fit_csv_dict[idx]['xmax']

            if fitted_wave is None or fitted_flux is None:
                ax.set_title(f"Line {idx+1}: Fit Error", fontsize=9, pad=2)
                continue

            line_color = 'lime' if self.fit_csv_dict[idx]['Fit_det'] else 'red'
            try:
                SpectrumPanel.plot_gaussian_fit(
                    ax, gauss_fit, fitted_wave, fitted_flux,
                    color=line_color, uncertainty_sigma=self.fit_line_uncertainty,
                )
                ax.set_title(
                    f"{self.fit_csv_dict[idx]['species']} "
                    f"{self.fit_csv_dict[idx]['lam']:.2f}",
                    fontsize=9, pad=2,
                )
                # y-limits relative to the observed flux in the fit window
                panel_flux = panel.flux_data
                if len(panel_flux) > 0:
                    y_min = np.min(panel_flux) - 0.1 * np.abs(np.min(panel_flux))
                    y_max = np.max(panel_flux) + 0.1 * np.abs(np.max(panel_flux))
                    ax.set_ylim(y_min, y_max)
            except Exception as e:
                ax.set_title("Plot Error", fontsize=9, pad=2)
                print(f"Error plotting line {idx+1}: {e}")

            # Compact tick labels to prevent overlap
            ax.tick_params(axis='both', labelsize=7, pad=1)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

            ax.set_xlabel("Wavelength (μm)", fontsize=7, labelpad=1)

            col = idx % self.cols
            if col == 0:
                ax.set_ylabel("Flux (Jy)", fontsize=7, labelpad=1)

        # Hide unused subplots
        for idx in range(n_plots, self.rows * self.cols):
            row = idx // self.cols
            col = idx % self.cols
            self.axs[row, col].set_visible(False)
    
    def plot(self):
        self.generate_plot()
        self.show(block=False)