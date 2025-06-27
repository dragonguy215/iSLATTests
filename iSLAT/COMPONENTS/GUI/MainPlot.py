import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import SpanSelector
import numpy as np
from lmfit.models import GaussianModel
import pandas as pd

class iSLATPlot:
    def __init__(self, master, input_spectrum_data, molecules_data, xp1, xp2, wave_data, flux_data, config):
        self.master = master
        self.input_spectrum_data = input_spectrum_data
        self.molecules_data = molecules_data
        self.xp1 = xp1
        self.xp2 = xp2
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.config = config
        self.theme = config["user_settings"]["theme"]

        self.selected_wave = None
        self.selected_flux = None
        self.fit_result = None

        self.init_plot()

    def init_plot(self):
        self.fig = plt.Figure(figsize=(10, 7))
        gs = GridSpec(2, 2, height_ratios=[2, 1], figure=self.fig)
        self.ax1 = self.fig.add_subplot(gs[0, :])
        self.ax2 = self.fig.add_subplot(gs[1, 0])
        self.ax3 = self.fig.add_subplot(gs[1, 1])

        self.ax1.set_title("Full Spectrum")
        self.ax2.set_title("Line Zoom")
        self.ax3.set_title("Population Diagram")

        self.ax1.plot(self.wave_data, self.flux_data, color=self.theme["foreground"])
        self.ax1.set_xlim(self.xp1, self.xp2)
        self.ax1.set_ylabel("Flux")

        # SpanSelector for selecting lines
        self.span = SpanSelector(self.ax1, self.onselect, 'horizontal', useblit=True,
                                 props=dict(alpha=0.5, facecolor=self.theme["selection_color"]))

    def onselect(self, xmin, xmax):
        mask = (self.wave_data >= xmin) & (self.wave_data <= xmax)
        if np.any(mask):
            self.selected_wave = self.wave_data[mask]
            self.selected_flux = self.flux_data[mask]
            self.update_zoom_plot()
            self.update_population_diagram()
        else:
            self.ax2.clear()
            self.ax2.set_title("Line Zoom")
            self.ax3.clear()
            self.ax3.set_title("Population Diagram")
        self.fig.canvas.draw_idle()

    def update_zoom_plot(self):
        self.ax2.clear()
        self.ax2.plot(self.selected_wave, self.selected_flux, color=self.theme["foreground"])
        self.ax2.set_title("Line Zoom")
        self.ax2.set_xlabel("Wavelength")
        self.ax2.set_ylabel("Flux")
        self.fit_gaussian()

    def fit_gaussian(self):
        model = GaussianModel()
        params = model.guess(self.selected_flux, x=self.selected_wave)
        self.fit_result = model.fit(self.selected_flux, params, x=self.selected_wave)
        self.ax2.plot(self.selected_wave, self.fit_result.best_fit, 'r--')

    def update_population_diagram(self):
        # For demonstration: fake data. Replace with real molecular energy levels from your molecules_data
        energy_levels = np.linspace(0, 5000, 10)
        populations = np.exp(-energy_levels / 1500)
        self.ax3.clear()
        self.ax3.scatter(energy_levels, populations, color=self.theme["foreground"])
        self.ax3.set_title("Population Diagram")
        self.ax3.set_xlabel("Energy Level (K)")
        self.ax3.set_ylabel("Relative Population")

    def embed(self, frame):
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.draw()