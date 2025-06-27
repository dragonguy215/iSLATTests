import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lmfit.models import GaussianModel
import numpy as np

class iSLATPlot:
    def __init__(self, window, input_spectrum_data, molecules_data, xp1, xp2, wave_data, flux_data, config):
        self.window = window
        self.input_spectrum_data = input_spectrum_data
        self.molecules_data = molecules_data
        self.xp1 = xp1
        self.xp2 = xp2
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.config = config

        self.selected_line = False
        self.line2save = None
        self.linesavepath = "LINESAVES/saved_lines.csv"  # or from config

        self.init_plot()

    def init_plot(self):
        user_settings = self.config['user_settings']
        self.fig = plt.Figure(figsize=(10, 6))
        gs = GridSpec(2, 1, height_ratios=[1, 1.5], figure=self.fig)
        self.ax1 = self.fig.add_subplot(gs[0])
        self.ax2 = self.fig.add_subplot(gs[1])

        self.ax1.plot(self.wave_data, self.flux_data, color=user_settings["theme"]["foreground"])
        self.ax1.set_xlim(self.xp1, self.xp2)
        self.ax1.set_ylabel("Flux density (Jy)")

        self.ax2.set_xlabel("Wavelength (μm)")
        self.ax2.set_ylabel("Flux density (Jy)")

    def embed(self, frame):
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        self.canvas.draw()

    def fit_line_selected(self):
        # Dummy: fits a gaussian on all data. Replace with your selection logic
        model = GaussianModel()
        params = model.guess(self.flux_data, x=self.wave_data)
        result = model.fit(self.flux_data, params, x=self.wave_data)
        center = np.round(result.params['center'].value, 5)
        fwhm = np.round(result.params['fwhm'].value, 2)
        area = np.round(result.params['amplitude'].value, 2)
        return ( (center, 0.01), (fwhm, 0.1), (area, 0.5) )

    def find_singles(self):
        # Dummy: finds peaks above threshold
        peaks = self.wave_data[self.flux_data > np.mean(self.flux_data) + 2*np.std(self.flux_data)]
        return peaks