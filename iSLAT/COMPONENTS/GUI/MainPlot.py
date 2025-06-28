import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector
import numpy as np

class iSLATPlot:
    def __init__(self, parent_frame, wave_data, flux_data, theme, islat_class_ref):
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.theme = theme
        self.islat = islat_class_ref

        self.fig = plt.Figure(figsize=(10, 7))
        gs = GridSpec(2, 2, height_ratios=[2, 1], figure=self.fig)
        self.ax1 = self.fig.add_subplot(gs[0, :])
        self.ax2 = self.fig.add_subplot(gs[1, 0])
        self.ax3 = self.fig.add_subplot(gs[1, 1])

        self.ax1.set_title("Full Spectrum")
        self.ax2.set_title("Line inspection plot")
        self.ax3.set_title("Population diagram")

        self.ax1.plot(self.wave_data, self.flux_data, color=self.theme["foreground"])
        self.ax1.set_ylabel("Flux")
        self.span = SpanSelector(
            self.ax1, self.onselect, "horizontal",
            useblit=True, interactive=True,
            props=dict(alpha=0.5, facecolor=self.theme["selection_color"])
        )

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent_frame)
        self.toolbar.pack(side="top", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.draw()

        self.selected_wave = None
        self.selected_flux = None
        self.fit_result = None

    def onselect(self, xmin, xmax):
        mask = (self.wave_data >= xmin) & (self.wave_data <= xmax)
        self.selected_wave = self.wave_data[mask]
        self.selected_flux = self.flux_data[mask]
        if len(self.selected_wave) < 5:
            self.ax2.clear()
            self.ax2.set_title("Line inspection plot")
            self.canvas.draw_idle()
            return

        # Perform the default single line fit to overlay on inspection
        self.fit_result = self.islat.fit_selected_line(xmin, xmax)
        self.update_zoom_line(self.selected_wave, self.selected_flux, self.fit_result)

        # Try to also update population diagram
        e, p = self.islat.generate_population_diagram()
        print("whats good home, heres the population diagram data, yo:", e, p)
        if e is not None:
            self.update_population_diagram(e, p)

    def update_zoom_line(self, wave, flux, fit_result=None):
        self.ax2.clear()
        self.ax2.plot(wave, flux, color=self.theme["foreground"], label="Data")
        if fit_result is not None:
            self.ax2.plot(wave, fit_result.best_fit, 'r--', label="Fit")
            self.ax2.legend()
        self.ax2.set_title("Line inspection plot")
        self.ax2.set_xlabel("Wavelength")
        self.ax2.set_ylabel("Flux")
        self.canvas.draw_idle()

    def update_population_diagram(self, energy_levels, populations):
        self.ax3.clear()
        self.ax3.scatter(energy_levels, populations, color=self.theme["foreground"])
        self.ax3.set_title("Population diagram")
        self.ax3.set_xlabel("E_upper (K)")
        self.ax3.set_ylabel("ln(Nu/gu)")
        self.canvas.draw_idle()