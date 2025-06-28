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

        self.ax1.set_title("Full Spectrum with Line Inspection")
        self.ax1.set_ylabel("Flux (Jy)")
        self.ax1.set_xlabel("Wavelength (μm)")
        self.line1, = self.ax1.plot(self.wave_data, self.flux_data, color='black', label='Observed Spectrum')
        self.ax1.legend()

        # Bottom left - line inspection
        self.ax2.set_title("Line inspection plot")
        self.ax2.set_xlabel("Wavelength (μm)")
        self.ax2.set_ylabel("Flux (Jy)")

        # Bottom right - population diagram
        self.ax3.set_title("Population diagram")
        self.ax3.set_xlabel(r"$E_u$ (K)")
        self.ax3.set_ylabel(r"$\ln(N_u/g_u)$")

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent_frame)
        self.toolbar.pack(side="top", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.draw()

        # SpanSelector
        self.span = SpanSelector(
            self.ax1, self.onselect, "horizontal",
            useblit=True, interactive=True,
            props=dict(alpha=0.5, facecolor=self.theme.get("selection_color", "green"))
        )

        self.selected_wave = None
        self.selected_flux = None
        self.fit_result = None

    def onselect(self, xmin, xmax):
        print(f"Selected range: {xmin:.3f} - {xmax:.3f}")
        mask = (self.wave_data >= xmin) & (self.wave_data <= xmax)
        self.selected_wave = self.wave_data[mask]
        self.selected_flux = self.flux_data[mask]

        if len(self.selected_wave) < 5:
            self.ax2.clear()
            self.ax2.set_title("Line inspection plot")
            self.canvas.draw_idle()
            return

        self.fit_result = self.islat.fit_selected_line(xmin, xmax)
        self.update_zoom_line(self.selected_wave, self.selected_flux, self.fit_result)

        # Try to also update population diagram
        e, p = self.islat.generate_population_diagram()
        if e is not None:
            self.update_population_diagram(e, p)

    def update_zoom_line(self, wave, flux, fit_result=None):
        self.ax2.clear()
        self.ax2.plot(wave, flux, color='black', label="Data")
        if fit_result is not None:
            self.ax2.plot(wave, fit_result.best_fit, 'r--', label="Fit")
            self.ax2.legend()
        self.ax2.set_title("Line inspection plot")
        self.ax2.set_xlabel("Wavelength (μm)")
        self.ax2.set_ylabel("Flux (Jy)")
        self.canvas.draw_idle()

    def update_population_diagram(self, energies, pops):
        self.ax3.clear()
        self.ax3.scatter(energies, pops, color='black')
        self.ax3.set_title("Population diagram")
        self.ax3.set_xlabel(r"$E_u$ (K)")
        self.ax3.set_ylabel(r"$\ln(N_u/g_u)$")
        self.canvas.draw_idle()

    def toggle_legend(self):
        leg = self.ax1.get_legend()
        if leg is not None:
            vis = not leg.get_visible()
            leg.set_visible(vis)
            self.canvas.draw_idle()