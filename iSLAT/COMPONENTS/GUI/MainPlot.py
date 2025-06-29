import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector
import numpy as np
from iSLAT.ir_model import Spectrum
from iSLAT.iSLATDefaultInputParms import dist, au, pc, ccum, hh

class iSLATPlot:
    def __init__(self, parent_frame, wave_data, flux_data, theme, islat_class_ref):
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.theme = theme
        self.islat = islat_class_ref

        self.fig = plt.Figure(figsize=(10, 7))
        gs = GridSpec(2, 2, height_ratios=[2, 1], figure=self.fig)
        self.ax1 = self.full_spectrum = self.fig.add_subplot(gs[0, :])
        self.ax2 = self.line_inspection = self.fig.add_subplot(gs[1, 0])
        self.ax3 = self.population_diagram = self.fig.add_subplot(gs[1, 1])

        self.ax1.set_title("Full Spectrum with Line Inspection")
        self.ax2.set_title("Line inspection plot")
        self.ax3.set_title("Population diagram")

        # Set default zoom similar to old behavior
        wmin, wmax = np.min(self.wave_data), np.max(self.wave_data)
        self.ax1.set_xlim(wmin, wmin + 0.25*(wmax - wmin))

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

        self.model_lines = []

        self.plot_model_lines()  # Initial plot of model lines
        self.plot_data_line(self.wave_data, self.flux_data, label="Observed Spectrum", color=self.theme["foreground"])
        self.compute_sum_flux()
        self.islat.update_model_spectrum()
        self.update_population_diagram()

    def clear_model_lines(self):
        # remove previously plotted lines
        for line in self.model_lines:
            line.remove()
        self.model_lines.clear()
        self.canvas.draw_idle()

    def update_model_plot(self):
        self.ax1.clear()
        # plot the data line
        self.ax1.plot(self.wave_data, self.flux_data, color=self.theme["foreground"], label="Data")

        # plot each molecule if it is turned on
        for mol in self.islat.molecules.values():
            if mol.is_active:
                model_flux = mol.get_flux(self.wave_data)
                self.ax1.fill_between(self.wave_data, model_flux, alpha=0.3, color=mol.color, label=mol.name)

        # plot sum of models if exists
        if hasattr(self.islat, 'sum_spectrum_flux'):
            self.ax1.plot(self.wave_data, self.islat.sum_spectrum_flux, 'r--', label='Sum of models')

        self.ax1.legend()
        self.canvas.draw_idle()

    def add_model_line(self, mol_name, temp, radius, density, color = None):
        """
        Adds a model spectrum line for given molecule parameters to the main plot.
        Assumes the islat.molecules dict has the MolData instances for calc.
        """
        molecule_obj = self.islat.molecules_dict[mol_name]
        model_flux = molecule_obj.spectrum.flux_jy

        if color is None:
            color = molecule_obj.line_color if hasattr(molecule_obj, "line_color") else self.theme["default_molecule_colors"][len(self.model_lines) % len(self.theme["default_molecule_colors"])]
        line, = self.ax1.plot(molecule_obj.spectrum.lamgrid, model_flux, linestyle='-', color=color, alpha=0.7, label=f"{mol_name}")

        self.model_lines.append(line)
        self.ax1.legend()
        self.canvas.draw_idle()

    def plot_model_lines(self):
        """
        Plots all model lines for the molecules in the islat.molecules_dict.
        Assumes the islat.molecules_dict has the MolData instances for calc.
        """
        self.clear_model_lines()
        for mol_name, molecule_obj in self.islat.molecules_dict.items():
            if molecule_obj.is_active:
                self.add_model_line(
                    mol_name,
                    temp=molecule_obj.temp,
                    radius=molecule_obj.radius,
                    density=molecule_obj.n_mol_init,
                )

    def plot_data_line(self, wave, flux, label=None, color=None):
        """
        Plots a data line on the main plot.
        """
        if label is None:
            label = "Data Line"
        if color is None:
            color = self.theme["foreground"]
        
        print("Plotting data line with wavelength and flux:")
        print("Wavelength:", wave)
        print("Flux:", flux)
        line, = self.ax1.plot(wave, flux, linestyle='-', color=color, alpha=0.7, label=label)
        self.model_lines.append(line)
        self.ax1.legend()
        self.canvas.draw_idle()
    
    def compute_sum_flux(self):
        """
        Computes the sum of all model fluxes and updates the sum line.
        """
        summed_flux = np.zeros_like(self.wave_data)
        for mol in self.islat.molecules_dict.values():
            if mol.is_active:
                mol_flux = mol.get_flux(self.wave_data)
                summed_flux += mol_flux
        return summed_flux

    def plot_sum_line(self, wave, flux, label=None, color=None, compute = True):
        """
        Plots the sum line on the main plot.
        """
        if compute:
            flux = self.compute_sum_flux()
        if label is None:
            label = "Sum Line"
        if color is None:
            color = self.theme["highlight"]
        
        line, = self.ax1.plot(wave, flux, linestyle='--', color=color, alpha=0.7, label=label)
        self.model_lines.append(line)
        self.ax1.legend()
        self.canvas.draw_idle()

    '''def draw_plot(self):
        self.ax1.clear()
        self.ax1.set_title("Full Spectrum")
        self.ax1.set_ylabel("Flux")

        # Plot data
        #self.ax1.plot(self.wave_data, self.flux_data, color=self.theme["foreground"], label="Data")
        self.plot_data_line(self.wave_data, self.flux_data, label="Observed Spectrum", color=self.theme["foreground"])

        # Plot sum line
        if hasattr(self.islat, 'sum_spectrum_flux'):
            self.ax1.plot(self.wave_data, self.islat.sum_spectrum_flux, color="orange", label="Sum Model")

        # Plot each molecule
        for mol in self.islat.molecules_dict.values():
            if mol.is_active:
                flux = mol.get_flux(self.wave_data)
                color = self.theme["default_molecule_colors"][len(self.model_lines) % len(self.theme["default_molecule_colors"])]
                self.ax1.plot(self.wave_data, flux, color=color, alpha=0.6, lw=1)
                self.ax1.fill_between(self.wave_data, flux, color=color, alpha=0.3, lw=0)

        self.ax1.legend()
        self.canvas.draw_idle()'''

    def onselect(self, xmin, xmax):
        mask = (self.wave_data >= xmin) & (self.wave_data <= xmax)
        self.selected_wave = self.wave_data[mask]
        self.selected_flux = self.flux_data[mask]
        self.islat.selected_wave = self.selected_wave
        self.islat.selected_flux = self.selected_flux

        if len(self.selected_wave) < 5:
            self.ax2.clear()
            self.canvas.draw_idle()
            return

        # Fit and update line
        self.fit_result = self.islat.fit_selected_line(xmin, xmax)
        #self.update_line_inspection_plot(self.selected_wave, self.selected_flux, self.fit_result)
        self.update_line_inspection_plot()

        self.update_population_diagram()

    def update_line_inspection_plot(self):
        self.ax2.clear()

        # Calculate indices of selection in global wavelength grid
        wave = self.wave_data
        model_sum = self.islat.sum_spectrum_flux
        data_flux = self.flux_data

        if hasattr(self.islat, 'selected_wave'):
            xmin = np.min(self.islat.selected_wave)
            xmax = np.max(self.islat.selected_wave)

            # Find indices
            idx_min = np.searchsorted(wave, xmin) - 1
            idx_max = np.searchsorted(wave, xmax) + 1
            idx_min = max(idx_min, 0)
            idx_max = min(idx_max, len(wave)-1)

            # Slice
            model_region_x = wave[idx_min:idx_max]
            model_region_y = model_sum[idx_min:idx_max]
            data_region_y = data_flux[idx_min:idx_max]

            # Determine max y
            max_y = max(np.nanmax(model_region_y), np.nanmax(data_region_y))
            self.ax2.set_ylim(0, max_y)

            # Plot summed model
            self.ax2.plot(model_region_x, model_region_y, 'k--', lw=2, label="Sum model")

            # Plot data
            self.ax2.plot(model_region_x, data_region_y, 'r-', lw=1.2, label="Data")

            # Optional: fill by molecule contributions
            for mol in self.islat.molecules_dict.values():
                if mol.is_active:
                    mol_flux = mol.spectrum.flux_jy[idx_min:idx_max]
                    color = mol.color
                    self.ax2.fill_between(model_region_x, mol_flux, alpha=0.4, color=color, label=mol.displaylabel)

            self.ax2.set_xlim(model_region_x[0], model_region_x[-1])
            self.ax2.set_xlabel("Wavelength (μm)")
            self.ax2.set_ylabel("Flux")
            self.ax2.legend()
            self.canvas.draw_idle()

    def update_population_diagram(self):
        self.ax3.clear()
        self.ax3.set_ylabel(r'ln(4πF/(hν$A_{u}$$g_{u}$))')
        self.ax3.set_xlabel(r'$E_{u}$ (K)')
        self.ax3.set_title('Population diagram', fontsize='medium')

        molecule_obj = self.islat.molecules_dict["H2O"]
        int_pars = molecule_obj.intensity.get_table
        int_pars.index = range(len(int_pars.index))

        # Parsing the components of the lines in int_pars
        wl = int_pars['lam']
        intens_mod = int_pars['intens']
        Astein_mod = int_pars['a_stein']
        gu = int_pars['g_up']
        eu = int_pars['e_up']

        # Calculating the y-axis for the population diagram for each line in int_pars
        area = np.pi * (molecule_obj.radius * au * 1e2) ** 2  # In cm^2
        Dist = dist * pc
        beam_s = area / Dist ** 2
        F = intens_mod * beam_s
        freq = ccum / wl
        rd_yax = np.log(4 * np.pi * F / (Astein_mod * hh * freq * gu))
        threshold = np.nanmax(F) / 100

        self.ax3.set_ylim(np.nanmin(rd_yax[F > threshold]), np.nanmax(rd_yax) + 0.5)
        self.ax3.set_xlim(np.nanmin(eu) - 50, np.nanmax(eu[F > threshold]))

        # Populating the population diagram graph with the lines
        self.ax3.scatter(eu, rd_yax, s=0.5, color='#838B8B')
        self.canvas.draw_idle()

    def toggle_legend(self):
        leg = self.ax1.get_legend()
        if leg is not None:
            vis = not leg.get_visible()
            leg.set_visible(vis)
            self.canvas.draw_idle()