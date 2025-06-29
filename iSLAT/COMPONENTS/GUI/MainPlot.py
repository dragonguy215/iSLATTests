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

        #self.ax1.plot(self.wave_data, self.flux_data, color=self.theme["foreground"], label="Observed Spectrum")
        #self.ax1.set_ylabel("Flux (Jy)")

        '''print("Here is the wave and flux data for the main plot:")
        print("Wavelength:", self.wave_data)
        print("Flux:", self.flux_data)'''

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
        self.plot_sum_line(self.wave_data, np.zeros_like(self.wave_data), label="Sum Line", color=self.theme["highlight"])
        self.islat.update_model_spectrum()
        #self.draw_plot()  # Initial draw of the main plot
        #print("ax1 lines:", self.ax1.lines)
        self.update_population_diagram()

    '''@property
    def model_lines(self):
        """
        Returns the list of model lines currently plotted on the main plot.
        """
        return self._model_lines
    
    @model_lines.setter
    def model_lines(self, value):
        """
        reloads the plot whenever model_lines is set.
        """
        self._model_lines = value
        self.canvas.'''

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
        #model_flux = mol_obj.calculate_spectrum(temp, radius, density, self.wave_data)  # Ensure the method name matches the actual implementation in MolData
        #model_flux = Spectrum(lam_min=self.wave_data[0], lam_max=self.wave_data[-1], dlambda=model_pixel_res, R=model_line_width, distance=dist).flux_jy
        #(self, lam_min=None, lam_max=None, dlambda=None, R=None, distance=None):
        model_flux = molecule_obj.spectrum.flux_jy

        '''if molecule_obj.name == "H2O":
            print("hey bro, here are the current values for h20:")
            print(f"temp: {molecule_obj.temp}, radius: {molecule_obj.radius}")
            print(f"n_mol: {molecule_obj.n_mol}, scale_exponent: {molecule_obj.scale_exponent}, scale_number: {molecule_obj.scale_number}")
            print(f"t_kin: {molecule_obj.t_kin}, intrinsic_line_width: {molecule_obj.intrinsic_line_width}, model_pixel_res: {molecule_obj.model_pixel_res}")
            print(f"model_line_width: {molecule_obj.model_line_width}, distance: {molecule_obj.distance}, wavelength_range: {molecule_obj.wavelength_range}")'''

        if color is None:
            color = molecule_obj.line_color if hasattr(molecule_obj, "line_color") else self.theme["default_molecule_colors"][len(self.model_lines) % len(self.theme["default_molecule_colors"])]
        #line, = self.ax1.plot(self.wave_data, model_flux, linestyle='-', color=color, alpha=0.7, label=f"{mol_name}")
        line, = self.ax1.plot(molecule_obj.spectrum.lamgrid, model_flux, linestyle='-', color=color, alpha=0.7, label=f"{mol_name}")
        '''print("Here is the lamgrid for the model line:")
        print(molecule_obj.spectrum.lamgrid)
        print("And here is the model flux:")
        print(model_flux)'''

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
    
    def plot_sum_line(self, wave, flux, label=None, color=None):
        """
        Plots a summed line on the main plot.
        """
        if label is None:
            label = "Sum Line"
        if color is None:
            color = self.theme["highlight"]
        
        line, = self.ax1.plot(wave, flux, linestyle='-', color=color, alpha=0.7, label=label)
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
            #self.ax2.set_title("Line inspection plot")
            self.canvas.draw_idle()
            return

        # Fit and update line
        self.fit_result = self.islat.fit_selected_line(xmin, xmax)
        self.update_zoom_line(self.selected_wave, self.selected_flux, self.fit_result)

        self.update_population_diagram()

    def update_zoom_line(self, wave, flux, fit_result=None):
        self.ax2.clear()
        self.ax2.plot(wave, flux, color=self.theme["foreground"], label="Data")
        if fit_result is not None:
            self.ax2.plot(wave, fit_result.best_fit, 'r--', label="Fit")
            try:
                conf = fit_result.eval_uncertainty(sigma=self.islat.user_settings["fit_line_uncertainty"])
                self.ax2.fill_between(wave, fit_result.best_fit - conf, fit_result.best_fit + conf,
                                      color=self.theme["fit_line_color"], alpha=0.3, label=f"±{self.islat.user_settings["fit_line_uncertainty"]}σ")
            except:
                pass
            self.ax2.legend()
        self.ax2.set_xlabel("Wavelength (μm)")
        self.ax2.set_ylabel("Flux (Jy)")
        self.ax2.set_title("Line inspection plot")
        self.canvas.draw_idle()

    def update_population_diagram(self):
        self.ax3.clear()
        self.ax3.set_ylabel(r'ln(4πF/(hν$A_{u}$$g_{u}$))')
        self.ax3.set_xlabel(r'$E_{u}$ (K)')
        self.ax3.set_title('Population diagram', fontsize='medium')

        #molecule_obj = self.islat.molecules_dict[mol_name]
        #print("Whats good, here is the molecules dict:")
        #print(self.islat.molecules_dict)
        molecule_obj = self.islat.molecules_dict["H2O"]
        #print("And here is the molecule object:")
        #print(molecule_obj)
        int_pars = molecule_obj.intensity.get_table
        #print("Hey man heres that table you wanted:")
        #print(int_pars)
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