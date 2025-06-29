import numpy as np
import tkinter as tk
#import pandas as pd
#import os
#from tkinter import ttk
import iSLAT.iSLATFileHandling as ifh
from .GUIFunctions import create_button

class BottomOptions:
    def __init__(self, master, islat, theme, main_plot, data_field, config):
        self.master = master
        self.islat = islat
        self.theme = theme
        self.main_plot = main_plot
        self.data_field = data_field
        self.config = config
        #self.theme = self.master.theme

        # Create the frame for top options
        self.frame = tk.Frame(master, borderwidth=2, relief="groove")
        #self.frame.grid(row=1, column=0, columnspan=2, sticky="ew")

        # Create buttons for top options
        create_button(self.frame, self.theme, "Save Line", self.save_line, 0, 0)
        create_button(self.frame, self.theme, "Show Saved Lines", self.show_saved_lines, 0, 1)
        create_button(self.frame, self.theme, "Fit Line", self.fit_selected_line, 0, 2)
        create_button(self.frame, self.theme, "Fit Saved Lines", self.fit_saved_lines, 0, 3)
        create_button(self.frame, self.theme, "Find Single Lines", self.find_single_lines, 0, 4)
        create_button(self.frame, self.theme, "Line De-blender", lambda: self.fit_selected_line(deblend=True), 0, 5)
        create_button(self.frame, self.theme, "Single Slab Fit", self.single_slab_fit, 0, 6)
    
    def save_line(self):
        if self.main_plot.selected_wave is None:
            self.data_field.insert_text("No line selected to save.\n")
            return
        line_info = {
            #"wavelength": np.mean(self.main_plot.selected_wave),
            #"flux": np.max(self.main_plot.selected_flux)
            "species": self.islat.active_molecule.name,
            "lev_up": self.islat.active_molecule.name,
            "lev_low": self.islat.active_molecule.name,
            "lam": np.mean(self.main_plot.selected_wave),
            "tau": np.max(self.main_plot.selected_flux),
            "intens": np.max(self.main_plot.selected_flux) * self.main_plot.ax1.get_xlim()[1] * 1e-6,  # Convert to correct units
            "a_stein": self.islat.active_molecule.name,
            "e_up": self.islat.active_molecule.name,
            "g_up": self.islat.active_molecule.name,
            "xmin": self.main_plot.ax1.get_xlim()[0],
            "xmax": self.main_plot.ax1.get_xlim()[1],
            #"lev_up": self.islat.molecules_dict[self.islat.active_molecule].lev_up,
            #"lev_low": self.islat.molecules_dict[self.islat.active_molecule].lev_low,
            #"lam": np.mean(self.main_plot.selected_wave),
            #"tau": np.max(self.main_plot.selected_flux),
            #"intens": np.max(self.main_plot.selected_flux) * self.main_plot.ax1.get_xlim()[1] * 1e-6,  # Convert to correct units
            #"a_stein": self.islat.molecules_dict[self.islat.active_molecule].a_stein,
            #"e_up": self.islat.molecules_dict[self.islat.active_molecule].e_up,
            #"g_up": self.islat.molecules_dict[self.islat.active_molecule].g_up,
            #xmin": self.main_plot.ax1.get_xlim()[0],
            #xmax": self.main_plot.ax1.get_xlim()[1],
        }
        ifh.save_line(line_info)
        self.data_field.insert_text(f"Saved line at ~{line_info['lam']:.4f} μm\n")

    def show_saved_lines(self): ####
        #self.data_field.clear()

        '''try:
            linelistfile = ifh.read_line_saves()
        except AttributeError:
            self.data_field.insert_text("Input line list is not defined!\n")
            return'''

        '''# Initialize variables
        self.main_plot.default_line = None
        self.main_plot.green_lines = []
        self.main_plot.green_scatter = []'''

        # Load saved lines
        try:
            svd_lns = ifh.read_line_saves()
        except Exception as e:
            self.data_field.insert_text(f"Error loading saved lines: {e}\n")
            return

        self.main_plot.plot_saved_lines(svd_lns)

        '''svd_lamb = np.array(svd_lns['lam'])
        x_min = np.array(svd_lns['xmin']) if 'xmin' in svd_lns else None
        x_max = np.array(svd_lns['xmax']) if 'xmax' in svd_lns else None

        # Plot vertical lines in the main plot for saved lines
        for i in range(len(svd_lamb)):
            self.main_plot.ax1.vlines(svd_lamb[i], -2, 10, linestyles='dashed', color='red')
            if x_min is not None and x_max is not None:
                self.main_plot.ax1.vlines(x_min[i], -2, 10, color='coral', alpha=0.5)
                self.main_plot.ax1.vlines(x_max[i], -2, 10, color='coral', alpha=0.5)

        self.data_field.insert_text("Saved lines retrieved from file.\n")
        self.main_plot.canvas.draw()'''

    def fit_selected_line(self, deblend=False):
        self.data_field.clear()
        if self.main_plot.selected_wave is None:
            self.data_field.insert_text("No region selected for fitting.\n")
            return

        self.main_plot.compute_fit_line(deblend=deblend)
        fit_result = self.main_plot.fit_result
        self.main_plot.update_line_inspection_plot()

        if fit_result:
            self.data_field.insert_text(fit_result)
        else:
            self.data_field.insert_text("Fit failed or insufficient data.\n")

    def fit_saved_lines(self): ########
        self.data_field.clear()
        saved_lines = ifh.read_line_saves()
        if not saved_lines:
            self.data_field.insert_text("No saved lines to fit.\n")
            return

        self.data_field.insert_text("Fitting saved lines...\n")
        fit_results = []
        for line in saved_lines:
            self.main_plot.selected_wave = [line['lam']]
            self.main_plot.selected_flux = [line['tau']]
            fit_result = self.main_plot.compute_fit_line()
            if fit_result:
                fit_results.append(fit_result)
                self.data_field.insert_text(f"Fit for line at {line['wavelength']:.4f} μm: {fit_result}\n")
            else:
                self.data_field.insert_text(f"Fit failed for line at {line['wavelength']:.4f} μm.\n")

        if fit_results:
            self.main_plot.update_line_inspection_plot()

    def find_single_lines(self):
        self.main_plot.find_single_lines()
        lines = self.main_plot.single_lines_list
        if lines:
            self.data_field.insert_text(f"Found {len(lines)} isolated lines.\n")
            self.main_plot.plot_single_lines()
            #for line in lines:
            #    self.data_field.insert_text(f"Line at {line['wavelength']:.4f} μm with flux {line['flux']:.4f}\n")
        else:
            self.data_field.insert_text("No isolated lines found.\n")

    def single_slab_fit(self):
        result_text = self.islat.run_single_slab_fit()
        self.data_field.insert_text(f"Slab fit results:\n{result_text}\n")

    def export_models(self):
        self.data_field.insert_text("Exporting current models...\n")
        if not self.islat.slab_model:
            self.data_field.insert_text("No slab model yet, running slab fit first...\n")
            self.islat.run_single_slab_fit()

        try:
            out_files = self.islat.slab_model.export_results()
            for f in out_files:
                self.data_field.insert_text(f"Exported to: {f}\n")
        except Exception as e:
            self.data_field.insert_text(f"Error exporting models: {e}\n")