import numpy as np
import tkinter as tk
from tkinter import ttk
from .GUIFunctions import create_button
from lmfit.models import GaussianModel

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
        create_button(self.frame, self.theme, "Fit Line", self.fit_selected_line, 0, 1)
        create_button(self.frame, self.theme, "Find Single Lines", self.find_single_lines, 0, 2)
        create_button(self.frame, self.theme, "Line De-blender", lambda: self.fit_selected_line(deblend=True), 0, 3)
        create_button(self.frame, self.theme, "Single Slab Fit", self.single_slab_fit, 0, 4)
    
    def save_line(self):
        if self.main_plot.selected_wave is None:
            self.data_field.insert_text("No line selected to save.\n")
            return
        line_info = {
            "wavelength": np.mean(self.main_plot.selected_wave),
            "flux": np.max(self.main_plot.selected_flux)
        }
        self.islat.save_line(line_info)
        self.data_field.insert_text(f"Saved line at ~{line_info['wavelength']:.4f} μm\n")

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

    def find_single_lines(self):
        lines = self.islat.find_single_lines()
        if lines:
            self.data_field.insert_text(f"Found {len(lines)} isolated lines.\n")
        else:
            self.data_field.insert_text("No lines found.\n")

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