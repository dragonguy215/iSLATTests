import numpy as np

class GUIHandlers:
    """
    Class to handle GUI button actions. Keeps references to plot and data field.
    """
    def __init__(self, main_plot, data_field, config):
        self.main_plot = main_plot
        self.data_field = data_field
        self.config = config
        self.user_settings = config["user_settings"]

    def save_line(self):
        if not self.main_plot.selected_line:
            self.data_field.insert_text("No Line Selected!")
            return
        try:
            line2save = self.main_plot.line2save
            line2save.to_csv(self.main_plot.linesavepath, mode='a', index=False, header=False)
            self.data_field.insert_text("Line Saved!")
        except Exception as e:
            self.data_field.insert_text(f"Error saving line: {e}")

    def fit_selected_line(self):
        self.data_field.clear()
        if not self.main_plot.selected_line:
            self.data_field.insert_text("No Line Selected!")
            return
        fit_result = self.main_plot.fit_line_selected()
        if fit_result:
            centroid, fwhm, area = fit_result
            self.data_field.insert_text(f"Gaussian fit results:\n"
                f"Centroid (μm) = {centroid[0]} ± {centroid[1]}\n"
                f"FWHM (km/s) = {fwhm[0]} ± {fwhm[1]}\n"
                f"Area (erg/s/cm²) = {area[0]} ± {area[1]}")
        else:
            self.data_field.insert_text("Fit failed or incomplete.")

    def find_single_lines(self):
        self.data_field.insert_text("Finding single lines...")
        singles = self.main_plot.find_singles()
        self.data_field.insert_text(f"Found {len(singles)} isolated lines.")

    def single_slab_fit(self):
        self.data_field.insert_text("Running single slab fit (placeholder).")
        # Example of connecting to further fit routines