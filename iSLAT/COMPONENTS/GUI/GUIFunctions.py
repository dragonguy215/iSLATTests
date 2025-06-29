import numpy as np
import tkinter as tk
from .Tooltips import CreateToolTip

def create_button(frame, theme, text, command, row, column):
        btn_theme = theme["buttons"].get(
            text.replace(" ", ""), theme["buttons"]["DefaultBotton"]
        )
        btn = tk.Button(
            frame, text=text,
            bg=btn_theme["background"],
            fg=theme["foreground"],
            activebackground=btn_theme["active_background"],
            command=command
        )
        btn.grid(row=row, column=column, padx=2, pady=2, sticky="nsew")
        CreateToolTip(btn, f"{text} button")
        return btn

class GUIHandlers:
    """
    Handles GUI button actions, connects to MainPlot for visualization and
    calls into iSLATClass for actual data processing (fits, CSV saves, slab fitting).
    """
    def __init__(self, main_plot, data_field, config, islat_class):
        self.main_plot = main_plot
        self.data_field = data_field
        self.config = config
        self.islat_class = islat_class
        self.user_settings = config
        self.legend_on = True

    def save_line(self):
        if self.main_plot.selected_wave is None:
            self.data_field.insert_text("No line selected to save.\n")
            return
        line_info = {
            "wavelength": np.mean(self.main_plot.selected_wave),
            "flux": np.max(self.main_plot.selected_flux)
        }
        self.islat_class.save_line(line_info)
        self.data_field.insert_text(f"Saved line at ~{line_info['wavelength']:.4f} μm\n")

    def fit_selected_line(self, deblend=False):
        self.data_field.clear()
        if self.main_plot.selected_wave is None:
            self.data_field.insert_text("No region selected for fitting.\n")
            return
        xmin, xmax = np.min(self.main_plot.selected_wave), np.max(self.main_plot.selected_wave)
        result = self.islat_class.fit_selected_line(xmin, xmax, deblend=deblend)
        if result:
            self.main_plot.update_zoom_line(self.main_plot.selected_wave, self.main_plot.selected_flux, result)
            self.data_field.insert_text(result.fit_report())
        else:
            self.data_field.insert_text("Fit failed or insufficient data.\n")

    def find_single_lines(self):
        lines = self.islat_class.find_single_lines()
        if lines:
            self.data_field.insert_text(f"Found {len(lines)} isolated lines.\n")
        else:
            self.data_field.insert_text("No lines found.\n")

    def single_slab_fit(self):
        result_text = self.islat_class.run_single_slab_fit()
        self.data_field.insert_text(f"Slab fit results:\n{result_text}\n")

    def toggle_legend(self):
        self.legend_on = not self.legend_on
        for ax in [self.main_plot.ax1, self.main_plot.ax2, self.main_plot.ax3]:
            if self.legend_on:
                ax.legend()
            else:
                leg = ax.get_legend()
                if leg: leg.set_visible(False)
        self.main_plot.canvas.draw_idle()

    def export_models(self):
        self.data_field.insert_text("Exporting current models...\n")
        if not self.islat_class.slab_model:
            self.data_field.insert_text("No slab model yet, running slab fit first...\n")
            self.islat_class.run_single_slab_fit()

        try:
            out_files = self.islat_class.slab_model.export_results()
            for f in out_files:
                self.data_field.insert_text(f"Exported to: {f}\n")
        except Exception as e:
            self.data_field.insert_text(f"Error exporting models: {e}\n")