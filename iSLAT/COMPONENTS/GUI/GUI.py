import tkinter as tk
from tkinter import ttk
from .MainPlot import iSLATPlot
#from MainPlot import iSLATPlot

class GUI:
    def __init__(self, master, isotopologue_data, data_field, mols, basem, isot, xp1, xp2, config):
        self.master = master
        self.isotopologue_data = isotopologue_data
        self.data_field = data_field
        self.mols = mols
        self.basem = basem
        self.isot = isot
        self.config = config
        #self.MainPlot = iSLATPlot(master, isotopologue_data, data_field, mols, basem, isot, config)
        self.MainPlot = iSLATPlot(master, isotopologue_data, data_field, xp1, xp2, data_field, data_field, config)

        # Create the main window for downloading HITRAN data
        self.create_window()

    def create_window(self):
        self.window = tk.Toplevel(self.master)

    def start(self):
        # Start the main plot
        self.MainPlot.create_plot()

        # Start the Tkinter main loop
        self.window.mainloop()


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    gui = GUI(root, None, None, None, None, None, None, None, config = {"iSLAT_version": "2.0", "user_settings":{'first_startup': False, 'reload_default_files': False, 'theme': {'theme_name': 'Light Theme', 'description': 'A light theme for iSLAT with a white background and black text.', 'foreground': '#000000', 'background': '#FFFFFF', 'toolbar': '#000000', 'graph_fill_color': '#D3D3D3', 'selection_color': '#00FF00', 'uncertainty_band_color': '#ABABAB', 'toolbar_highlight_background': '#000000', 'toolbar_highlight_color': '#000000', 'toolbar_active_color': '#FFFFFF', 'buttons': {'DefaultBotton': {'background': '#D3D3D3', 'active_background': '#808080', 'description': 'If a button is not specified, this will be used as the default button style.'}, 'ToggleLegend': {'background': '#D3D3D3', 'active_background': '#808080'}}}}})
    gui.start()