import tkinter as tk
from tkinter import filedialog
from .MainPlot import iSLATPlot
from .Data_field import DataField
from .MoleculeWindow import MoleculeWindow
#from .Tooltips import CreateToolTip
#from .GUIFunctions import GUIHandlers
from .ControlPanel import ControlPanel
from .TopOptions import TopOptions
from .BottomOptions import BottomOptions
import os

class GUI:
    def __init__(self, master, molecule_data, wave_data, flux_data, config, islat_class_ref):
        self.master = master
        self.molecule_data = molecule_data
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.config = config
        self.theme = config["theme"]
        self.islat_class = islat_class_ref

    def build_left_panel(self, parent):
        # Top control buttons
        self.top_options = TopOptions(parent, self.islat_class, theme=self.theme)
        self.top_options.frame.pack(fill="x")

        # Molecule table
        self.molecule_table = MoleculeWindow("Molecule Table", parent, self.molecule_data, self.plot, self.config, self.islat_class)
        self.molecule_table.frame.pack(fill="both", expand=True, padx=5, pady=5)

        # Spectrum file selector
        file_frame = tk.LabelFrame(parent, text="Spectrum File")
        file_frame.pack(fill="x", padx=5, pady=5)
        self.file_label = tk.Label(file_frame, text="Loaded: File")
        self.file_label.pack()
        tk.Button(file_frame, text="Load Spectrum", command=self.islat_class.load_spectrum).pack()

        # Control panel for input parameters
        control_panel_frame = tk.LabelFrame(parent, text="Control Panel")
        control_panel_frame.pack(fill="both", expand=True, padx=5, pady=5)
        self.control_panel = ControlPanel(control_panel_frame, self.islat_class)

        # Main data field
        self.data_field = DataField("Main Data Field", "", parent)
        self.data_field.frame.pack(fill="both", expand=True)

    def create_window(self):
        self.window = self.master
        self.window.title("iSLAT Version 5.00.00")
        self.window.columnconfigure(0, weight=1)
        self.window.columnconfigure(1, weight=3)
        self.window.rowconfigure(0, weight=1)
        self.window.rowconfigure(1, weight=0)

        # Right side: plots
        right_frame = tk.Frame(self.window)
        right_frame.grid(row=0, column=1, sticky="nsew")
        self.plot = iSLATPlot(right_frame, self.wave_data, self.flux_data, self.theme, self.islat_class)

        # Left side: all controls
        left_frame = tk.Frame(self.window)
        left_frame.grid(row=0, column=0, rowspan=2, sticky="nsew")
        self.build_left_panel(left_frame)

        # Bottom function buttons
        self.bottom_options = BottomOptions(self.window, self.islat_class, self.theme, self.plot, self.data_field, self.config)
        self.bottom_options.frame.grid(row=1, column=0, columnspan=2, sticky="ew")

    def start(self):
        self.create_window()
        self.window.mainloop()

    '''def load_spectrum_file(self):
        file_path = filedialog.askopenfilename(title="Select spectrum CSV")
        if file_path:
            self.islat_class.load_spectrum(file_path)
            self.file_label.config(text=f"Loaded: {os.path.basename(file_path)}")
            self.plot.wave_data = self.islat_class.wave_data
            self.plot.flux_data = self.islat_class.flux_data
            self.plot.ax1.clear()
            self.plot.ax1.plot(self.plot.wave_data, self.plot.flux_data, color=self.theme["foreground"])
            self.plot.canvas.draw_idle()'''