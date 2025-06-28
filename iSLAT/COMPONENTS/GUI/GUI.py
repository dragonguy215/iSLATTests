import tkinter as tk
from tkinter import filedialog
from .MainPlot import iSLATPlot
from .Data_field import DataField
from .MoleculeWindow import MoleculeWindow
from .Tooltips import CreateToolTip
from .GUIFunctions import GUIHandlers
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

    def create_button(self, frame, text, command, row, column):
        btn_theme = self.theme["buttons"].get(
            text.replace(" ", ""), self.theme["buttons"]["DefaultBotton"]
        )
        btn = tk.Button(
            frame, text=text,
            bg=btn_theme["background"],
            fg=self.theme["foreground"],
            activebackground=btn_theme["active_background"],
            command=command
        )
        btn.grid(row=row, column=column, padx=2, pady=2, sticky="nsew")
        CreateToolTip(btn, f"{text} button")
        return btn

    def build_left_panel(self, parent):
        # Configure the left panel
        control_frame = tk.Frame(parent)
        control_frame.pack(fill="x")
        self.create_button(control_frame, "Default Molecules", self.default_molecules, 0, 0)
        self.create_button(control_frame, "Load Parameters", self.load_parameters, 0, 1)
        self.create_button(control_frame, "Save Parameters", self.save_parameters, 0, 2)
        self.create_button(control_frame, "HITRAN Query", self.hitran_query, 1, 0)
        self.create_button(control_frame, "Export Models", self.export_models, 1, 1)
        self.create_button(control_frame, "Toggle Legend", self.toggle_legend, 1, 2)

        # Spectrum file selector
        file_frame = tk.LabelFrame(parent, text="Spectrum File")
        file_frame.pack(fill="x", padx=5, pady=5)
        self.file_label = tk.Label(file_frame, text="Loaded: File")
        self.file_label.pack()
        tk.Button(file_frame, text="Load Spectrum", command=self.load_spectrum_file).pack()

        # Molecule table
        self.molecule_table = MoleculeWindow("Molecule Table", parent, self.molecule_data, self.config, self.islat_class)
        self.molecule_table.table_frame.pack(fill="both", expand=True, padx=5, pady=5)

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

        # Left side: all controls
        left_frame = tk.Frame(self.window)
        left_frame.grid(row=0, column=0, rowspan=2, sticky="nsew")
        self.build_left_panel(left_frame)

        # Right side: plots
        right_frame = tk.Frame(self.window)
        right_frame.grid(row=0, column=1, sticky="nsew")
        self.plot = iSLATPlot(right_frame, self.wave_data, self.flux_data, self.theme, self.islat_class)

        # Bottom function buttons
        func_frame = tk.Frame(self.window)
        func_frame.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.handlers = GUIHandlers(self.plot, self.data_field, self.config, self.islat_class)
        self.create_button(func_frame, "Save Line", self.handlers.save_line, 0, 0)
        self.create_button(func_frame, "Fit Line", self.handlers.fit_selected_line, 0, 1)
        self.create_button(func_frame, "Find Single Lines", self.handlers.find_single_lines, 0, 2)
        self.create_button(func_frame, "Line De-blender", lambda: self.handlers.fit_selected_line(deblend=True), 0, 3)
        self.create_button(func_frame, "Single Slab Fit", self.handlers.single_slab_fit, 0, 4)

    def start(self):
        self.create_window()
        self.window.mainloop()

    # Callbacks for top buttons
    def default_molecules(self):
        print("Default molecules loaded")

    def load_parameters(self):
        print("Load parameters from file")

    def save_parameters(self):
        print("Save parameters to file")

    def hitran_query(self):
        print("Perform HITRAN query")

    def export_models(self):
        print("Export models to file")

    def toggle_legend(self):
        print("Toggled legend on plot")
        self.plot.toggle_legend()

    def load_spectrum_file(self):
        file_path = filedialog.askopenfilename(title="Select spectrum CSV")
        if file_path:
            self.islat_class.load_spectrum(file_path)
            self.file_label.config(text=f"Loaded: {os.path.basename(file_path)}")
            self.plot.wave_data = self.islat_class.wave_data
            self.plot.flux_data = self.islat_class.flux_data
            self.plot.ax1.clear()
            self.plot.ax1.plot(self.plot.wave_data, self.plot.flux_data, color=self.theme["foreground"])
            self.plot.canvas.draw_idle()