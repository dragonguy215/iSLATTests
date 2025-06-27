import tkinter as tk
from tkinter import ttk
from .MainPlot import iSLATPlot
from .Data_field import DataField
from .Tooltips import CreateToolTip
from .GUIFunctions import GUIHandlers

class GUI:
    def __init__(self, master, isotopologue_data, input_spectrum_data, wave_data, flux_data, mols, basem, isot, xp1, xp2, config):
        self.master = master
        self.isotopologue_data = isotopologue_data
        self.input_spectrum_data = input_spectrum_data
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.mols = mols
        self.basem = basem
        self.isot = isot
        self.xp1 = xp1
        self.xp2 = xp2
        self.config = config
        self.theme = config["user_settings"]["theme"]
        self.user_settings = config["user_settings"]

    def create_button(self, frame, text, command, row, column):
        btn_theme = self.theme["buttons"].get(
            text.replace(" ", ""), 
            self.theme["buttons"]["DefaultBotton"]
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

    def create_window(self):
        self.window = self.master
        self.window.title("iSLAT - Infrared Spectral Line Analysis Tool")
        self.window.columnconfigure(0, weight=3)
        self.window.columnconfigure(1, weight=1)
        self.window.rowconfigure(0, weight=1)

        # LEFT side (plot and molecule table)
        left_frame = tk.Frame(self.window)
        left_frame.grid(row=0, column=0, sticky="nsew")
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(0, weight=4)
        left_frame.rowconfigure(1, weight=2)

        # RIGHT side (data field)
        right_frame = tk.Frame(self.window)
        right_frame.grid(row=0, column=1, sticky="nsew")
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)

        # Plot area on top
        plot_frame = tk.Frame(left_frame, borderwidth=2, relief="sunken")
        plot_frame.grid(row=0, column=0, sticky="nsew")
        self.main_plot = iSLATPlot(self.window, self.input_spectrum_data, self.isotopologue_data, self.xp1, self.xp2, self.wave_data, self.flux_data, self.config)
        self.main_plot.embed(plot_frame)

        # Molecule table (bottom left)
        table_frame = tk.Frame(left_frame)
        table_frame.grid(row=1, column=0, sticky="nsew")
        self.build_molecule_table(table_frame)

        # Data field on right
        self.data_field = DataField("Main Data Field", "", right_frame)
        self.data_field.frame.grid(row=0, column=0, sticky="nsew")

    def build_molecule_table(self, parent):
        headers = ['Molecule', 'Temp.', 'Radius', 'Col. Dens', 'On', 'Del.', 'Color']
        for col, text in enumerate(headers):
            tk.Label(parent, text=text, bg=self.theme["background"], fg=self.theme["foreground"]).grid(row=0, column=col)

        # Example: add one row per molecule
        for i, mol in enumerate(self.mols):
            tk.Label(parent, text=mol, bg=self.theme["background"], fg=self.theme["foreground"]).grid(row=i+1, column=0)
            tk.Entry(parent).grid(row=i+1, column=1)
            tk.Entry(parent).grid(row=i+1, column=2)
            tk.Entry(parent).grid(row=i+1, column=3)
            tk.Checkbutton(parent).grid(row=i+1, column=4)
            tk.Checkbutton(parent).grid(row=i+1, column=5)
            tk.Button(parent, text="Pick").grid(row=i+1, column=6)

    def start(self):
        # Create function buttons after window is initialized
        self.create_window()

        func_frame = tk.Frame(self.window)
        func_frame.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.handlers = GUIHandlers(self.main_plot, self.data_field, self.config)

        self.create_button(func_frame, "Save Line", self.handlers.save_line, 0, 0)
        self.create_button(func_frame, "Fit Line", self.handlers.fit_selected_line, 0, 1)
        self.create_button(func_frame, "Find Single Lines", self.handlers.find_single_lines, 0, 2)
        self.create_button(func_frame, "Single Slab Fit", self.handlers.single_slab_fit, 0, 3)

        self.window.mainloop()