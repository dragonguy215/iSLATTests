import tkinter as tk
from tkinter import ttk
from .MainPlot import iSLATPlot
from .Data_field import DataField
from .Tooltips import CreateToolTip
from .GUIFunctions import GUIHandlers

class GUI:
    def __init__(self, master, isotopologue_data, input_spectrum_data, wave_data, flux_data, mols, basem, isot, xp1, xp2, config):
        self.master = self.window = master
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

        self.nb_of_columns = 10

        # Set up main plot and data field
        self.main_plot = iSLATPlot(master, input_spectrum_data, isotopologue_data, xp1, xp2, wave_data, flux_data, config)
        #self.data_field = DataField("Main Data Field", "", None)

        # Set up GUI function handlers with references
        #self.handlers = GUIHandlers(self.main_plot, self.data_field, self.config)

    def create_button(self, frame, text, command, row, column, width=13, height=1):
        button_theme = self.config["user_settings"]["theme"]["buttons"].get(
            text.replace(" ", ""), 
            self.config["user_settings"]["theme"]["buttons"]["DefaultBotton"]
        )
        button = tk.Button(
            frame,
            text=text,
            bg=button_theme["background"],
            activebackground=button_theme["active_background"],
            fg=self.config["user_settings"]["theme"]["foreground"],
            command=command,
            width=width,
            height=height
        )
        button.grid(row=row, column=column, padx=2, pady=2, sticky="nsew")
        CreateToolTip(button, f"{text} button")
        return button

    def create_window(self):
        #self.window = tk.Toplevel(self.master)
        #self.window.title("iSLAT GUI")

        # Configure window resizing
        self.window.columnconfigure(0, weight=3)
        self.window.columnconfigure(1, weight=1)
        self.window.rowconfigure(0, weight=1)

        # Left side frame (plot + functions)
        self.left_frame = tk.Frame(self.window)
        self.left_frame.grid(row=0, column=0, sticky="nsew")
        self.left_frame.columnconfigure(0, weight=1)
        self.left_frame.rowconfigure(0, weight=4)
        self.left_frame.rowconfigure(1, weight=1)

        # Right side frame (data field)
        self.right_frame = tk.Frame(self.window)
        self.right_frame.grid(row=0, column=1, sticky="nsew")
        self.right_frame.columnconfigure(0, weight=1)
        self.right_frame.rowconfigure(0, weight=1)

        # Plot area
        self.plot_frame = tk.Frame(self.left_frame, borderwidth=2, relief="sunken")
        self.plot_frame.grid(row=0, column=0, sticky="nsew")
        self.main_plot.embed(self.plot_frame)

        # Now create data field in right frame
        self.data_field = DataField("Main Data Field", "", self.right_frame)
        self.data_field.frame.grid(row=0, column=0, sticky="nsew")

        self.handlers = GUIHandlers(self.main_plot, self.data_field, self.config)

        # Functions area
        self.functions_frame = tk.Frame(self.left_frame, borderwidth=2, relief="groove")
        self.functions_frame.grid(row=1, column=0, sticky="nsew")
        self.functions_frame.columnconfigure([0,1], weight=1)

        self.create_button(self.functions_frame, "Save Line", self.handlers.save_line, 0, 0)
        self.create_button(self.functions_frame, "Fit Line", self.handlers.fit_selected_line, 0, 1)
        self.create_button(self.functions_frame, "Find Single Lines", self.handlers.find_single_lines, 1, 0)
        self.create_button(self.functions_frame, "Single Slab Fit", self.handlers.single_slab_fit, 1, 1)

    def start(self):
        # Start the Tkinter main loop
        self.create_window()
        self.window.mainloop()