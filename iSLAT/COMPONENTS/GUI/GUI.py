import tkinter as tk
from tkinter import ttk
from .MainPlot import iSLATPlot
#from MainPlot import iSLATPlot

class GUI:
    def __init__(self, master, isotopologue_data, input_spectrum_data, data_field, wave_data, flux_data, mols, basem, isot, xp1, xp2, config):
        self.master = master
        self.isotopologue_data = isotopologue_data
        self.data_field = data_field
        self.mols = mols
        self.basem = basem
        self.isot = isot
        self.config = config
        self.nb_of_columns = 10  # Adjust as needed
        self.MainPlot = iSLATPlot(window=master,
                                  input_spectrum_data=input_spectrum_data,
                                  molecules_data=isotopologue_data,
                                  xp1=xp1,
                                  xp2=xp2,
                                  wave_data=wave_data,
                                  flux_data=flux_data,
                                  config=config)

    def create_button(self, frame, text, command, row, column, width=13, height=1):
        button_theme = self.config["user_settings"]["theme"]["buttons"].get(text.replace(" ", ""), self.config["user_settings"]["theme"]["buttons"]["DefaultBotton"])
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
        button.grid(row=row, column=column)
        return button

    def create_window(self):
        self.window = tk.Toplevel(self.master)
        self.window.title("iSLAT GUI")

        # Title frame
        self.title_frame = tk.Frame(self.window, bg="gray")
        self.title_frame.grid(row=0, column=0, columnspan=self.nb_of_columns, sticky='ew')

        # Outer frame
        self.outer_frame = tk.Frame(self.window)
        self.outer_frame.grid(row=self.title_frame.grid_info()['row'] + self.title_frame.grid_info()['rowspan'],
                              column=0, rowspan=10, columnspan=5, sticky="nsew")

        # Files frame
        self.files_frame = tk.Frame(self.window, borderwidth=2, relief="groove")
        self.files_frame.grid(row=self.outer_frame.grid_info()['row'] + self.outer_frame.grid_info()['rowspan'],
                              column=0, rowspan=7, columnspan=5, sticky="nsew")

        # Plot parameters frame
        self.plotparams_frame = tk.Frame(self.window, borderwidth=2, relief="groove")
        self.plotparams_frame.grid(row=self.files_frame.grid_info()['row'] + self.files_frame.grid_info()['rowspan'],
                                   column=0, rowspan=6, columnspan=5, sticky="nsew")

        # Functions frame
        self.functions_frame = tk.Frame(self.window, borderwidth=2, relief="groove")
        self.functions_frame.grid(row=self.plotparams_frame.grid_info()["row"] + self.plotparams_frame.grid_info()["rowspan"],
                      column=0, rowspan=6, columnspan=5, sticky="nsew")

        # Add buttons to the functions frame
        save_button = self.create_button(self.functions_frame, "Save Line", Save, 0, 0)
        fit_button = self.create_button(self.functions_frame, "Fit Line", fit_onselect, 0, 1)
        savedline_button = self.create_button(self.functions_frame, "Show Saved Lines", print_saved_lines, 1, 0)
        fitsavedline_button = self.create_button(self.functions_frame, "Fit Saved Lines", fit_saved_lines, 1, 1)
        autofind_button = self.create_button(self.functions_frame, "Find Single Lines", single_finder, 2, 0)
        atomlines_button = self.create_button(self.functions_frame, "Show Atomic Lines", print_atomic_lines, 2, 1)
        slabfit_button = self.create_button(self.functions_frame, "Single Slab Fit", run_slabfit, 3, 0)
        deblender_button = self.create_button(self.functions_frame, "Line De-blender", fitmulti_onselect, 3, 1)

        # Text frame
        self.text_frame = tk.Frame(self.window)
        self.text_frame.grid(row=self.functions_frame.grid_info()['row'] + self.functions_frame.grid_info()['rowspan'],
                             column=0, columnspan=5, sticky='nsew')

        # Data field
        self.data_field = tk.Text(self.text_frame, wrap="word", height=13, width=24)
        self.data_field.pack(fill="both", expand=True)

        # Add matplotlib canvas to the plot parameters frame
        plot_canvas = self.MainPlot.create_plot()
        plot_canvas.get_tk_widget().pack(padx=5, pady=5, fill="both", expand=True)

    def start(self):
        # Start the Tkinter main loop
        self.create_window()
        self.window.mainloop()