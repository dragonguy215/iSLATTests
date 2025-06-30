import tkinter as tk
from tkinter import ttk

class ControlPanel:
    def __init__(self, master, islat):
        self.master = master
        self.islat = islat

        # Create the control panel frame
        self.frame = tk.Frame(master, borderwidth=2, relief="groove")
        self.frame.pack(side="left", fill="y")

        # Create and place entry fields
        #self.create_entry("Plot start:", 0, 0, "ax1_starting_x", self.update_xp1_rng)
        #self.create_entry("Plot range:", 0, 2, "ax1_range_x", self.update_xp1_rng)
        self.create_plot_start(0, 0)
        self.create_plot_range(0, 2)
        self.create_wavelength_range(1, 0, 1, 2)
        
        #self.create_entry("Min. Wave:", 1, 0, "min_wavelength", self.update_initvals)
        #self.create_entry("Max. Wave:", 1, 2, "max_wavelength", self.update_initvals)
        self.create_entry("Distance:", 2, 0, "distance", self.update_initvals)
        self.create_entry("Stellar RV:", 2, 2, "star_rv", self.update_initvals)
        self.create_entry("FWHM:", 3, 0, "fwhm", self.update_initvals)
        self.create_entry("Broadening:", 3, 2, "intrinsic_line_width", self.update_initvals)

        self.create_molecule_dropdown(4, 0)
        self.reload_molecule_dropdown()

    def create_molecule_dropdown(self, row, column):
        label = tk.Label(self.frame, text="Molecule:")
        label.grid(row=row, column=column, padx=5, pady=5)

        dropdown_options = list(self.islat.molecules_dict.keys()) + ["SUM", "ALL"]
        self.molecule_var = tk.StringVar(self.frame)
        self.molecule_var.set(dropdown_options[0])  # Default to the first option

        self.dropdown = ttk.Combobox(self.frame, textvariable=self.molecule_var, values=dropdown_options)
        self.dropdown.grid(row=row, column=column + 1, padx=5, pady=5)

        # Update self.islat.active_molecule when a new molecule is selected
        self.dropdown.bind("<<ComboboxSelected>>", lambda event: setattr(self.islat, 'active_molecule', self.molecule_var.get()))

    def reload_molecule_dropdown(self):
        dropdown_options = list(self.islat.molecules_dict.keys()) + ["SUM", "ALL"]
        self.dropdown['values'] = dropdown_options
        if self.molecule_var.get() not in dropdown_options:
            self.molecule_var.set(dropdown_options[0])

    def create_entry(self, label_text, row, column, attribute_name, callback):
        label = tk.Label(self.frame, text=label_text)
        label.grid(row=row, column=column, padx=5, pady=5)
        
        entry = tk.Entry(self.frame, bg='lightgray', width=8)
        entry.grid(row=row, column=column + 1, padx=5, pady=5)
        entry.bind("<Return>", lambda event: callback())
        
        # Tie the entry field to a property
        setattr(self, f"_{attribute_name}_entry", entry)
        setattr(self, attribute_name, property(
            lambda self: self._get_entry_value(attribute_name),
            lambda self, value: self._set_entry_value(attribute_name, value, callback)
        ))

    def _get_entry_value(self, attribute_name):
        entry = getattr(self, f"_{attribute_name}_entry")
        return entry.get()

    def _set_entry_value(self, attribute_name, value, callback):
        entry = getattr(self, f"_{attribute_name}_entry")
        entry.delete(0, tk.END)
        entry.insert(0, str(value))
        callback()

    def create_plot_start(self, row, column):
        label = tk.Label(self.frame, text="Plot start:")
        label.grid(row=row, column=column, padx=5, pady=5)

        self.ax1_starting_x = tk.Entry(self.frame, bg='lightgray', width=8)
        self.ax1_starting_x.grid(row=0, column=1, padx=5, pady=5)
        self.ax1_starting_x.insert(0, str(self.islat.display_range[0]))  # Load initial value from islat.display_range
        self.ax1_starting_x.bind("<Return>", lambda _: self.update_xp1_rng())

    def create_plot_range(self, row, column):
        label = tk.Label(self.frame, text="Plot range:")
        label.grid(row=row, column=column, padx=5, pady=5)

        self.ax1_range_x = tk.Entry(self.frame, bg='lightgray', width=8)
        self.ax1_range_x.grid(row=row, column=3, padx=5, pady=5)
        self.ax1_range_x.insert(0, str(self.islat.display_range[1] - self.islat.display_range[0]))  # Load initial value from islat.display_range
        self.ax1_range_x.bind("<Return>", lambda _: self.update_xp1_rng())

    def update_xp1_rng(self):
        try:
            xp1_value = float(self.ax1_starting_x.get())
            rng_value = float(self.ax1_range_x.get())
            self.islat.display_range = (xp1_value, xp1_value + rng_value)  # Update the display range with both values
            print(f"Updated display_range to start: {xp1_value}, range: {rng_value}")
        except ValueError:
            print("Invalid input for xp1 or rng")

    def create_wavelength_range(self, minrow, mincolumn, maxrow, maxcolumn):
        label = tk.Label(self.frame, text="Min. Wave:")
        label.grid(row=minrow, column=mincolumn, padx=5, pady=5)
        self.min_wavelength = tk.Entry(self.frame, bg='lightgray', width=8)
        self.min_wavelength.grid(row=minrow, column=mincolumn + 1, padx=5, pady=5)
        self.min_wavelength.insert(0, str(self.islat.wavelength_range[0]))
        self.min_wavelength.bind("<Return>", lambda _: self.update_wavelength_range())

        label = tk.Label(self.frame, text="Max. Wave:")
        label.grid(row=maxrow, column=maxcolumn, padx=5, pady=5)
        self.max_wavelength = tk.Entry(self.frame, bg='lightgray', width=8)
        self.max_wavelength.grid(row=maxrow, column=maxcolumn + 1, padx=5, pady=5)
        self.max_wavelength.insert(0, str(self.islat.wavelength_range[1]))
        self.max_wavelength.bind("<Return>", lambda _: self.update_wavelength_range())

    def update_wavelength_range(self):
        try:
            min_wave = float(self.min_wavelength.get())
            max_wave = float(self.max_wavelength.get())
            if min_wave < max_wave:
                self.islat.wavelength_range = (min_wave, max_wave)
                print(f"Updated wavelength range to: {self.islat.wavelength_range}")
            else:
                print("Min wavelength must be less than max wavelength")
        except ValueError:
            print("Invalid input for wavelength range")

    def update_initvals(self):
        # Placeholder for initialization values update logic
        print("Initialization values updated")