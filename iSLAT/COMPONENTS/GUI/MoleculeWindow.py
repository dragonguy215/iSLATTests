import tkinter as tk
from tkinter import colorchooser
from functools import partial
'''from tkinter import filedialog
from .MainPlot import iSLATPlot
from .Data_field import DataField
from .Tooltips import CreateToolTip
from .GUIFunctions import GUIHandlers 
'''

class MoleculeWindow:
    def __init__(self, name, parent, molecule_data, config, islat_class_ref):
        self.name = name
        self.parent = parent
        self.molecule_data = molecule_data
        self.config = config
        self.theme = config["theme"]
        self.islat_class = islat_class_ref

        # Molecule table
        self.table_frame = tk.LabelFrame(parent, text="Molecules")
        self.table_frame.pack(fill="both", expand=True, padx=5, pady=5)
        headers = ['Molecule', 'Temp.', 'Radius', 'Col. Dens', 'On', 'Del.', 'Color', 'Visibility']
        for col, text in enumerate(headers):
            tk.Label(self.table_frame, text=text, bg=self.theme["background"], fg=self.theme["foreground"]).grid(row=0, column=col)

        # Fill rows with default values from initial_values and theme
        for i, molecule in enumerate(self.molecule_data):
            mol_name = molecule["name"]
            defaults = self.islat_class.initial_values.get(mol_name, {
                "t_kin": "", "radius_init": "", "scale_number": "", "scale_exponent": "", "n_mol_init": ""
            })
            default_color = self.theme["default_molecule_colors"][i] if i < len(self.theme["default_molecule_colors"]) else "#FFFFFF"

            tk.Label(self.table_frame, text=mol_name,
                     bg=self.theme["background"], fg=self.theme["foreground"]).grid(row=i+1, column=0)

            t_entry = tk.Entry(self.table_frame)
            t_entry.insert(0, str(defaults["t_kin"]))
            t_entry.grid(row=i+1, column=1)

            r_entry = tk.Entry(self.table_frame)
            r_entry.insert(0, str(defaults["radius_init"]))
            r_entry.grid(row=i+1, column=2)

            n_entry = tk.Entry(self.table_frame)
            n_entry.insert(0, str(defaults["n_mol_init"]))
            n_entry.grid(row=i+1, column=3)

            # On button to toggle line profiles
            def toggle_line_profiles(mol_name, var):
                if var.get():
                    self.islat_class.enable_line_profiles(mol_name)
                    print(f"Line profiles enabled for {mol_name}")
                else:
                    self.islat_class.disable_line_profiles(mol_name)
                    print(f"Line profiles disabled for {mol_name}")

            on_var = tk.BooleanVar(value=False)
            on_button = tk.Checkbutton(self.table_frame, variable=on_var, command=lambda mol=mol_name, var=on_var: toggle_line_profiles(mol, var))
            on_button.grid(row=i+1, column=4)

            tk.Checkbutton(self.table_frame).grid(row=i+1, column=5)

            def open_color_selector(button, mol_name):
                color_code = colorchooser.askcolor(title="Choose Color")[1]
                if color_code:
                    button.config(bg=color_code)
                    print(f"Selected color for {mol_name}: {color_code}")

            color_button = tk.Button(self.table_frame, text="", bg=default_color)
            color_button.config(command=partial(open_color_selector, color_button, mol_name))
            color_button.grid(row=i+1, column=6)

            # Visibility Checkbutton
            visibility_var = tk.BooleanVar(value=(mol_name.lower() == 'h2o'))
            visibility_button = tk.Checkbutton(
                self.table_frame, text='', variable=visibility_var,
                command=lambda mn=mol_name.lower(): self.model_visible(mn, visibility_var.get())
            )
            visibility_button.grid(row=i+1, column=7)

    def model_visible(self, mol_name, is_visible):
        if is_visible:
            print(f"{mol_name} is now visible")
        else:
            print(f"{mol_name} is now hidden")