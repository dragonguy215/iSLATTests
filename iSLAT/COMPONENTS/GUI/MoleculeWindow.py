import tkinter as tk
from tkinter import colorchooser

class MoleculeWindow:
    def __init__(self, name, parent_frame, molecule_data, output_plot, config, islat):
        self.parent_frame = parent_frame
        self.name = name
        self.plot = output_plot
        self.config = config
        self.theme = config["theme"]
        self.islat = islat

        self.molecule_list = list(self.islat.molecules.keys())
        self.molecules = {}

        self.build_table()

    def build_table(self):
        frame = tk.LabelFrame(self.parent_frame, text="Molecules")
        frame.pack(fill="both", expand=True, padx=5, pady=5)

        headers = ["Molecule", "Temp", "Radius", "Density", "On", "Color"]
        for col, text in enumerate(headers):
            tk.Label(frame, text=text).grid(row=0, column=col)

        for i, mol in enumerate(self.molecule_list):
            init = self.islat.initial_values[mol]

            lbl = tk.Label(frame, text=mol)
            lbl.grid(row=i+1, column=0)

            temp_entry = tk.Entry(frame, width=6)
            temp_entry.insert(0, f"{init['t_kin']}")
            temp_entry.grid(row=i+1, column=1)

            rad_entry = tk.Entry(frame, width=6)
            rad_entry.insert(0, f"{init['radius_init']}")
            rad_entry.grid(row=i+1, column=2)

            dens_entry = tk.Entry(frame, width=6)
            dens_entry.insert(0, f"{init['n_mol_init']:.1e}")
            dens_entry.grid(row=i+1, column=3)

            on_var = tk.BooleanVar(value=True)
            on_btn = tk.Checkbutton(frame, variable=on_var, command=self.update_lines)
            on_btn.grid(row=i+1, column=4)

            color = self.theme["default_molecule_colors"][i]
            color_btn = tk.Button(frame, bg=color, width=4, command=lambda m=mol: self.pick_color(m))
            color_btn.grid(row=i+1, column=5)

            self.molecules[mol] = {
                "temp_entry": temp_entry,
                "rad_entry": rad_entry,
                "dens_entry": dens_entry,
                "on_var": on_var,
                "color": color
            }

        # initially draw
        self.update_lines()

    def pick_color(self, mol_name):
        color_code = colorchooser.askcolor(title=f"Pick color for {mol_name}")[1]
        if color_code:
            self.molecules[mol_name]["color"] = color_code
            self.update_lines()

    def update_lines(self):
        self.plot.clear_model_lines()
        for mol, props in self.molecules.items():
            if props["on_var"].get():
                temp = float(props["temp_entry"].get())
                rad = float(props["rad_entry"].get())
                dens = float(props["dens_entry"].get())
                color = props["color"]
                self.plot.add_model_line(mol, temp, rad, dens, color)
        self.plot.canvas.draw_idle()