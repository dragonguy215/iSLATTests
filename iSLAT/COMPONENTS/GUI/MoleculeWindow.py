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

        self.molecules_dict = self.islat.molecules_dict
        self.molecules = {}

        self.build_table()

    def build_table(self):
        self.frame = tk.LabelFrame(self.parent_frame, text="Molecules")
        self.frame.pack(fill="both", expand=True, padx=5, pady=5)

        headers = ["Molecule", "Temp", "Radius", "Density", "On", "Color"]
        for col, text in enumerate(headers):
            tk.Label(self.frame, text=text).grid(row=0, column=col)

        for i, mol_name in enumerate(self.molecules_dict.keys()):
            mol_data = self.molecules_dict[mol_name]

            lbl = tk.Label(self.frame, text=mol_name)
            lbl.grid(row=i+1, column=0)

            temp_entry = tk.Entry(self.frame, width=6)
            temp_entry.insert(0, f"{mol_data.temp}")
            temp_entry.grid(row=i+1, column=1)

            rad_entry = tk.Entry(self.frame, width=6)
            rad_entry.insert(0, f"{mol_data.radius}")
            rad_entry.grid(row=i+1, column=2)

            dens_entry = tk.Entry(self.frame, width=6)
            dens_entry.insert(0, f"{mol_data.n_mol_init:.1e}")
            dens_entry.grid(row=i+1, column=3)

            on_var = tk.BooleanVar(value=mol_data.is_active)
            on_btn = tk.Checkbutton(self.frame, variable=on_var, command=self.update_lines)
            on_btn.grid(row=i+1, column=4)

            color = self.theme["default_molecule_colors"][i]
            color_btn = tk.Button(self.frame, bg=color, width=4, command=lambda m=mol_name: self.pick_color(m))
            color_btn.grid(row=i+1, column=5)

            self.molecules[mol_name] = {
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
                # actually update molecule parameters
                m_obj = self.islat.molecules_dict[mol]
                m_obj.temp = temp
                m_obj.radius = rad
                m_obj.n_mol_init = dens
                m_obj.color = color
                self.plot.add_model_line(mol, temp, rad, dens, color)
        self.islat.update_model_spectrum()
        self.plot.update_line_inspection_plot()