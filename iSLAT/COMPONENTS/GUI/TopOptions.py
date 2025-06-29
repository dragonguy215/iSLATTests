import tkinter as tk
from tkinter import ttk
from .GUIFunctions import create_button

class TopOptions:
    def __init__(self, master, islat, theme):
        self.master = master
        self.islat = islat
        self.theme = theme
        #self.theme = self.master.theme

        # Create the frame for top options
        self.frame = tk.Frame(master, borderwidth=2, relief="groove")
        self.frame.pack(side="top", fill="x")

        # Create buttons for top options
        create_button(self.frame, self.theme, "HITRAN Query", self.hitran_query, 0, 0)
        create_button(self.frame, self.theme, "Spectra Browser", self.spectra_browser, 0, 1)
        create_button(self.frame, self.theme, "Default Molecules", self.default_molecules, 1, 0)
        create_button(self.frame, self.theme, "Add Molecule", self.add_molecule, 1, 1)
        create_button(self.frame, self.theme, "Save Parameters", self.save_parameters, 0, 2)
        create_button(self.frame, self.theme, "Load Parameters", self.load_parameters, 1, 2)
        create_button(self.frame, self.theme, "Export Models", self.export_models, 0, 3)
        create_button(self.frame, self.theme, "Toggle Legend", self.toggle_legend, 1, 3)
    
    def hitran_query(self):
        print("Perform HITRAN query")
    
    def spectra_browser(self):
        print("Open spectra browser")

    def default_molecules(self):
        print("Default molecules loaded")

    def add_molecule(self):
        

    def save_parameters(self):
        print("Save parameters to file")

    def load_parameters(self):
        print("Load parameters from file")

    def export_models(self):
        print("Export models to file")

    def toggle_legend(self):
        #print("Toggled legend on plot")
        self.islat.GUI.plot.toggle_legend()