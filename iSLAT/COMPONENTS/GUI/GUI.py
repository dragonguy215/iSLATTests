#from .GUI import *

class GUI:
    def __init__(self, master, isotopologue_data, data_field, mols, basem, isot):
        self.master = master
        self.isotopologue_data = isotopologue_data
        self.data_field = data_field
        self.mols = mols
        self.basem = basem
        self.isot = isot

        # Create the main window for downloading HITRAN data
        self.create_window()

    def create_window(self):
        import tkinter as tk
        from tkinter import ttk
        #from iSLAT.COMPONENTS.HITRAN import download_hitran_data

        # Initialize the window with the master widget
        self.window = None