iSLAT_version = 'v4.06.02'
print (' ')
print ('Loading iSLAT ' + iSLAT_version + ': Please Wait ...')

# Import necessary modules
import numpy as np
import pandas as pd
import warnings
import matplotlib
import os
import json

# matplotlib.use('Agg')
matplotlib.use ("TKAgg")
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, SpanSelector, TextBox, CheckButtons
from matplotlib.artist import Artist
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import sys
from astropy.io import ascii, fits
from astropy.table import vstack, Table
from astropy import stats
import lmfit
from lmfit.models import GaussianModel
from lmfit.models import PseudoVoigtModel
import tkinter as tk
from tkinter import filedialog, simpledialog, ttk, Toplevel, Label, LEFT, SOLID  # For ttk.Style
from tkinter import colorchooser
import inspect
# from PyQt5.QtWidgets import QApplication, QMainWindow
from datetime import datetime as dt
import time
import threading
from astroquery import hitran
from astropy import units as un
from scipy import constants as con
import datetime
import certifi
import ssl
import urllib
import webbrowser

context = ssl.create_default_context (cafile=certifi.where ())

from .ir_model import *
from .COMPONENTS.chart_window import MoleculeSelector
from .COMPONENTS.Hitran_data import get_Hitran_data
from .COMPONENTS.partition_function_writer import write_partition_function
from .COMPONENTS.line_data_writer import write_line_data
from .COMPONENTS.slabfit_config import *
from .COMPONENTS.slabfit_loader import *
from .COMPONENTS.slabfit_runner import *
from .iSLATDefaultInputParms import *
from .iSLATCSVHandling import *
from .COMPONENTS.GUI import *

class iSLAT:
    """
    iSLAT class to handle the iSLAT functionalities.
    This class is used to initialize the iSLAT application, load user settings, and manage the main functionalities.
    """

    def __init__(self):
        """
        Initialize the iSLAT application.
        """
        self.create_folders()

        # Load settings
        self.user_settings = self.load_user_settings()
        self.update_default_molecule_parameters()

        '''self.root = tk.Tk()
        self.root.title("iSLAT - Infrared Spectral Line Analysis Tool")
        self.root.geometry("800x600")
        self.root.resizable(True, True)'''
        
        self.mols = ["H2", "HD", "H2O", "H218O", "CO2", "13CO2", "CO", "13CO", "C18O", "CH4", "HCN", "H13CN", "NH3", "OH", "C2H2", "13CCH2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.basem = ["H2", "H2", "H2O", "H2O", "CO2", "CO2", "CO", "CO", "CO", "CH4", "HCN", "HCN", "NH3", "OH", "C2H2", "C2H2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.isot = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1]

        """self.min_wave = 0.3  # micron
        self.max_wave = 1000  # micron"""
        self.wave_range = (0.3, 1000)

        self.min_vu = 1 / (self.wave_range[0] / 1E6) / 100.
        self.max_vu = 1 / (self.wave_range[1] / 1E6) / 100.

        # Check for HITRAN files and download if necessary
        self.check_HITRAN()

        self.molecules_data_default = molecules_data.copy()
        self.deleted_molecules = []

        self.xp1 = None
        self.xp2 = None  

    def init_gui(self):
        """
        Initialize the GUI components of iSLAT.
        This function sets up the main window, menus, and other GUI elements.
        """
        if not hasattr(self, "root"):
            self.root = tk.Tk()
            self.root.title("iSLAT - Infrared Spectral Line Analysis Tool")
            self.root.geometry("800x600")
            self.root.resizable(True, True)

        #self.root.mainloop()

        if not hasattr(self, "GUI"):
            self.GUI = GUI(
                master=self.root,
                isotopologue_data=self.molecules_data_default,
                input_spectrum_data = self.input_spectrum_data if hasattr(self, 'input_spectrum_data') else None,
                data_field=None,  # Placeholder for data field, to be set later
                wave_data = self.wave_data if hasattr(self, 'wave_data') else None,
                flux_data = self.flux_data if hasattr(self, 'flux_data') else None,
                mols=self.mols,
                basem=self.basem,
                isot=self.isot,
                xp1=self.xp1,
                xp2=self.xp2,
                config={"iSLAT_version": "2.0", "user_settings":{'first_startup': False, 'reload_default_files': False, 'theme': {'theme_name': 'Light Theme', 'description': 'A light theme for iSLAT with a white background and black text.', 'foreground': '#000000', 'background': '#FFFFFF', 'toolbar': '#000000', 'graph_fill_color': '#D3D3D3', 'selection_color': '#00FF00', 'uncertainty_band_color': '#ABABAB', 'toolbar_highlight_background': '#000000', 'toolbar_highlight_color': '#000000', 'toolbar_active_color': '#FFFFFF', 'buttons': {'DefaultBotton': {'background': '#D3D3D3', 'active_background': '#808080', 'description': 'If a button is not specified, this will be used as the default button style.'}, 'ToggleLegend': {'background': '#D3D3D3', 'active_background': '#808080'}}}}},

            )
        
        self.GUI.start()

    def run(self):
        """
        Run the iSLAT application.
        This function starts the main event loop of the Tkinter application.
        """
        # Start the main event loop
        self.selectfileinit()
    
        self.init_gui()

    def load_user_settings(self):
        """ load_user_settings() loads the user settings from the UserSettings.json file."""
        user_settings_file = "CONFIG/UserSettings.json"
        if os.path.exists(user_settings_file):
            with open(user_settings_file, 'r') as f:
                user_settings = json.load(f)
        else:
            # If the file does not exist, return default settings and save them as a new json file
            default_settings = {
                "first_startup": True,
                "reload_default_files": True,
                "theme": "LightTheme"
            }
            with open(user_settings_file, 'w') as f:
                json.dump(default_settings, f, indent=4)
            user_settings = default_settings
        
        # append theme information to the user settings dictonary
        theme_file = f"CONFIG/GUIThemes/{user_settings['theme']}.json"
        if os.path.exists(theme_file):
            with open(theme_file, 'r') as f:
                theme_settings = json.load(f)
            user_settings["theme"] = theme_settings
        return user_settings

    def update_default_molecule_parameters(self):
        """
        update_default_molecule_parameters() updates the default molecule parameters from the DefaultMoleculeParameters.json file.
        """
        #global default_molecule_parameters
        with open("CONFIG/DefaultMoleculeParameters.json", 'r') as f:
            default_molecule_parameters = json.load(f)["default_initial_params"]
        self.default_initial_parameters = default_molecule_parameters
    
    def update_initial_molecule_parameters(self):
        """
        update_initial_molecule_parameters() updates the initial molecule parameters from the DefaultMoleculeParameters.json file.
        """
        #global initial_molecule_parameters
        with open("CONFIG/DefaultMoleculeParameters.json", 'r') as f:
            initial_molecule_parameters = json.load(f)["initial_parameters"]
        self.initial_molecule_parameters = initial_molecule_parameters

    def check_HITRAN(self, print_statments = True):
        """ check_HITRAN(print_statments=True) checks if the HITRAN files are present and downloads them if necessary.
        If print_statments is True, it will print the status of the HITRAN files to the console."""
        if print_statments: print('\nChecking for HITRAN files: ...')
        # If this is the first startup or reload_default_files is True, download the default HITRAN files
        if self.user_settings["first_startup"] or self.user_settings["reload_default_files"]:
            print('First startup or reload_default_files is True. Downloading default HITRAN files ...')
            for mol, bm, iso in zip(self.mols, self.basem, self.isot):
                save_folder = 'HITRANdata'
                file_path = os.path.join(save_folder, "data_Hitran_2020_{:}.par".format(mol))

                if os.path.exists(file_path):
                    print("File already exists for mol: {:}. Skipping.".format(mol))
                    continue

                print("Downloading data for mol: {:}".format(mol))
                Htbl, qdata, M, G = get_Hitran_data(bm, iso, self.min_vu, self.max_vu)

                with open(file_path, 'w') as fh:
                    fh.write("# HITRAN 2020 {:}; id:{:}; iso:{:};gid:{:}\n".format(mol, M, iso, G))
                    fh.write("# Downloaded from the Hitran website\n")
                    fh.write("# {:s}\n".format(str(datetime.date.today())))
                    fh = write_partition_function(fh, qdata)
                    fh = write_line_data(fh, Htbl)

                print("Data for Mol: {:} downloaded and saved.".format(mol))

            self.user_settings["first_startup"] = False
            self.user_settings["reload_default_files"] = False
            with open("UserSettings.json", 'w') as f:
                json.dump(self.user_settings, f, indent=4)
        else:
            print('Not the first startup and reload_default_files is False. Skipping HITRAN files download.')

    def create_folders(self): # see if we need this one and/or add config for directories
        """
        create_folders() creates the necessary folders for saving data and models.
        This is typically done at the first launch of iSLAT.
        """
        # Create necessary folders, if they don't exist
        save_folder = "SAVES"
        os.makedirs(save_folder, exist_ok=True)
        output_dir = "MODELS"
        os.makedirs(output_dir, exist_ok=True)
        linesave_folder = "LINESAVES"
        os.makedirs(linesave_folder, exist_ok=True)
        # create HITRAN folder, only needed for first start
        HITRAN_folder = "HITRANdata"
        os.makedirs(HITRAN_folder, exist_ok=True)

    def selectfileinit(self):
        filetypes = [('CSV Files', '*.csv')]
        spectra_directory = os.path.abspath ("EXAMPLE-data")
        # Ask the user to select a file
        infiles = filedialog.askopenfilename(multiple=True, title='Choose Spectrum Data File', filetypes=filetypes, initialdir=spectra_directory)

        if infiles:
            for file_path in infiles:
                # Process each selected file
                print("Selected file:", file_path)
                file_name = os.path.basename(file_path)

                # code to process each file
                self.input_spectrum_data = pd.read_csv(filepath_or_buffer=file_path, sep=',')
                self.wave_data = np.array(self.input_spectrum_data['wave'])
                self.wave_original = np.array(self.input_spectrum_data['wave'])
                self.flux_data = np.array(self.input_spectrum_data['flux'])
                if 'err' in self.input_spectrum_data:
                    err_data = np.array(self.input_spectrum_data['err'])
                else:
                    err_data = np.full_like(self.flux_data, np.nanmedian(self.flux_data)/100)  # assumed, if not present

                    # Set initial values of xp1 and rng
                fig_max_limit = np.nanmax(self.wave_data)
                fig_min_limit = np.nanmin(self.wave_data)
                self.xp1 = np.around(fig_min_limit + (fig_max_limit - fig_min_limit) / 2, decimals=2)
                rng = np.around((fig_max_limit - fig_min_limit) / 10, decimals=2)
                self.xp2 = self.xp1 + rng
        else:
            print("No files selected.")
    
    def selectfile(self):
        filetypes = [('CSV Files', '*.csv')]
        spectra_directory = os.path.abspath("EXAMPLE-data")
        infiles = filedialog.askopenfilename(multiple=True, title='Choose Spectrum Data File', filetypes=filetypes,
                                            initialdir=spectra_directory)

        if infiles:
            for file_path in infiles:
                # Process each selected file
                print ("Selected file:", file_path)
                file_name = os.path.basename(file_path)

                file_name_label.config(text=str (file_name))
                # filename_box_data.set_val(file_name)
                # Add your code to process each file
                # THIS IS THE OLD FILE SYSTEM (THIS WILL BE USED UNTIL THE NEW FILE SYSTEM IS DEVELOPED) USE THIS!!!!!
                input_spectrum_data = pd.read_csv(filepath_or_buffer=(file_path), sep=',')
                wave_data = np.array(input_spectrum_data['wave'])
                wave_original = np.array(input_spectrum_data['wave'])
                flux_data = np.array(input_spectrum_data['flux'])
                if 'err' in input_spectrum_data:
                    err_data = np.array(input_spectrum_data['err'])
                else:
                    err_data = np.full_like(flux_data, np.nanmedian(flux_data)/100)  # assumed, if not present

                # Set new values of xp1 and rng only if the new spectrum is in a different wave range
                fig_max_limit = np.nanmax(wave_data)
                fig_min_limit = np.nanmin(wave_data)
                xp1_current = float(xp1_entry.get ())
                if xp1_current > fig_max_limit or xp1_current < fig_min_limit:
                    xp1 = fig_min_limit + (fig_max_limit - fig_min_limit) / 2
                    rng = (fig_max_limit - fig_min_limit) / 10
                    xp2 = xp1 + rng
                    xp1_entry.delete(0, "end")
                    xp1_entry.insert(0, np.around (xp1, decimals=2))
                    rng_entry.delete(0, "end")
                    rng_entry.insert(0, np.around (rng, decimals=2))

                # now = dt.now()
                # dateandtime = now.strftime("%d-%m-%Y-%H-%M-%S")
                # print(dateandtime)
                # svd_line_file = f'savedlines-{dateandtime}.csv'

                self.update()

                data_field.delete('1.0', "end")
                data_field.insert('1.0', 'New spectrum loaded!')
        else:
            data_field.delete('1.0', "end")
            data_field.insert('1.0', 'No file selected.')

    def update(self):
        """
        update() updates the GUI components and data fields.
        This function is called after loading a new spectrum or making changes to the GUI.
        """
        # Update the GUI components
        if hasattr(self, 'GUI'):
            self.GUI.update_gui()
        
        # Update the data field with a message
        if hasattr(self, 'data_field'):
            self.data_field.delete('1.0', "end")
            self.data_field.insert('1.0', 'Spectrum data updated.')