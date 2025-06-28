iSLAT_version = 'v5.00.00'
print(' ')
print('Loading iSLAT ' + iSLAT_version + ': Please Wait ...')

# Import necessary modules
import numpy as np
import pandas as pd
import warnings
import matplotlib
import os
import json

# matplotlib.use('Agg')
matplotlib.use("TKAgg")
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

context = ssl.create_default_context(cafile=certifi.where ())

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
        self.update_initial_molecule_parameters()

        '''self.root = tk.Tk()
        self.root.title("iSLAT - Infrared Spectral Line Analysis Tool")
        self.root.geometry("800x600")
        self.root.resizable(True, True)'''
        
        self.mols = ["H2", "HD", "H2O", "H218O", "CO2", "13CO2", "CO", "13CO", "C18O", "CH4", "HCN", "H13CN", "NH3", "OH", "C2H2", "13CCH2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.basem = ["H2", "H2", "H2O", "H2O", "CO2", "CO2", "CO", "CO", "CO", "CH4", "HCN", "HCN", "NH3", "OH", "C2H2", "C2H2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.isot = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1]

        self.wave_range = (0.3, 1000)

        self.min_vu = 1 / (self.wave_range[0] / 1E6) / 100.
        self.max_vu = 1 / (self.wave_range[1] / 1E6) / 100.

        # Check for HITRAN files and download if necessary
        self.check_HITRAN()

        self.molecules_data_default = molecules_data.copy()
        self.deleted_molecules = []

        self.xp1 = self.xp2 = None

    def init_gui(self):
        """
        Initialize the GUI components of iSLAT.
        This function sets up the main window, menus, and other GUI elements.
        """
        if not hasattr(self, "root"):
            self.root = tk.Tk()
            self.root.title("iSLAT - Infrared Spectral Line Analysis Tool")
            self.root.resizable(True, True)

        #self.root.mainloop()

        if not hasattr(self, "GUI"):
            self.GUI = GUI(
                master=self.root,
                molecule_data=self.molecules_data_default,
                wave_data=self.wave_data,
                flux_data=self.flux_data,
                config=self.user_settings,
                islat_class_ref=self
            )
        
        self.GUI.start()

    def init_molecules(self):
        self.molecules = {}
        self.initial_values = {}

        for mol_entry in molecules_data:
            mol_name = mol_entry["name"]
            mol_filepath = mol_entry["file"]
            mol_label = mol_entry["label"]
            mol_name_lower = mol_name.lower()

            # Load the molecule line data
            mol_data = MolData(mol_name, mol_filepath)

            # Get precomputed parameters
            params = self.initial_molecule_parameters[mol_name]
            scale_exponent = params["scale_exponent"]
            scale_number = params["scale_number"]
            t_kin = params["t_kin"]
            radius_init = params["radius_init"]
            n_mol_init = float(scale_number * (10 ** scale_exponent))

            # Store in structured dicts
            self.molecules[mol_name] = mol_data
            self.initial_values[mol_name] = {
                "scale_exponent": scale_exponent,
                "scale_number": scale_number,
                "t_kin": t_kin,
                "radius_init": radius_init,
                "n_mol_init": n_mol_init
            }

            print(f"Molecule Initialized: {mol_name}")

    def run(self):
        """
        Run the iSLAT application.
        This function starts the main event loop of the Tkinter application.
        """
        # Start the main event loop
        self.load_spectrum()
        self.init_molecules()
        self.err_data = np.full_like(self.flux_data, np.nanmedian(self.flux_data)/100)
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

    def load_spectrum(self, file_path=None):
        filetypes = [('CSV Files', '*.csv')]
        spectra_directory = os.path.abspath("EXAMPLE-data")
        #file_path = filedialog.askopenfilename(title='Choose Spectrum Data File', filetypes=filetypes, initialdir=spectra_directory)
        if file_path is None:
            file_path = filedialog.askopenfilename(title='Choose Spectrum Data File', filetypes=filetypes, initialdir=spectra_directory)

        if file_path:
            df = pd.read_csv(file_path)
            self.wave_data = df['wave'].values
            self.flux_data = df['flux'].values
            #print(np.min(self.wave_data), np.max(self.wave_data))
            #print(type(self.wave_data), type(np.min(self.wave_data)))
            lam_min = np.min(self.wave_data)
            lam_max = np.max(self.wave_data)
            self.input_spectrum = Spectrum(
                lam_min = lam_min,
                lam_max = lam_max,
                dlambda= (lam_max - lam_min) / len(self.wave_data),
                R= 1 / (lam_min / 1E6) / 100.,  # Convert to cm-1
            )
            print(f"Loaded spectrum from {file_path}")
        else:
            print("No file selected.")
    
    def generate_population_diagram(self):
        if not hasattr(self, "selected_lines") or not self.selected_lines:
            print("No lines selected yet to build population diagram.")
            return None, None

        energies, pops = [], []
        for wave in self.selected_lines:
            for mol in self.molecules.values():
                closest_line = mol.find_closest_line(wave)
                if closest_line:
                    e_upper = closest_line['E_up']
                    g_u = closest_line['g_u']
                    n_u = Intensity.calculate_population(e_upper, g_u)
                    energies.append(e_upper)
                    pops.append(np.log(n_u / g_u))

        return np.array(energies), np.array(pops)
    
    def run_single_slab_fit(self):
        loader = DataLoader(self.molecule_data)
        self.slab_model = ModelFitting(loader)
        result = self.slab_model.fit()
        return result.summary()

    def find_single_lines(self):
        threshold = line_threshold * np.max(self.flux_data)
        sep = 0.1
        found = []
        for i, val in enumerate(self.flux_data):
            if val > threshold:
                if len(found) == 0 or (self.wave_data[i] - found[-1]) > sep:
                    found.append(self.wave_data[i])
        self.selected_lines = found
        print(f"Found {len(found)} lines.")
        return found

    def save_line(self, line_info):
        df = pd.DataFrame([line_info])
        file = os.path.join("SAVES", "lines_saved.csv")
        df.to_csv(file, mode='a', header=not os.path.exists(file), index=False)
        print(f"Line saved to {file}")

    '''def fit_single_gaussian(self, x, y, err):
        model = GaussianModel()
        params = model.guess(y, x=x)
        fit_result = model.fit(y, params, x=x, weights=1/err, nan_policy='omit')
        print(fit_result.fit_report())
        return fit_result

    def fit_multi_gaussian(self, x, y, err):
        # try to fit 2 components as an example of deblend
        g1 = GaussianModel(prefix='g1_')
        g2 = GaussianModel(prefix='g2_')
        model = g1 + g2
        params = g1.guess(y, x=x) + g2.guess(y, x=x)
        fit_result = model.fit(y, params, x=x, weights=1/err, nan_policy='omit')
        print(fit_result.fit_report())
        return fit_result'''

    def fit_selected_line(self, xmin, xmax, deblend=False):
        x_fit = self.wave_data[(self.wave_data >= xmin) & (self.wave_data <= xmax)]
        y_fit = self.flux_data[(self.wave_data >= xmin) & (self.wave_data <= xmax)]
        err = self.err_data[(self.wave_data >= xmin) & (self.wave_data <= xmax)]

        if len(x_fit) < 5:
            print("Not enough data points to fit.")
            return None

        print(f"Fitting line in range: {xmin:.4f}-{xmax:.4f}, points: {len(x_fit)}")

        if deblend:
            g1 = GaussianModel(prefix='g1_')
            g2 = GaussianModel(prefix='g2_')
            model = g1 + g2
            params = g1.guess(y_fit, x=x_fit) + g2.guess(y_fit, x=x_fit)
        else:
            model = GaussianModel()
            params = model.guess(y_fit, x=x_fit)

        fit_result = model.fit(y_fit, params, x=x_fit, weights=1/err, nan_policy='omit')
        print(fit_result.fit_report())

        return fit_result