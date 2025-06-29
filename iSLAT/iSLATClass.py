iSLAT_version = 'v5.00.00'
print(' ')
print('Loading iSLAT ' + iSLAT_version + ': Please Wait ...')

# Import necessary modules
import numpy as np
import pandas as pd
import os
import json

from lmfit.models import GaussianModel
import tkinter as tk
from tkinter import filedialog
import ssl
import certifi
import datetime

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
from .COMPONENTS.Molecule import Molecule
from .COMPONENTS.MoleculeDict import MoleculeDict

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
        
        self.mols = ["H2", "HD", "H2O", "H218O", "CO2", "13CO2", "CO", "13CO", "C18O", "CH4", "HCN", "H13CN", "NH3", "OH", "C2H2", "13CCH2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.basem = ["H2", "H2", "H2O", "H2O", "CO2", "CO2", "CO", "CO", "CO", "CH4", "HCN", "HCN", "NH3", "OH", "C2H2", "C2H2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.isot = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1]

        self.wave_range = wavelength_range

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
                molecule_data=self.molecules_dict,
                wave_data=self.wave_data,
                flux_data=self.flux_data,
                config=self.user_settings,
                islat_class_ref=self
            )
        
        self.GUI.start()

    def init_molecules(self):
        self.molecules_dict = MoleculeDict()
        #self.initial_values = {}

        self.molecules_dict.load_molecules_data(molecules_data=self.molecules_data_default,
                                                initial_molecule_parameters=self.initial_molecule_parameters,
                                                save_file_data = self.savedata,
                                                wavelength_range = self.wave_range, 
                                           intrinsic_line_width = intrinsic_line_width, 
                                           model_pixel_res = model_pixel_res, 
                                           model_line_width = model_line_width,
                                           dist = dist)

    def run(self):
        """
        Run the iSLAT application.
        This function starts the main event loop of the Tkinter application.
        """
        # Start the main event loop
        self.get_save_data()
        self.load_spectrum()
        self.init_molecules()
        #self.err_data = np.full_like(self.flux_data, np.nanmedian(self.flux_data)/100)
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

    def get_save_data(self):
        """
        get_save_data() loads the save data from the SAVES folder.
        It returns a list of dictionaries with the save data.
        """
        save_file = os.path.join("SAVES", "molecules_list.csv")
        if os.path.exists(save_file):
            try:
                df = pd.read_csv(save_file)
                self.savedata = {row['Molecule Name']: {col: row[col] for col in df.columns if col != 'Molecule Name'} for _, row in df.iterrows()}
                #print("self.savedata:", self.savedata)
            except Exception as e:
                print(f"Error reading save file: {e}")
                self.savedata = []
        else:
            print("No save file found.")
            self.savedata = []

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
        os.makedirs("SAVES", exist_ok=True)
        os.makedirs("MODELS", exist_ok=True)
        os.makedirs("LINESAVES", exist_ok=True)
        os.makedirs("HITRANdata", exist_ok=True)

    def load_spectrum(self, file_path=None):
        filetypes = [('CSV Files', '*.csv')]
        spectra_directory = os.path.abspath("EXAMPLE-data")
        if file_path is None:
            file_path = filedialog.askopenfilename(title='Choose Spectrum Data File', filetypes=filetypes, initialdir=spectra_directory)

        if file_path:
            df = pd.read_csv(file_path)
            self.wave_data = np.array(df['wave'].values)
            self.wave_data_original = self.wave_data.copy()
            self.flux_data = np.array(df['flux'].values)
            self.err_data = np.array(df['err'].values)
            self.continuum_data = np.array(df['cont'].values)
            print(f"Loaded spectrum from {file_path}")
        else:
            print("No file selected.")

    def update_model_spectrum(self):
        summed_flux = np.zeros_like(self.wave_data)
        for mol in self.molecules_dict.values():
            if mol.is_active:
                mol_flux = mol.get_flux(self.wave_data)
                summed_flux += mol_flux
        self.sum_spectrum_flux = summed_flux

    def run_single_slab_fit(self):
        loader = DataLoader(self.molecules_data_default)
        self.slab_model = ModelFitting(loader)
        result = self.slab_model.fit()
        return result.summary()

    def find_single_lines(self):
        """ find_single_lines() identifies spectral lines in the flux data based on a threshold.
        It returns a list of wavelengths where the flux exceeds the threshold.
        The threshold is calculated as a fraction of the maximum flux value."""
        threshold = self.user_settings.get("line_threshold", 0.03) * np.max(self.flux_data)
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