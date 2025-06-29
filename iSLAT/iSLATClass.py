iSLAT_version = 'v5.00.00'
print(' ')
print('Loading iSLAT ' + iSLAT_version + ': Please Wait ...')

# Import necessary modules
import numpy as np
import pandas as pd
import os
import json

#from lmfit.models import GaussianModel
import tkinter as tk
from tkinter import filedialog
import ssl
import certifi
import datetime

context = ssl.create_default_context(cafile=certifi.where ())

from .iSLATFileHandling import load_user_settings, read_from_csv, read_default_csv, read_from_user_csv, read_default_molecule_parameters, read_initial_molecule_parameters, read_save_data

from .ir_model import *
from .COMPONENTS.chart_window import MoleculeSelector
from .COMPONENTS.Hitran_data import get_Hitran_data
from .COMPONENTS.partition_function_writer import write_partition_function
from .COMPONENTS.line_data_writer import write_line_data
from .COMPONENTS.slabfit_config import *
from .COMPONENTS.slabfit_loader import *
from .COMPONENTS.slabfit_runner import *
from .iSLATDefaultInputParms import *
#from .iSLATFileHandling import *
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
        self.directorypath = os.path.dirname(os.path.abspath(__file__))
        #print(f"iSLAT directory path: {self.directorypath}")
        self.create_folders()

        # Load settings
        #self.user_settings = self.load_user_settings()
        self.user_settings = load_user_settings()
        #self.update_default_molecule_parameters()
        #self.update_initial_molecule_parameters()
        self.initial_molecule_parameters = read_initial_molecule_parameters()
        self.molecules_data_default = read_default_molecule_parameters()
        
        self.mols = ["H2", "HD", "H2O", "H218O", "CO2", "13CO2", "CO", "13CO", "C18O", "CH4", "HCN", "H13CN", "NH3", "OH", "C2H2", "13CCH2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.basem = ["H2", "H2", "H2O", "H2O", "CO2", "CO2", "CO", "CO", "CO", "CH4", "HCN", "HCN", "NH3", "OH", "C2H2", "C2H2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.isot = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1]

        self.wavelength_range = wavelength_range
        self._display_range = (23.52, 25.41)

        self.min_vu = 1 / (self.wavelength_range[0] / 1E6) / 100.
        self.max_vu = 1 / (self.wavelength_range[1] / 1E6) / 100.

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

        self.molecules_dict.load_molecules_data(molecules_data=self.molecules_data_default,
                                                initial_molecule_parameters=self.initial_molecule_parameters,
                                                save_file_data = self.savedata,
                                                wavelength_range = self.wavelength_range, 
                                           intrinsic_line_width = intrinsic_line_width, 
                                           model_pixel_res = model_pixel_res, 
                                           model_line_width = model_line_width,
                                           dist = dist)
        
        # Initialize the active molecule based on user settings
        active_molecule_name = self.user_settings.get("default_active_molecule", "H2O")
        if active_molecule_name in self.molecules_dict:
            self._active_molecule = self.molecules_dict[active_molecule_name]
        else:
            print(f"Active molecule '{active_molecule_name}' not found in the dictionary. Defaulting to 'H2O'.")
            self._active_molecule = self.molecules_dict.get("H2O", None)

    def run(self):
        """
        Run the iSLAT application.
        This function starts the main event loop of the Tkinter application.
        """
        # Start the main event loop
        self.savedata = read_save_data()
        self.load_spectrum()
        self.init_molecules()
        #self.err_data = np.full_like(self.flux_data, np.nanmedian(self.flux_data)/100)
        self.init_gui()

    '''def check_HITRAN(self, print_statments = True):
        """ check_HITRAN(print_statments=True) checks if the HITRAN files are present and downloads them if necessary.
        If print_statments is True, it will print the status of the HITRAN files to the console."""
        if print_statments: print('\nChecking for HITRAN files: ...')
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
            print('Not the first startup and reload_default_files is False. Skipping HITRAN files download.')'''
    
    def check_HITRAN(self, print_statments = True):
        """ check_HITRAN(print_statments=True) checks if the HITRAN files are present and downloads them if necessary.
        If print_statments is True, it will print the status of the HITRAN files to the console."""
        if print_statments: print('\nChecking for HITRAN files: ...')
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
            if mol.is_visible:
                mol_flux = mol.get_flux(self.wave_data)
                summed_flux += mol_flux
        self.sum_spectrum_flux = summed_flux

    def run_single_slab_fit(self):
        loader = DataLoader(self.molecules_data_default)
        self.slab_model = ModelFitting(loader)
        result = self.slab_model.fit()
        return result.summary()

    def get_line_data_in_range(self, xmin, xmax):
        selected_mol = self.active_molecule
        if not selected_mol:
            return None
        lines_df = selected_mol.intensity.get_table
        subset = lines_df[(lines_df['lam'] >= xmin) & (lines_df['lam'] <= xmax)]
        if subset.empty:
            return None
        return (subset['lam'].values,
                subset['intens'].values,
                subset['e_up'].values,
                subset['a_stein'].values,
                subset['g_up'].values)
    
    @property
    def active_molecule(self):
        return self._active_molecule
    
    @active_molecule.setter
    def active_molecule(self, molecule):
        """
        Sets the active molecule based on the provided name or object.
        If the molecule is not found, it throws an error and does not update the active molecule.
        """
        try:
            if isinstance(molecule, Molecule):
                self._active_molecule = molecule
            elif isinstance(molecule, str):
                if molecule in self.molecules_dict:
                    self._active_molecule = self.molecules_dict[molecule]
                else:
                    raise ValueError(f"Molecule '{molecule}' not found in the dictionary.")
            else:
                raise TypeError("Active molecule must be a Molecule object or a string representing the molecule name.")
            
            if hasattr(self, "GUI") and hasattr(self.GUI, "plot"):
                self.GUI.plot.update_population_diagram()
                #self.GUI.plot.update_line_inspection_plot()

        except (ValueError, TypeError) as e:
            print(f"Error setting active molecule: {e}")
        
    @property
    def display_range(self):
        """tuple: Display range for the spectrum plot."""
        return self._display_range
    
    @display_range.setter
    def display_range(self, value):
        """
        Sets the display range for the spectrum plot.
        The value should be a tuple of two floats representing the start and end wavelengths.
        """
        if isinstance(value, tuple) and len(value) == 2:
            self._display_range = value
            if hasattr(self, "GUI") and hasattr(self.GUI, "plot"):
                self.GUI.plot.match_display_range()
        else:
            raise ValueError("Display range must be a tuple of two floats (start, end).")