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

from .iSLATFileHandling import load_user_settings, read_default_molecule_parameters, read_initial_molecule_parameters, read_save_data, read_HITRAN_data

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
        #self.directorypath = os.path.dirname(os.path.abspath(__file__))
        #print(f"iSLAT directory path: {self.directorypath}")
        self._hitran_data = {}
        self.create_folders()

        # Load settings
        #self.user_settings = self.load_user_settings()
        self.user_settings = load_user_settings()
        #self.update_default_molecule_parameters()
        #self.update_initial_molecule_parameters()
        self.initial_molecule_parameters = read_initial_molecule_parameters()
        self.molecules_parameters_default = read_default_molecule_parameters()
        
        self.mols = ["H2", "HD", "H2O", "H218O", "CO2", "13CO2", "CO", "13CO", "C18O", "CH4", "HCN", "H13CN", "NH3", "OH", "C2H2", "13CCH2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.basem = ["H2", "H2", "H2O", "H2O", "CO2", "CO2", "CO", "CO", "CO", "CH4", "HCN", "HCN", "NH3", "OH", "C2H2", "C2H2", "C2H4", "C4H2", "C2H6", "HC3N"]
        self.isot = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1]

        self.wavelength_range = wavelength_range
        self._display_range = (23.52, 25.41)

        self.min_vu = 1 / (self.wavelength_range[0] / 1E6) / 100.
        self.max_vu = 1 / (self.wavelength_range[1] / 1E6) / 100.

        # Check for HITRAN files and download if necessary
        self.check_HITRAN()
        self.load_default_HITRAN_data()

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
        if not hasattr(self, "molecules_dict"):
            self.molecules_dict = MoleculeDict()
        #self.molecules_dict = MoleculeDict()
        #print("Hey whats up heres that dict", self.molecules_dict)

        new_molecules = []
        for mol_name, mol_data in self._hitran_data.items():
            if mol_name not in self.molecules_dict:
                new_molecule = Molecule(
                    name=mol_name,
                    filepath=mol_data.get("file_path", None),
                    initial_molecule_parameters=self.initial_molecule_parameters.get(mol_name, self.molecules_parameters_default),
                    hitran_data=mol_data,
                    intrinsic_line_width=intrinsic_line_width,
                    wavelength_range=self.wavelength_range,
                    model_pixel_res=model_pixel_res,
                    model_line_width=model_line_width,
                    distance=dist
                )
                new_molecules.append(new_molecule)

        #print("Hey whats up heres that dict again man", self.molecules_dict)

        if new_molecules:
            #print("Adding new molecules:")
            #print(new_molecules)
            self.molecules_dict.add_molecules(new_molecules)
        
        #print("Hey whats up heres that dict again again man", self.molecules_dict)

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
        self.init_molecules()
        self.load_spectrum()
        #self.err_data = np.full_like(self.flux_data, np.nanmedian(self.flux_data)/100)
        self.init_gui()
    
    def load_HITRAN_data(self, files=None):
        """
        Loads HITRAN data for the specified molecules or files.
        If no files are provided, it opens a file dialog to select files.
        The data is stored in self.hitran_data.
        """
        if files is None:
            filetypes = [('PAR Files', '*.par')]
            hitran_directory = os.path.abspath("HITRANdata")
            files = filedialog.askopenfilenames(title='Choose HITRAN Data Files', filetypes=filetypes, initialdir=hitran_directory)

        if not files:
            print("No files selected.")
            return

        hitran_data_list = []
        for file_path in files:
            mol_name = os.path.basename(file_path).split('_')[-1].split('.')[0]  # Extract molecule name from file name
            if not os.path.exists(file_path):
                print(f"WARNING: HITRAN file not found at {file_path}")
                continue

            lines = read_HITRAN_data(file_path)
            if lines:
                print(f"Preparing HITRAN data for {mol_name} with {len(lines)} lines.")
                try:
                    index = self.mols.index(mol_name)
                    base_molecule = self.basem[index]
                    isotope = self.isot[index]
                except ValueError:
                    base_molecule = None
                    isotope = None
                    print(f"WARNING: Molecule {mol_name} not found in predefined lists. Base molecule and isotope set to None.")

                hitran_data_list.append({
                    "lines": lines,
                    "base_molecule": base_molecule,
                    "isotope": isotope,
                    "file_path": file_path
                })
            else:
                print(f"WARNING: HITRAN file for {mol_name} could not be parsed.")

        if hitran_data_list:
            self.update_hitran_data_from_list(hitran_data_list)

    def check_HITRAN(self):
        """
        Checks that all expected HITRAN files are present,
        loads them using read_HITRAN_data, and optionally
        stores them in self.hitran_data.
        """
        print("\nChecking HITRAN files:")

        if self.user_settings["first_startup"] or self.user_settings["reload_default_files"]:
            print('First startup or reload_default_files is True. Downloading default HITRAN files ...')
            self.hitran_data = {}  # option: central dict holding lines by molecule name

            for mol, bm, iso in zip(self.mols, self.basem, self.isot):
                hitran_file = f"HITRANdata/data_Hitran_2020_{mol}.par"
                if not os.path.exists(hitran_file):
                    print(f"WARNING: HITRAN file for {mol} not found at {hitran_file}")
                    self.hitran_data[mol] = []
                    continue

                lines = read_HITRAN_data(hitran_file)
                if lines:
                    #print(f"Loaded HITRAN file for {mol}: {len(lines)} lines.")
                    self.hitran_data[mol] = {"lines": lines, "base_molecule": bm, "isotope": iso, "file_path": hitran_file}
                else:
                    #print(f"WARNING: HITRAN file for {mol} could not be parsed.")
                    self.hitran_data[mol] = []
        else:
            print('Not the first startup and reload_default_files is False. Skipping HITRAN files download.')

        print("Finished HITRAN file check.\n")

    def load_default_HITRAN_data(self, reset=False):
        """
        Loads default HITRAN data for the predefined molecules.
        """
        print("Loading default HITRAN data...")
        # Do not clear self.hitran_data here!
        if reset:
            #self._hitran_data = {}
            if hasattr(self, "molecules_dict"):
                self.molecules_dict.clear()
            self.hitran_data = {}
        valid_molecules = [
            mol for mol in self.mols
            if mol in self.initial_molecule_parameters
        ]
        hitran_files = [
            f"HITRANdata/data_Hitran_2020_{mol}.par"
            for mol in valid_molecules
        ]
        self.load_HITRAN_data(files=hitran_files)

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
    
    @property
    def hitran_data(self):
        """dict: Dictionary containing HITRAN data for molecules."""
        return self._hitran_data
    
    @hitran_data.setter
    def hitran_data(self, value):
        """
        Updates the HITRAN data and reloads relevant values when a new molecule is added.
        This includes updating the molecule table, spectrum data, and any dependent GUI components.
        """
        self._hitran_data = value

        # Reload molecule data
        self.init_molecules()

        # Update GUI components if they exist
        if hasattr(self, "GUI"):
            if hasattr(self.GUI, "molecule_table"):
                self.GUI.molecule_table.update_table()
            if hasattr(self.GUI, "control_panel"):
                self.GUI.control_panel.reload_molecule_dropdown()
            if hasattr(self.GUI, "plot"):
                self.GUI.plot.update_population_diagram()
                self.GUI.plot.update_line_inspection_plot()

    def update_hitran_data_from_list(self, hitran_data_list):
        for hitran_entry in hitran_data_list:
            if not isinstance(hitran_entry, dict):
                print("Invalid HITRAN data entry. Skipping...")
                continue
            mol_name = hitran_entry.get("base_molecule")
            if not mol_name:
                print("Base molecule name missing in HITRAN data entry. Skipping...")
                continue
            self._hitran_data[mol_name] = hitran_entry  # accumulate

        # Refresh molecules & GUI
        self.init_molecules()
        if hasattr(self, "GUI"):
            if hasattr(self.GUI, "molecule_table"):
                self.GUI.molecule_table.update_table()
            if hasattr(self.GUI, "control_panel"):
                self.GUI.control_panel.reload_molecule_dropdown()
            if hasattr(self.GUI, "plot"):
                self.GUI.plot.update_population_diagram()
                self.GUI.plot.update_line_inspection_plot()