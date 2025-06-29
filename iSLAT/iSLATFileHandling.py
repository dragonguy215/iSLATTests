import os
import csv
import json
import pandas as pd
from .iSLATDefaultInputParms import molecules_data

save_folder_path = "SAVES"
user_configuration_file_path = config_file_path = "CONFIG"
theme_file_path = "CONFIG/GUIThemes"
user_configuration_file_name = "UserSettings.json"

molsave_file_name = "molsave.csv"
defaults_file_name = "default.csv"
molecule_list_file_name = "molecules_list.csv"

default_molecule_parameters_file_name = "DefaultMoleculeParameters.json"
default_initial_parameters_file_name = "DefaultMoleculeParameters.json"

line_saves_file_name = "saved_lines.csv"

def load_user_settings(file_path=user_configuration_file_path, file_name=user_configuration_file_name, theme_file_path=theme_file_path):
    """ load_user_settings() loads the user settings from the UserSettings.json file."""
    file = os.path.join(file_path, file_name)
    if os.path.exists(file):
        with open(file, 'r') as f:
            user_settings = json.load(f)
    else:
        # If the file does not exist, return default settings and save them as a new json file
        default_settings = {
            "first_startup": True,
            "reload_default_files": True,
            "theme": "LightTheme"
        }
        with open(file, 'w') as f:
            json.dump(default_settings, f, indent=4)
        user_settings = default_settings
    
    # append theme information to the user settings dictonary
    theme_file = f"{theme_file_path}/{user_settings['theme']}.json"
    if os.path.exists(theme_file):
        with open(theme_file, 'r') as f:
            theme_settings = json.load(f)
        user_settings["theme"] = theme_settings
    return user_settings

def read_from_csv(file_path=save_folder_path, file_name=molsave_file_name):
    file = os.path.join(file_path, file_name)
    if os.path.exists(file):
        try:
            with open(file, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_default_csv(file_path=save_folder_path, file_name=defaults_file_name):
    file = os.path.join(file_path, file_name)
    if os.path.exists(file):
        try:
            with open(file, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_from_user_csv(file_path=save_folder_path, file_name=molecule_list_file_name):
    file = os.path.join(file_path, file_name)
    if os.path.exists(file):
        try:
            with open(file, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_default_molecule_parameters(file_path=config_file_path, file_name=default_molecule_parameters_file_name):
    """
    read_default_molecule_parameters() updates the default molecule parameters from the DefaultMoleculeParameters.json file.
    """
    file = os.path.join(file_path, file_name)
    with open(file, 'r') as f:
        default_molecule_parameters = json.load(f)["default_initial_params"]
    return default_molecule_parameters

def read_initial_molecule_parameters(file_path=config_file_path, file_name=default_initial_parameters_file_name):
    """
    read_initial_molecule_parameters() updates the initial molecule parameters from the DefaultMoleculeParameters.json file.
    """
    file = os.path.join(file_path, file_name)
    with open(file, 'r') as f:
        initial_molecule_parameters = json.load(f)["initial_parameters"]
    return initial_molecule_parameters

def read_save_data(file_path = save_folder_path, file_name=molecule_list_file_name):
    """
    read_save_data() loads the save data from the SAVES folder.
    It returns a list of dictionaries with the save data.
    """
    #save_file = os.path.join("SAVES", "molecules_list.csv")
    file = os.path.join(file_path, file_name)
    if os.path.exists(file):
        try:
            df = pd.read_csv(file)
            savedata = {row['Molecule Name']: {col: row[col] for col in df.columns if col != 'Molecule Name'} for _, row in df.iterrows()}
            return savedata
        except Exception as e:
            print(f"Error reading save file: {e}")
            savedata = {}
            return savedata
    else:
        print("No save file found.")
        savedata = {}
        return savedata

def read_HITRAN_data(file_path):
    """
    read_HITRAN_data(file_path) reads the HITRAN .par file at the given path.
    Returns the contents as a list of lines (or processes to DataFrame if needed).
    """
    if not os.path.exists(file_path):
        #print(f"HITRAN file '{file_path}' does not exist.")
        return []

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        #print(f"Successfully read HITRAN data from {file_path}")
        return lines
    except Exception as e:
        #print(f"Failed to read HITRAN file '{file_path}': {e}")
        return []

def read_line_saves(file_path=save_folder_path, file_name=line_saves_file_name):
    filename = os.path.join(file_path, file_name)
    if os.path.exists(file_path):
        try:
            with open(filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return []

def save_line(line_info, file_path=save_folder_path, file_name=line_saves_file_name):
    """Save a line to the line saves file."""
    filename = os.path.join(file_path, file_name)
    #print(f"Saving line to {filename}")
    df = pd.DataFrame([line_info])
    
    # Ensure the directory exists
    os.makedirs(file_path, exist_ok=True)
    
    # Save the line to the CSV file
    df.to_csv(filename, mode='a', header=not os.path.exists(filename), index=False)
    print(f"Saved line at ~{line_info['lam']:.4f} μm to {filename}")