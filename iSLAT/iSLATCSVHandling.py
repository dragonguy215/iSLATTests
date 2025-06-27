import os
import csv

# read more molecules if saved by the user in a previous iSLAT session
def read_from_csv():
    #global file_name
    filename = os.path.join(save_folder, f"{file_name}-molsave.csv")

    if os.path.exists (filename):
        try:
            with open (filename, 'r') as csvfile:
                reader = csv.reader (csvfile)
                next (reader)  # Skip the header row
                return [tuple (row[:3]) for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_default_csv():
    #global file_name
    filename = os.path.join(save_folder, f"default.csv")

    if os.path.exists (filename):
        try:
            with open (filename, 'r') as csvfile:
                reader = csv.reader (csvfile)
                next (reader)  # Skip the header row
                return [tuple (row[:3]) for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

# read more molecules if saved by the user in a previous iSLAT session
def read_from_user_csv():
    #global file_name
    filename = os.path.join(save_folder, f"molecules_list.csv")

    if os.path.exists (filename):
        try:
            with open (filename, 'r') as csvfile:
                reader = csv.reader (csvfile)
                next (reader)  # Skip the header row
                return [tuple (row[:3]) for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data