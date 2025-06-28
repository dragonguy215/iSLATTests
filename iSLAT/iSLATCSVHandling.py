import os
import csv
from .iSLATDefaultInputParms import molecules_data

save_folder = "SAVES"

def read_from_csv():
    filename = os.path.join(save_folder, "molsave.csv")
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_default_csv():
    filename = os.path.join(save_folder, "default.csv")
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data

def read_from_user_csv():
    filename = os.path.join(save_folder, "molecules_list.csv")
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                return [row for row in reader]
        except FileNotFoundError:
            pass
    return molecules_data