import tkinter as tk
from tkinter import ttk

class DataField:
    """
    A class to represent a data field in a GUI application.
    
    Attributes:
        name (str): The name of the data field.
        value (any): The value of the data field.
    """

    def __init__(self, name: str, value: any, window: tk.Tk):
        """
        Initializes the DataField with a name and value.

        Args:
            name (str): The name of the data field.
            value (any): The value of the data field.
        """
        # create buttons for top of GUI
        self.name = name
        self.value = value

    def __repr__(self):
        return f"DataField(name={self.name}, value={self.value})"