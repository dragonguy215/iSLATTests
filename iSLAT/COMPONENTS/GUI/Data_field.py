import tkinter as tk
from tkinter import ttk

class DataField:
    """
    A class to represent a text area in a GUI application.
    No layout here — the caller does grid/pack for .frame.
    """
    def __init__(self, name: str, value: any, master: tk.Widget):
        self.name = name
        self.value = value

        # Entire frame to hold label + text + scrollbar
        self.frame = ttk.Frame(master)

        # Label
        self.label = ttk.Label(self.frame, text=self.name)
        self.label.pack(fill="x")

        # Text widget with scrollbar
        self.text = tk.Text(self.frame, height=20, width=60, wrap="word")
        self.scrollbar = ttk.Scrollbar(self.frame, orient="vertical", command=self.text.yview)
        self.text.configure(yscrollcommand=self.scrollbar.set)

        # Place text & scrollbar
        self.text.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

    def insert_text(self, content):
        self.text.insert("end", str(content) + "\n")
        self.text.see("end")

    def clear(self):
        self.text.delete("1.0", "end")

    def delete(self, start="1.0", end="end"):
        self.text.delete(start, end)

    def __repr__(self):
        return f"DataField(name={self.name}, value={self.value})"