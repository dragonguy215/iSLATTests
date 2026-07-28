"""SortMoleculesWindow — popout window for sorting the control-panel molecules.

Lets the user reorder the molecules shown in the control panel by any of their
numeric parameters (temperature, radius, column density, …) or by name, in
ascending or descending order.

Usage
-----
    from iSLAT.Modules.GUI.Widgets.SortMoleculesWindow import SortMoleculesWindow
    SortMoleculesWindow.open(parent_widget, islat_instance, control_panel, theme)
"""
from __future__ import annotations

import tkinter as tk
import tkinter.ttk as ttk
from typing import Any, Dict, List, Optional, Tuple

# Ordered mapping of display label -> (molecule attribute, is_numeric).
# ``is_numeric`` selects the sort key/normalisation used for that property.
SORT_PROPERTIES: List[Tuple[str, str, bool]] = [
    ("Name", "name", False),
    ("Temperature", "temp", True),
    ("Radius", "radius", True),
    ("Column Density", "n_mol", True),
    ("Distance", "distance", True),
    ("Instrumental FWHM", "fwhm", True),
    ("Keplerian FWHM", "keplerian_fwhm", True),
    ("Broadening", "broad", True),
    ("RV Shift", "rv_shift", True),
]

class SortMoleculesWindow:
    """Singleton popout window that sorts the control-panel molecules."""

    # Class-level registry: one window per (islat instance id)
    _instances: Dict[int, "SortMoleculesWindow"] = {}

    # ------------------------------------------------------------------
    # Public factory
    # ------------------------------------------------------------------
    @classmethod
    def open(
        cls,
        parent: tk.Widget,
        islat: Any,
        control_panel: Any,
        theme: Optional[Dict[str, Any]] = None,
    ) -> "SortMoleculesWindow":
        """Open (or raise) the sort window for *islat*."""
        key = id(islat)
        existing = cls._instances.get(key)
        if existing is not None:
            try:
                if existing._win.winfo_exists():
                    existing._win.lift()
                    existing._win.focus_force()
                    return existing
            except Exception:
                pass
            del cls._instances[key]

        instance = cls(parent, islat, control_panel, theme or {})
        cls._instances[key] = instance
        return instance

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------
    def __init__(
        self,
        parent: tk.Widget,
        islat: Any,
        control_panel: Any,
        theme: Dict[str, Any],
    ) -> None:
        self._islat = islat
        self._control_panel = control_panel
        self._theme = theme

        self._win = tk.Toplevel(parent)
        self._win.title("Sort Molecules")
        self._win.resizable(False, False)
        self._win.protocol("WM_DELETE_WINDOW", self._on_close)

        self._property_var = tk.StringVar(value=SORT_PROPERTIES[0][0])
        self._order_var = tk.StringVar(value="ascending")

        self._build_ui()

        # Centre over the parent window.
        self._win.update_idletasks()
        x = self._win.winfo_screenwidth() // 2 - self._win.winfo_width() // 2
        y = self._win.winfo_screenheight() // 2 - self._win.winfo_height() // 2
        self._win.geometry(f"+{x}+{y}")

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        win = self._win

        content = ttk.Frame(win)
        content.pack(fill="both", expand=True, padx=10, pady=10)

        ttk.Label(content, text="Sort by:").grid(
            row=0, column=0, sticky="w", padx=(0, 6), pady=(0, 6)
        )
        self._property_combobox = ttk.Combobox(
            content,
            textvariable=self._property_var,
            values=[label for label, _attr, _num in SORT_PROPERTIES],
            state="readonly",
            width=20,
        )
        self._property_combobox.grid(row=0, column=1, sticky="ew", pady=(0, 6))

        order_frame = ttk.LabelFrame(content, text="Order")
        order_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        ttk.Radiobutton(
            order_frame,
            text="Ascending",
            value="ascending",
            variable=self._order_var,
        ).pack(side="left", padx=6, pady=4)
        ttk.Radiobutton(
            order_frame,
            text="Descending",
            value="descending",
            variable=self._order_var,
        ).pack(side="left", padx=6, pady=4)

        btn_frame = ttk.Frame(content)
        btn_frame.grid(row=2, column=0, columnspan=2, sticky="e")
        ttk.Button(btn_frame, text="Apply", command=self._on_apply).pack(
            side="right", padx=(4, 0)
        )
        ttk.Button(btn_frame, text="Cancel", command=self._on_close).pack(
            side="right", padx=4
        )

        content.columnconfigure(1, weight=1)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------
    def _selected_property(self) -> Tuple[str, bool]:
        """Return (attribute_name, is_numeric) for the chosen property."""
        label = self._property_var.get()
        for lbl, attr, is_numeric in SORT_PROPERTIES:
            if lbl == label:
                return attr, is_numeric
        return SORT_PROPERTIES[0][1], SORT_PROPERTIES[0][2]

    def _on_apply(self) -> None:
        mol_dict = getattr(self._islat, "molecules_dict", None)
        if not mol_dict:
            self._on_close()
            return

        attr, is_numeric = self._selected_property()
        descending = self._order_var.get() == "descending"
        sorted_names = sort_molecule_names(mol_dict, attr, is_numeric, descending)

        try:
            mol_dict.reorder(sorted_names)
        except Exception:
            self._on_close()
            return

        if self._control_panel is not None and hasattr(
            self._control_panel, "refresh_from_molecules_dict"
        ):
            self._control_panel.refresh_from_molecules_dict()

        self._on_close()

    def _on_close(self) -> None:
        SortMoleculesWindow._instances.pop(id(self._islat), None)
        try:
            self._win.destroy()
        except Exception:
            pass

def sort_molecule_names(
    mol_dict: Any,
    attr: str,
    is_numeric: bool,
    descending: bool,
) -> List[str]:
    """Return the molecule names sorted by *attr*.

    Numeric properties sort by float value (missing/invalid values treated as
    0.0); the name property sorts case-insensitively. The returned list always
    contains exactly the same names as ``mol_dict``.
    """
    names = list(mol_dict.keys())

    def numeric_key(name: str) -> float:
        try:
            return float(getattr(mol_dict[name], attr, 0.0))
        except (TypeError, ValueError):
            return 0.0

    def text_key(name: str) -> str:
        value = getattr(mol_dict[name], attr, name)
        return str(value).casefold()

    key = numeric_key if is_numeric else text_key
    return sorted(names, key=key, reverse=descending)