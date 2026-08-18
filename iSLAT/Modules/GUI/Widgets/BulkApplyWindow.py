"""BulkApplyWindow - popout window for bulk-applying properties to molecules.

Lets the user set one or more molecule properties (temperature, radius, column
density, ...) on many molecules at once, optionally restricting the change to only
the molecules that are currently visible in the control panel.

Usage
-----
    from iSLAT.Modules.GUI.Widgets.BulkApplyWindow import BulkApplyPropertiesWindow
    BulkApplyPropertiesWindow.open(parent_widget, islat_instance, control_panel, theme)
"""
from __future__ import annotations
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import messagebox
from typing import Any, Dict, List, Optional, Tuple
from iSLAT.Modules.DataProcessing.InstrumentalProfiles import PROFILE_DISPLAY_NAMES

# Ordered mapping of display label -> molecule attribute. Each property is a
# numeric physical parameter accepted by ``MoleculeDict.bulk_update_parameters``.
BULK_PROPERTIES: List[Tuple[str, str]] = [
    ("Temperature", "temp"),
    ("Radius", "radius"),
    ("Column Density", "n_mol"),
    ("Distance", "distance"),
    ("Instrumental FWHM", "fwhm"),
    ("Keplerian FWHM", "keplerian_fwhm"),
    ("Broadening", "broad"),
    ("RV Shift", "rv_shift"),
]

# Sentinel shown in choice drop-downs meaning "do not change this property".
UNCHANGED_CHOICE = "(unchanged)"

# Ordered mapping of display label -> (molecule attribute, {display -> value}).
# These are non-numeric properties selected from a drop-down.
BULK_CHOICE_PROPERTIES: List[Tuple[str, str, Dict[str, str]]] = [
    (
        "Inst. Profile",
        "instrumental_profile_key",
        {label: key for key, label in PROFILE_DISPLAY_NAMES.items()},
    ),
]

class BulkApplyPropertiesWindow:
    """Popout window that bulk-applies properties to molecules."""
    # Class-level registry: one window per (islat instance id)
    _instances: Dict[int, "BulkApplyPropertiesWindow"] = {}
    # Public factory
    @classmethod
    def open(
        cls,
        parent: tk.Widget,
        islat: Any,
        control_panel: Any,
        theme: Optional[Dict[str, Any]] = None,
    ) -> "BulkApplyPropertiesWindow":
        """Open (or raise) the bulk-apply window for *islat*."""
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
        self._win.title("Bulk Apply Properties")
        self._win.resizable(False, False)
        self._win.protocol("WM_DELETE_WINDOW", self._on_close)

        self._entry_vars: Dict[str, tk.StringVar] = {
            attr: tk.StringVar() for _label, attr in BULK_PROPERTIES
        }
        self._choice_vars: Dict[str, tk.StringVar] = {
            attr: tk.StringVar(value=UNCHANGED_CHOICE)
            for _label, attr, _options in BULK_CHOICE_PROPERTIES
        }
        self._visible_only_var = tk.BooleanVar(value=False)

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

        ttk.Label(
            content,
            text=f"Leave a field blank (or set to {UNCHANGED_CHOICE})\n"
            "to keep that property unchanged.",
            justify="left",
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))

        row = 0
        for row, (label, attr) in enumerate(BULK_PROPERTIES, start=1):
            ttk.Label(content, text=f"{label}:").grid(
                row=row, column=0, sticky="w", padx=(0, 6), pady=2
            )
            ttk.Entry(
                content, textvariable=self._entry_vars[attr], width=18
            ).grid(row=row, column=1, sticky="ew", pady=2)

        for label, attr, options in BULK_CHOICE_PROPERTIES:
            row += 1
            ttk.Label(content, text=f"{label}:").grid(
                row=row, column=0, sticky="w", padx=(0, 6), pady=2
            )
            ttk.Combobox(
                content,
                textvariable=self._choice_vars[attr],
                values=[UNCHANGED_CHOICE, *options.keys()],
                state="readonly",
                width=16,
            ).grid(row=row, column=1, sticky="ew", pady=2)

        ttk.Checkbutton(
            content,
            text="Apply to visible molecules only",
            variable=self._visible_only_var,
        ).grid(
            row=row + 1,
            column=0,
            columnspan=2,
            sticky="w",
            pady=(8, 8),
        )

        btn_frame = ttk.Frame(content)
        btn_frame.grid(row=row + 2, column=0, columnspan=2, sticky="e")
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
    def _on_apply(self) -> None:
        mol_dict = getattr(self._islat, "molecules_dict", None)
        if not mol_dict:
            self._on_close()
            return

        raw_values = {attr: var.get() for attr, var in self._entry_vars.items()}
        param_dict, invalid = build_parameter_dict(raw_values)
        param_dict.update(
            build_choice_dict(
                {attr: var.get() for attr, var in self._choice_vars.items()}
            )
        )

        if invalid:
            invalid_labels = [
                label for label, attr in BULK_PROPERTIES if attr in invalid
            ]
            messagebox.showerror(
                "Invalid value",
                "Please enter numeric values for:\n" + ", ".join(invalid_labels),
                parent=self._win,
            )
            return

        if not param_dict:
            messagebox.showinfo(
                "Nothing to apply",
                "Enter a value for at least one property.",
                parent=self._win,
            )
            return

        visible_only = self._visible_only_var.get()
        target_names = resolve_target_names(mol_dict, visible_only)

        if not target_names:
            messagebox.showinfo(
                "No molecules",
                "There are no visible molecules to apply properties to."
                if visible_only
                else "There are no molecules to apply properties to.",
                parent=self._win,
            )
            return

        try:
            mol_dict.bulk_update_parameters(param_dict, molecules=target_names)
        except Exception as exc:
            messagebox.showerror(
                "Error", f"Failed to apply properties:\n{exc}", parent=self._win
            )
            return

        self._refresh_gui()
        self._on_close()

    def _refresh_gui(self) -> None:
        """Refresh the plot and control panel after a bulk change."""
        gui = getattr(self._islat, "GUI", None)
        plot = getattr(gui, "plot", None) if gui is not None else None
        if plot is not None and hasattr(plot, "update_model_plot"):
            try:
                plot.update_model_plot()
            except Exception as exc:
                print(f"BulkApplyPropertiesWindow: plot refresh warning — {exc}")

        if self._control_panel is not None and hasattr(
            self._control_panel, "refresh_from_molecules_dict"
        ):
            try:
                self._control_panel.refresh_from_molecules_dict()
            except Exception as exc:
                print(
                    f"BulkApplyPropertiesWindow: control-panel refresh warning — {exc}"
                )

    def _on_close(self) -> None:
        BulkApplyPropertiesWindow._instances.pop(id(self._islat), None)
        try:
            self._win.destroy()
        except Exception:
            pass

def build_parameter_dict(
    raw_values: Dict[str, str],
) -> Tuple[Dict[str, float], List[str]]:
    """Convert raw entry strings into a parameter dict.

    Blank/whitespace-only fields are skipped. Non-empty fields are parsed as
    floats. Returns ``(param_dict, invalid_attrs)`` where ``invalid_attrs`` lists
    the attributes whose values could not be parsed as numbers.
    """
    param_dict: Dict[str, float] = {}
    invalid: List[str] = []
    for attr, raw in raw_values.items():
        text = (raw or "").strip()
        if not text:
            continue
        try:
            param_dict[attr] = float(text)
        except (TypeError, ValueError):
            invalid.append(attr)
    return param_dict, invalid

def build_choice_dict(raw_values: Dict[str, str]) -> Dict[str, str]:
    """Convert drop-down selections into a parameter dict.

    Selections equal to :data:`UNCHANGED_CHOICE` (or blank) are skipped; the
    remaining display labels are mapped back to their stored attribute values.
    """
    option_maps = {attr: options for _label, attr, options in BULK_CHOICE_PROPERTIES}
    param_dict: Dict[str, str] = {}
    for attr, raw in raw_values.items():
        text = (raw or "").strip()
        if not text or text == UNCHANGED_CHOICE:
            continue
        value = option_maps.get(attr, {}).get(text)
        if value is not None:
            param_dict[attr] = value
    return param_dict

def resolve_target_names(mol_dict: Any, visible_only: bool) -> List[str]:
    """Return the molecule names to update.

    When *visible_only* is True, only currently-visible molecules are returned;
    otherwise every molecule in *mol_dict* is returned.
    """
    if visible_only:
        if hasattr(mol_dict, "get_visible_molecules"):
            return list(mol_dict.get_visible_molecules())
        return [
            name
            for name in mol_dict.keys()
            if getattr(mol_dict[name], "is_visible", False)
        ]
    return list(mol_dict.keys())