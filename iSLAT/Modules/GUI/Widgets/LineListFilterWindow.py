"""
LineListFilterWindow — interactive line-list filter popup.

Opens a Toplevel window bound to a Molecule's MoleculeLineList.  The user
can apply numeric range filters, quantum-label filters, and vibrational-band
filters (e.g. "v1", "v1-0", "v2-1") using the LineListMaker fluent API.
Results can be exported to a CSV file.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from typing import Optional, TYPE_CHECKING

from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule

# ---------------------------------------------------------------------------
# LineListFilterWindow
# ---------------------------------------------------------------------------

class LineListFilterWindow(tk.Toplevel):
    """Interactive popup for filtering a molecule's line list.

    Parameters
    ----------
    parent : tk widget
        Parent window (typically the ControlPanel or root).
    mol_obj : Molecule
        The molecule whose ``line_list`` will be filtered.
    data_field : optional
        iSLAT data field for status messages (``insert_text`` interface).
    islat : optional
        iSLAT instance.  Required for "Apply to Molecule" and
        "Duplicate & Apply" GUI refresh.
    """

    def __init__(self, parent, mol_obj: "Molecule", data_field=None, islat=None):
        super().__init__(parent)
        self.mol_obj = mol_obj
        self.data_field = data_field
        self._islat = islat

        mol_name = getattr(mol_obj, "name", "Molecule")
        self.title(f"Filter Line List: {mol_name}")
        self.resizable(True, True)
        self._constrain_to_screen()

        # Build the maker from the molecule's line list
        line_list = mol_obj.line_list
        if line_list is None:
            self.destroy()
            raise ValueError(f"Molecule '{mol_name}' has no line list loaded.")
        self._maker = LineListMaker(line_list)
        self._total_lines = len(self._maker)

        # Build UI
        self._build_ui()
        self._refresh_status()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _constrain_to_screen(self):
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()
        win_w = min(520, int(screen_w * 0.45))
        win_h = min(680, int(screen_h * 0.80))
        pos_x = max((screen_w - win_w) // 2, 0)
        pos_y = max((screen_h - win_h) // 2, 0)
        self.geometry(f"{win_w}x{win_h}+{pos_x}+{pos_y}")
        self.minsize(420, 480)

    def _build_ui(self):
        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(0, weight=0)  # filter controls
        outer.rowconfigure(1, weight=1)  # active filters listbox
        outer.rowconfigure(2, weight=0)  # action bar

        # ── Filter controls ──────────────────────────────────────────
        ctrl_lf = ttk.LabelFrame(outer, text="Filters", padding=6)
        ctrl_lf.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        ctrl_lf.columnconfigure(1, weight=1)
        ctrl_lf.columnconfigure(2, weight=1)

        # Column headers
        ttk.Label(ctrl_lf, text="Property", font=("", 9, "bold")).grid(
            row=0, column=0, sticky="w", padx=4)
        ttk.Label(ctrl_lf, text="Min / Value", font=("", 9, "bold")).grid(
            row=0, column=1, sticky="w", padx=4)
        ttk.Label(ctrl_lf, text="Max / Option", font=("", 9, "bold")).grid(
            row=0, column=2, sticky="w", padx=4)
        ttk.Separator(ctrl_lf, orient="horizontal").grid(
            row=1, column=0, columnspan=3, sticky="ew", pady=2)

        self._entries = {}  # name -> (min_var, max_var) or special

        numeric_fields = [
            ("Wavelength (µm)",   "wavelength"),
            ("E_up (K)",          "eup"),
            ("E_low (K)",         "elow"),
            ("Einstein A (s⁻¹)", "astein"),
            ("g_up",              "gup"),
            ("g_low",             "glow"),
        ]
        for r, (label, key) in enumerate(numeric_fields, start=2):
            ttk.Label(ctrl_lf, text=label).grid(
                row=r, column=0, sticky="w", padx=4, pady=2)
            min_var = tk.StringVar()
            max_var = tk.StringVar()
            ttk.Entry(ctrl_lf, textvariable=min_var, width=10).grid(
                row=r, column=1, padx=4, pady=2, sticky="ew")
            ttk.Entry(ctrl_lf, textvariable=max_var, width=10).grid(
                row=r, column=2, padx=4, pady=2, sticky="ew")
            self._entries[key] = (min_var, max_var)

        sep_row = len(numeric_fields) + 2
        ttk.Separator(ctrl_lf, orient="horizontal").grid(
            row=sep_row, column=0, columnspan=3, sticky="ew", pady=4)

        # Quantum label rows
        for r_offset, (label, key) in enumerate([
            ("Quantum lev_up",  "lev_up"),
            ("Quantum lev_low", "lev_low"),
        ], start=sep_row + 1):
            ttk.Label(ctrl_lf, text=label).grid(
                row=r_offset, column=0, sticky="w", padx=4, pady=2)
            text_var = tk.StringVar()
            contains_var = tk.BooleanVar(value=True)
            ttk.Entry(ctrl_lf, textvariable=text_var, width=14).grid(
                row=r_offset, column=1, padx=4, pady=2, sticky="ew")
            ttk.Checkbutton(
                ctrl_lf, text="contains", variable=contains_var,
            ).grid(row=r_offset, column=2, padx=4, pady=2, sticky="w")
            self._entries[key] = (text_var, contains_var)

        vib_sep_row = sep_row + 3
        ttk.Separator(ctrl_lf, orient="horizontal").grid(
            row=vib_sep_row, column=0, columnspan=3, sticky="ew", pady=4)

        # Vibrational band row
        ttk.Label(ctrl_lf, text="Vib. Band").grid(
            row=vib_sep_row + 1, column=0, sticky="w", padx=4, pady=2)
        vib_var = tk.StringVar()
        vib_entry = ttk.Entry(ctrl_lf, textvariable=vib_var, width=10)
        vib_entry.grid(row=vib_sep_row + 1, column=1, padx=4, pady=2, sticky="ew")
        # Tooltip-style placeholder hint via focus events
        _hint = "e.g. v1, v1-0, v2-1"
        vib_entry.insert(0, _hint)
        vib_entry.config(foreground="grey")

        def _on_focus_in(e):
            if vib_var.get() == _hint:
                vib_entry.delete(0, tk.END)
                vib_entry.config(foreground="")

        def _on_focus_out(e):
            if not vib_var.get():
                vib_entry.insert(0, _hint)
                vib_entry.config(foreground="grey")

        vib_entry.bind("<FocusIn>", _on_focus_in)
        vib_entry.bind("<FocusOut>", _on_focus_out)

        ttk.Label(ctrl_lf, text="# modes").grid(
            row=vib_sep_row + 1, column=2, sticky="w", padx=(4, 0), pady=2)
        n_modes_var = tk.IntVar(value=3)
        ttk.Spinbox(ctrl_lf, from_=2, to=6, textvariable=n_modes_var,
                    width=4).grid(row=vib_sep_row + 1, column=2,
                                  padx=(52, 4), pady=2, sticky="w")
        self._entries["vib_band"] = (vib_var, n_modes_var)
        self._vib_hint = _hint

        # ── Active filters listbox ───────────────────────────────────
        af_lf = ttk.LabelFrame(outer, text="Active Filters", padding=6)
        af_lf.grid(row=1, column=0, sticky="nsew", pady=(0, 6))
        af_lf.columnconfigure(0, weight=1)
        af_lf.rowconfigure(0, weight=1)

        self._filter_listbox = tk.Listbox(af_lf, height=6, selectmode=tk.SINGLE)
        self._filter_listbox.grid(row=0, column=0, sticky="nsew")
        af_scroll = ttk.Scrollbar(af_lf, orient="vertical",
                                   command=self._filter_listbox.yview)
        af_scroll.grid(row=0, column=1, sticky="ns")
        self._filter_listbox.configure(yscrollcommand=af_scroll.set)

        ttk.Button(af_lf, text="Remove Selected",
                   command=self._remove_selected_filter).grid(
            row=1, column=0, columnspan=2, sticky="ew", pady=(4, 0))

        # ── Status bar ───────────────────────────────────────────────
        self._status_var = tk.StringVar()
        ttk.Label(outer, textvariable=self._status_var, relief="sunken",
                  anchor="w").grid(row=2, column=0, sticky="ew", pady=(0, 4))

        # ── Action buttons ───────────────────────────────────────────
        btn_frame = ttk.Frame(outer)
        btn_frame.grid(row=3, column=0, sticky="ew")
        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        btn_frame.columnconfigure(2, weight=1)

        ttk.Button(btn_frame, text="Apply Filters",
                   command=self._apply_filters).grid(
            row=0, column=0, padx=4, pady=4, sticky="ew")
        ttk.Button(btn_frame, text="Reset",
                   command=self._reset_all).grid(
            row=0, column=1, padx=4, pady=4, sticky="ew")
        ttk.Button(btn_frame, text="Export CSV",
                   command=self._export_csv).grid(
            row=0, column=2, padx=4, pady=4, sticky="ew")

        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=1, column=0, columnspan=3, sticky="ew", pady=(2, 0))

        apply_mol_frame = ttk.Frame(btn_frame)
        apply_mol_frame.grid(row=2, column=0, columnspan=3, sticky="ew")
        apply_mol_frame.columnconfigure(0, weight=1)
        apply_mol_frame.columnconfigure(1, weight=1)

        ttk.Button(apply_mol_frame, text="Apply to Molecule",
                   command=self._apply_to_molecule).grid(
            row=0, column=0, padx=4, pady=4, sticky="ew")
        ttk.Button(apply_mol_frame, text="Duplicate & Apply",
                   command=self._duplicate_and_apply).grid(
            row=0, column=1, padx=4, pady=4, sticky="ew")

    # ------------------------------------------------------------------
    # Logic
    # ------------------------------------------------------------------

    def _parse_float(self, s: str) -> Optional[float]:
        s = s.strip()
        if not s:
            return None
        try:
            return float(s)
        except ValueError:
            return None

    def _apply_filters(self):
        """Reset maker, then re-apply all non-empty fields."""
        self._maker.reset()

        # Numeric range filters
        numeric_map = {
            "wavelength": self._maker.filter_wavelength,
            "eup":        self._maker.filter_eup,
            "elow":       self._maker.filter_elow,
            "astein":     self._maker.filter_astein,
            "gup":        self._maker.filter_gup,
            "glow":       self._maker.filter_glow,
        }
        for key, fn in numeric_map.items():
            min_var, max_var = self._entries[key]
            lo = self._parse_float(min_var.get())
            hi = self._parse_float(max_var.get())
            if lo is not None or hi is not None:
                fn(min_val=lo, max_val=hi)

        # Quantum label filters
        for key in ("lev_up", "lev_low"):
            text_var, contains_var = self._entries[key]
            text = text_var.get().strip()
            if text:
                kwargs = {key: text, "contains": contains_var.get()}
                self._maker.filter_quantum(**kwargs)

        # Vibrational band filter
        vib_var, n_modes_var = self._entries["vib_band"]
        vib_spec = vib_var.get().strip()
        if vib_spec and vib_spec != self._vib_hint:
            try:
                n_modes = int(n_modes_var.get())
                self._maker.filter_vib_band(vib_spec, n_modes)
            except (ValueError, Exception) as e:
                messagebox.showerror(
                    "Band Parse Error",
                    f"Could not parse vibrational band '{vib_spec}':\n{e}",
                    parent=self,
                )
                self._maker.reset()
                return

        self._refresh_status()

    def _reset_all(self):
        """Clear all entries and reset the maker to original data."""
        self._maker.reset()
        # Clear numeric entries
        for key in ("wavelength", "eup", "elow", "astein", "gup", "glow"):
            min_var, max_var = self._entries[key]
            min_var.set("")
            max_var.set("")
        # Clear quantum entries
        for key in ("lev_up", "lev_low"):
            text_var, _ = self._entries[key]
            text_var.set("")
        # Clear vibrational band (restore hint)
        vib_var, _ = self._entries["vib_band"]
        vib_var.set("")
        self._refresh_status()

    def _remove_selected_filter(self):
        sel = self._filter_listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        try:
            self._maker.remove_filter(idx)
        except IndexError:
            pass
        self._refresh_status()

    def _refresh_status(self):
        current = len(self._maker)
        self._status_var.set(f"Lines: {current} / {self._total_lines}")
        # Rebuild listbox
        self._filter_listbox.delete(0, tk.END)
        for name, kwargs in self._maker.filters:
            kw_str = ", ".join(
                f"{k}={v!r}" for k, v in kwargs.items()
                if v is not None and k != "label"
            )
            display = f"{name}({kw_str})" if kw_str else name
            self._filter_listbox.insert(tk.END, display)

    # ------------------------------------------------------------------
    # Line-list assignment helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _assign_line_list(mol_obj, new_ll) -> None:
        """Replace *mol_obj*'s line list and invalidate all cached results.

        Setting ``mol_obj.lines`` directly bypasses the molecule's cache
        invalidation logic, so we must mark every dirty flag and clear the
        flux / intensity / spectrum caches explicitly.
        """
        mol_obj.lines = new_ll
        # Nullify already-computed results so they are recalculated
        mol_obj.intensity = None
        mol_obj.spectrum = None
        # Mark all dirty flags
        try:
            mol_obj._dirty_flags['intensity'] = True
            mol_obj._dirty_flags['spectrum'] = True
            mol_obj._dirty_flags['flux'] = True
        except (AttributeError, KeyError):
            pass
        # Clear the flux cache (dict)
        try:
            mol_obj._flux_cache.clear()
        except AttributeError:
            pass
        # Clear intensity / spectrum caches
        try:
            mol_obj._intensity_cache['data'] = None
            mol_obj._intensity_cache['hash'] = None
        except (AttributeError, KeyError):
            pass
        try:
            mol_obj._spectrum_cache['data'] = None
            mol_obj._spectrum_cache['hash'] = None
        except (AttributeError, KeyError):
            pass

    def _apply_to_molecule(self):
        """Replace the molecule's line list with the current filtered result."""
        if len(self._maker) == 0:
            messagebox.showwarning(
                "Empty Result",
                "The filtered line list is empty — not applied.",
                parent=self,
            )
            return
        mol_name = getattr(self.mol_obj, "name", "this molecule")
        confirm = messagebox.askyesno(
            "Apply to Molecule",
            f"Replace '{mol_name}' line list with {len(self._maker)} filtered lines?\n"
            "This overwrites the in-memory line list and cannot be undone.",
            parent=self,
        )
        if not confirm:
            return
        new_ll = self._maker.to_linelist()
        self._assign_line_list(self.mol_obj, new_ll)
        self._total_lines = len(new_ll)
        self._trigger_refresh(full_rebuild=False)
        self._refresh_status()
        if self.data_field is not None:
            self.data_field.insert_text(
                f"Applied filtered line list to '{mol_name}': {len(new_ll)} lines."
            )

    def _duplicate_and_apply(self):
        """Duplicate the molecule and apply the filtered line list to the copy."""
        if len(self._maker) == 0:
            messagebox.showwarning(
                "Empty Result",
                "The filtered line list is empty — not applied.",
                parent=self,
            )
            return
        if self._islat is None:
            messagebox.showerror(
                "No iSLAT Reference",
                "Cannot duplicate: no iSLAT instance was provided to this window.",
                parent=self,
            )
            return
        mol_name = getattr(self.mol_obj, "name", None)
        if not mol_name:
            messagebox.showerror("Error", "Molecule has no name.", parent=self)
            return
        new_name = self._islat.molecules_dict.duplicate(mol_name)
        if new_name is None:
            messagebox.showerror(
                "Duplicate Failed",
                f"Could not duplicate molecule '{mol_name}'.",
                parent=self,
            )
            return
        new_mol = self._islat.molecules_dict.get(new_name)
        if new_mol is not None:
            self._assign_line_list(new_mol, self._maker.to_linelist())
        self._trigger_refresh(full_rebuild=True)
        msg = (
            f"Created '{new_name}' with {len(self._maker)} filtered lines."
        )
        if self.data_field is not None:
            self.data_field.insert_text(msg)
        else:
            messagebox.showinfo("Done", msg, parent=self)

    def _trigger_refresh(self, full_rebuild: bool = False):
        """Attempt to refresh the iSLAT GUI after a line-list change."""
        try:
            if self._islat is None:
                return
            gui = getattr(self._islat, "GUI", None)
            if gui is None:
                return
            plot = getattr(gui, "plot", None)
            if plot is not None:
                plot.update_model_plot()
            cp = getattr(gui, "control_panel", None)
            if cp is not None:
                if full_rebuild:
                    cp.refresh_from_molecules_dict()
                else:
                    cp._update_molecule_parameter_fields()
        except Exception as e:
            print(f"LineListFilterWindow: GUI refresh error — {e}")

    def _export_csv(self):
        mol_name = getattr(self.mol_obj, "name", "molecule")
        path = filedialog.asksaveasfilename(
            parent=self,
            title="Export Filtered Line List",
            defaultextension=".csv",
            initialfile=f"{mol_name}_filtered.csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            out = self._maker.to_csv(Path(path))
            msg = f"Exported {len(self._maker)} lines → {out.name}"
            if self.data_field is not None:
                self.data_field.insert_text(msg)
            else:
                messagebox.showinfo("Export Complete", msg, parent=self)
        except Exception as e:
            messagebox.showerror("Export Error", str(e), parent=self)
