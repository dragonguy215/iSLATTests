"""
LineListViewWindow — read-only popout table for a molecule's line list.

Opens a Toplevel window displaying the MoleculeLineList as a sortable,
searchable ttk.Treeview table.  Columns are auto-discovered from the
DataFrame returned by ``MoleculeLineList.get_pandas_table()``.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk, filedialog
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_FLOAT_FMT = "{:.6g}"   # compact scientific/decimal formatting for floats


def _fmt(value) -> str:
    """Format a single cell value for display."""
    try:
        import math
        if isinstance(value, float):
            if math.isnan(value):
                return "NaN"
            return _FLOAT_FMT.format(value)
        if hasattr(value, "item"):          # numpy scalar → Python primitive
            return _fmt(value.item())
    except Exception:
        pass
    return str(value)


# ---------------------------------------------------------------------------
# LineListViewWindow
# ---------------------------------------------------------------------------

class LineListViewWindow(tk.Toplevel):
    """Read-only popout table showing a molecule's line list DataFrame.

    Parameters
    ----------
    parent : tk.Widget
        Parent window (typically the ControlPanel or root).
    mol_obj : Molecule
        The molecule whose ``line_list`` will be displayed.
    data_field : optional
        iSLAT data field for status messages (``insert_text`` interface).
    """

    def __init__(self, parent, mol_obj: "Molecule", data_field=None):
        super().__init__(parent)
        self.mol_obj = mol_obj
        self.data_field = data_field

        mol_name = getattr(mol_obj, "name", "Molecule")
        self.title(f"Line List: {mol_name}")
        self.resizable(True, True)
        self._constrain_to_screen()

        self._all_rows: list[tuple] = []    # full unfiltered row data
        self._columns: list[str] = []
        self._sort_col: str | None = None
        self._sort_reverse: bool = False

        self._build_ui()
        self._load_data()

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def _constrain_to_screen(self) -> None:
        sw = self.winfo_screenwidth()
        sh = self.winfo_screenheight()
        w = min(1100, int(sw * 0.80))
        h = min(700, int(sh * 0.75))
        x = max((sw - w) // 2, 0)
        y = max((sh - h) // 2, 0)
        self.geometry(f"{w}x{h}+{x}+{y}")
        self.minsize(500, 300)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        # ── Top bar: search ────────────────────────────────────────────
        top = ttk.Frame(self, padding=(6, 4))
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(1, weight=1)

        ttk.Label(top, text="Search:").grid(row=0, column=0, padx=(0, 4))
        self._search_var = tk.StringVar()
        self._search_var.trace_add("write", self._on_search_change)
        search_entry = ttk.Entry(top, textvariable=self._search_var)
        search_entry.grid(row=0, column=1, sticky="ew", padx=(0, 4))

        try:
            from iSLAT.Modules.GUI.Tooltips import CreateToolTip
            CreateToolTip(
                search_entry,
                "Filter rows by typing any text.\n"
                "All columns are searched (case-insensitive).",
            )
        except Exception:
            pass

        clear_btn = ttk.Button(top, text="✕", width=3, command=self._clear_search)
        clear_btn.grid(row=0, column=2, padx=(0, 4))

        refresh_btn = ttk.Button(top, text="Refresh", command=self._refresh)
        refresh_btn.grid(row=0, column=3, padx=(0, 8))
        try:
            from iSLAT.Modules.GUI.Tooltips import CreateToolTip
            CreateToolTip(
                refresh_btn,
                "Reload the line list from the molecule.\n"
                "Useful if the line list has been filtered or changed\n"
                "since this window was opened.",
            )
        except Exception:
            pass

        self._export_btn = ttk.Button(top, text="Export CSV…", command=self._export_csv)
        self._export_btn.grid(row=0, column=4)
        try:
            from iSLAT.Modules.GUI.Tooltips import CreateToolTip
            self._export_tooltip = CreateToolTip(
                self._export_btn,
                "Export selected rows to CSV.\n"
                "If no rows are selected, all visible rows are exported.\n"
                "Ctrl/Cmd+click or Shift+click to select multiple rows.\n"
                "Keyboard shortcut: Ctrl+S / Cmd+S.",
            )
        except Exception:
            pass

        # ── Treeview ───────────────────────────────────────────────────
        table_frame = ttk.Frame(self)
        table_frame.grid(row=1, column=0, sticky="nsew", padx=4, pady=(0, 2))
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        self._tree = ttk.Treeview(table_frame, show="headings", selectmode="extended")
        self._tree.grid(row=0, column=0, sticky="nsew")
        self._tree.bind("<<TreeviewSelect>>", self._on_selection_change)

        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self._tree.yview)
        vsb.grid(row=0, column=1, sticky="ns")
        hsb = ttk.Scrollbar(table_frame, orient="horizontal", command=self._tree.xview)
        hsb.grid(row=1, column=0, sticky="ew")
        self._tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # ── Status bar ─────────────────────────────────────────────────
        self._status_var = tk.StringVar(value="Loading…")
        status = ttk.Label(self, textvariable=self._status_var, anchor="w")
        status.grid(row=2, column=0, sticky="ew", padx=6, pady=(0, 4))

        # Ctrl+S / Cmd+S → export
        self.bind("<Control-s>", lambda _e: self._export_csv())
        self.bind("<Command-s>", lambda _e: self._export_csv())

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _refresh(self) -> None:
        """Reload the line list from the molecule, preserving the current search."""
        current_search = self._search_var.get()
        self._all_rows = []
        self._columns = []
        self._sort_col = None
        self._sort_reverse = False
        self._tree.delete(*self._tree.get_children())
        self._status_var.set("Refreshing…")
        self.update_idletasks()
        self._load_data()
        # Re-apply any active search after reload
        if current_search:
            self._search_var.set(current_search)

    def _load_data(self) -> None:
        """Fetch the line list DataFrame and populate the Treeview."""
        try:
            ll = self.mol_obj.line_list
        except Exception as exc:
            self._status_var.set(f"Error accessing line list: {exc}")
            return

        if ll is None:
            self._status_var.set("No line list loaded for this molecule.")
            return

        try:
            df = ll.get_pandas_table()
        except Exception as exc:
            self._status_var.set(f"Error building table: {exc}")
            return

        if df is None or df.empty:
            self._status_var.set("Line list is empty.")
            return

        self._columns = list(df.columns)
        self._setup_columns()

        # Convert all rows to formatted strings up-front for fast filtering
        self._all_rows = [
            tuple(_fmt(v) for v in row)
            for row in df.itertuples(index=False, name=None)
        ]
        self._populate_tree(self._all_rows)

    def _setup_columns(self) -> None:
        """Configure Treeview columns and headings."""
        self._tree["columns"] = self._columns
        for col in self._columns:
            self._tree.heading(
                col,
                text=col,
                command=lambda c=col: self._sort_by(c),
            )
            # Estimate a reasonable minimum column width
            width = max(70, min(160, len(col) * 9 + 20))
            self._tree.column(col, width=width, minwidth=50, stretch=True)

    def _populate_tree(self, rows: list[tuple]) -> None:
        """Clear and refill the Treeview with *rows*."""
        self._tree.delete(*self._tree.get_children())
        for row in rows:
            self._tree.insert("", "end", values=row)
        self._update_status()

    def _update_status(self) -> None:
        """Refresh the status bar to reflect visible and selected counts."""
        n_shown = len(self._tree.get_children())
        n_total = len(self._all_rows)
        n_sel = len(self._tree.selection())
        if n_sel:
            self._status_var.set(
                f"{n_sel:,} selected  ·  "
                f"{n_shown:,} visible  ·  {n_total:,} total"
            )
        elif n_shown == n_total:
            self._status_var.set(f"{n_total:,} lines  ·  Ctrl/Cmd+click or Shift+click to select")
        else:
            self._status_var.set(
                f"{n_shown:,} of {n_total:,} lines (filtered)  ·  "
                "Ctrl/Cmd+click or Shift+click to select"
            )

    def _on_selection_change(self, _event=None) -> None:
        """Called whenever the Treeview selection changes."""
        self._update_status()

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def _on_search_change(self, *_) -> None:
        query = self._search_var.get().lower()
        if not query:
            self._populate_tree(self._all_rows)
            return
        filtered = [
            row for row in self._all_rows
            if any(query in cell.lower() for cell in row)
        ]
        self._populate_tree(filtered)

    def _clear_search(self) -> None:
        self._search_var.set("")

    # ------------------------------------------------------------------
    # Column sorting
    # ------------------------------------------------------------------

    def _sort_by(self, col: str) -> None:
        """Sort Treeview rows by the clicked column header."""
        if self._sort_col == col:
            self._sort_reverse = not self._sort_reverse
        else:
            self._sort_col = col
            self._sort_reverse = False

        col_idx = self._columns.index(col)

        # Try numeric sort, fall back to lexicographic
        try:
            current = [
                self._tree.item(iid, "values")
                for iid in self._tree.get_children()
            ]
            sorted_rows = sorted(
                current,
                key=lambda r: float(r[col_idx]),
                reverse=self._sort_reverse,
            )
        except (ValueError, TypeError):
            sorted_rows = sorted(
                [self._tree.item(iid, "values") for iid in self._tree.get_children()],
                key=lambda r: r[col_idx],
                reverse=self._sort_reverse,
            )

        self._tree.delete(*self._tree.get_children())
        for row in sorted_rows:
            self._tree.insert("", "end", values=row)

        # Update heading indicator
        for c in self._columns:
            self._tree.heading(c, text=c)
        arrow = " ▲" if not self._sort_reverse else " ▼"
        self._tree.heading(col, text=col + arrow)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_csv(self) -> None:
        """Save selected rows (or all visible rows if none selected) to a CSV file."""
        selected = self._tree.selection()
        if selected:
            rows = [self._tree.item(iid, "values") for iid in selected]
        else:
            rows = [self._tree.item(iid, "values") for iid in self._tree.get_children()]
        if not rows:
            return

        mol_name = getattr(self.mol_obj, "name", "molecule")
        filepath = filedialog.asksaveasfilename(
            parent=self,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile=f"{mol_name}_line_list",
        )
        if not filepath:
            return

        try:
            import csv
            with open(filepath, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(self._columns)
                writer.writerows(rows)
            msg = f"Line list exported to {filepath}"
            if self.data_field is not None:
                self.data_field.insert_text(msg, console_print=True)
            else:
                print(msg)
        except Exception as exc:
            tk.messagebox.showerror("Export failed", str(exc), parent=self)
