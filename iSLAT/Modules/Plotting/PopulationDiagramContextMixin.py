"""
PopulationDiagramContextMixin — shared right-click context menu for the
population (Boltzmann) diagram panel.

Both :class:`ThreePanelView` and :class:`PopulationDiagramView` inherit
this mixin so that the *Color By* and *Axis Settings* dialogs, and the
menu-building logic, live in exactly one place.

Usage::

    class MyView(PlotView, PopulationDiagramContextMixin):
        def build_context_menu(self, event, canvas_widget):
            pdp = ...  # get the PopulationDiagramPlot instance
            draw_idle = ...  # callable that flushes the canvas
            return self._build_population_diagram_menu(pdp, canvas_widget, draw_idle)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Optional

class PopulationDiagramContextMixin:
    """Mixin that supplies the population-diagram right-click menu and dialogs.

    Methods intentionally receive ``pdp`` (a :class:`PopulationDiagramPlot`)
    and ``draw_idle`` (a zero-argument callable that flushes the canvas) as
    explicit parameters so the mixin stays decoupled from any particular view
    structure.
    """

    # ------------------------------------------------------------------
    # Public helper — build the tk.Menu
    # ------------------------------------------------------------------

    def _build_population_diagram_menu(
        self,
        pdp: Any,
        canvas_widget: Any,
        draw_idle: Callable[[], None],
    ) -> Optional[Any]:
        """Return a populated ``tk.Menu`` for the population diagram, or ``None``.

        Parameters
        ----------
        pdp :
            The active :class:`PopulationDiagramPlot` instance, or ``None``
            if no diagram has been rendered yet.
        canvas_widget :
            The Tk widget that acts as the menu's parent.
        draw_idle :
            Zero-argument callable that triggers a canvas redraw
            (typically ``canvas.draw_idle``).
        """
        try:
            import tkinter as tk
        except ImportError:
            return None

        menu = tk.Menu(canvas_widget, tearoff=0)

        # ---- helpers --------------------------------------------------

        def _color_by_dialog():
            self._open_color_by_dialog(pdp, canvas_widget, draw_idle)

        def _clear_color_mapping():
            if pdp is not None:
                pdp.clear_color_mapping(regenerate=True)
                draw_idle()

        def _toggle_log_scale():
            if pdp is None or getattr(pdp, '_color_mapping', None) is None:
                return
            current = pdp._color_mapping.get("log_scale", False)
            pdp._color_mapping["log_scale"] = not current
            pdp.generate_plot()
            draw_idle()

        def _show_all_active_molecules():
            if pdp is None:
                return
            islat = getattr(self, '_islat', None)
            molecules_dict = getattr(islat, 'molecules_dict', None) if islat else None
            if molecules_dict is None:
                return
            pdp.set_molecules(molecules_dict)
            draw_idle()

        def _show_active_molecule_only():
            if pdp is None:
                return
            islat = getattr(self, '_islat', None)
            active_mol = getattr(islat, 'active_molecule', None) if islat else None
            if active_mol is None:
                return
            if isinstance(active_mol, str):
                molecules_dict = getattr(islat, 'molecules_dict', None)
                if molecules_dict is not None and active_mol in molecules_dict:
                    active_mol = molecules_dict[active_mol]
                else:
                    return
            pdp.set_molecule(active_mol)
            draw_idle()

        # ---- state indicators -----------------------------------------

        _all_mode = pdp is not None and getattr(pdp, '_all_molecules_mode', False)
        _has_mapping = pdp is not None and getattr(pdp, '_color_mapping', None) is not None
        _log_scale = _has_mapping and pdp._color_mapping.get("log_scale", False)

        # ---- menu items -----------------------------------------------

        menu.add_command(label="Color By\u2026", command=_color_by_dialog)
        menu.add_command(label="Clear Color Mapping", command=_clear_color_mapping)
        menu.add_command(
            label=("\u2713 " if _log_scale else "  ") + "Log Scale Colorbar",
            command=_toggle_log_scale,
            state="normal" if _has_mapping else "disabled",
        )
        menu.add_separator()
        menu.add_command(
            label="Axis Settings\u2026",
            command=lambda: self._open_axis_settings_dialog(pdp, canvas_widget, draw_idle),
        )
        menu.add_separator()
        menu.add_command(
            label=("\u2713 " if _all_mode else "  ") + "Show All Active Molecules",
            command=_show_all_active_molecules,
        )
        menu.add_command(
            label=("\u2713 " if not _all_mode else "  ") + "Show Active Molecule Only",
            command=_show_active_molecule_only,
        )
        menu.add_separator()

        # --- Save active lines as line list ----------------------------
        _active_info = self._get_active_lines_info(pdp)

        def _save_line_list():
            self._save_active_lines_as_line_list(canvas_widget)

        menu.add_command(
            label="Save Active Lines as Line List\u2026",
            command=_save_line_list,
            state="normal" if _active_info else "disabled",
        )

        # ------------------------------------------------------------------
        # View switching
        # ------------------------------------------------------------------
        pm = getattr(self, '_pm', None)
        if pm is not None and hasattr(pm, 'switch_view'):
            menu.add_separator()

            active_name = getattr(pm, 'active_view_name', None)

            if active_name == "Three Panel":
                # Inside the three-panel layout — offer to pop out to standalone views
                def _to_population_diagram():
                    pm.switch_view("Population Diagram")

                def _to_line_inspection():
                    pm.switch_view("Line Inspection")

                menu.add_command(
                    label="Switch to Population Diagram View",
                    command=_to_population_diagram,
                )
                menu.add_command(
                    label="Switch to Line Inspection View",
                    command=_to_line_inspection,
                )
            else:
                # Already in a standalone view — offer to go back or cross-navigate
                def _to_three_panel():
                    pm.switch_view("Three Panel")

                def _to_line_inspection():
                    pm.switch_view("Line Inspection")

                menu.add_command(
                    label="Switch to Three Panel View",
                    command=_to_three_panel,
                )
                menu.add_command(
                    label="Switch to Line Inspection View",
                    command=_to_line_inspection,
                )

        return menu

    # ------------------------------------------------------------------
    # Save active lines as line list
    # ------------------------------------------------------------------

    def _get_active_lines_info(self, pdp: Any) -> list:
        """Return the list of ``info_dict`` entries for currently active lines.

        Tries, in order:
        1. ``self.active_lines`` — populated by :class:`ThreePanelView`.
        2. ``pdp._active_lines_cache`` — populated by the standalone
           :class:`PopulationDiagramView`.

        Returns an empty list when no active lines are available.
        """
        # ThreePanelView path
        own_lines = getattr(self, "active_lines", None)
        if own_lines:
            return [entry[3] for entry in own_lines if entry[3] is not None]

        # Standalone PopulationDiagramView path
        if pdp is not None:
            cache = getattr(pdp, "_active_lines_cache", None)
            if cache:
                cached_entries = cache.get("active_lines", [])
                if cached_entries:
                    return [entry[3] for entry in cached_entries if entry[3] is not None]

        return []

    def _save_active_lines_as_line_list(self, parent_widget: Any) -> None:
        """Open a save-as dialog and write the currently active lines to a
        CSV line list file compatible with iSLAT's fitting workflow.

        Each row contains the columns produced by
        :meth:`LineSaveService.format_line_for_save`:
        ``species``, ``lam``, ``xmin``, ``xmax``, ``lev_up``,
        ``lev_low``, ``tau``, ``intens``, ``a_stein``,
        ``e_up``, ``g_up``, ``e_low``, ``g_low``.
        """
        try:
            import tkinter as tk
            from tkinter import filedialog, messagebox
            import pandas as pd
        except ImportError:
            return

        # Resolve the active pdp so we can read its cache
        pdp = getattr(self, "_plot", None)
        info_list = self._get_active_lines_info(pdp)

        if not info_list:
            messagebox.showinfo(
                "No Active Lines",
                "No highlighted lines are available to save.\n"
                "Select a wavelength range first to activate lines.",
                parent=parent_widget,
            )
            return

        # Build rows — derive xmin/xmax as a ±0.015 µm window if not present
        _DEFAULT_HALF_WIN = 0.015
        rows = []
        for info in info_list:
            lam = info.get("lam")
            if lam is None:
                continue
            lam = float(lam)
            rows.append({
                "species":  info.get("molecule_name", ""),
                "lam":      lam,
                "xmin":     float(info.get("xmin", lam - _DEFAULT_HALF_WIN)),
                "xmax":     float(info.get("xmax", lam + _DEFAULT_HALF_WIN)),
                "lev_up":   info.get("up_lev", ""),
                "lev_low":  info.get("low_lev", ""),
                "tau":      info.get("tau", 0.0),
                "intens":   info.get("intensity", info.get("inten", 0.0)),
                "a_stein":  info.get("a_stein", 0.0),
                "e_up":     info.get("e_up", 0.0),
                "g_up":     info.get("g_up", 1.0),
                "e_low":    info.get("e_low", 0.0),
                "g_low":    info.get("g_low", 1.0),
            })

        if not rows:
            messagebox.showinfo(
                "No Valid Lines",
                "None of the active lines contained wavelength data.",
                parent=parent_widget,
            )
            return

        df = pd.DataFrame(rows)

        filepath = filedialog.asksaveasfilename(
            title="Save Active Lines as Line List",
            defaultextension=".csv",
            filetypes=[
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
            initialfile="active_lines.csv",
            parent=parent_widget,
        )
        if not filepath:
            return

        try:
            df.to_csv(filepath, index=False)
            messagebox.showinfo(
                "Saved",
                f"Saved {len(df)} line(s) to:\n{filepath}",
                parent=parent_widget,
            )
        except Exception as exc:
            messagebox.showerror(
                "Save Failed",
                f"Could not write file:\n{exc}",
                parent=parent_widget,
            )

    # ------------------------------------------------------------------
    # Color By dialog
    # ------------------------------------------------------------------

    def _open_color_by_dialog(
        self,
        pdp: Any,
        parent_widget: Any,
        draw_idle: Callable[[], None],
    ) -> None:
        """Open the *Color By* dialog for *pdp*."""
        try:
            import tkinter as tk
            import tkinter.ttk as ttk
        except ImportError:
            return

        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine as _ML

        _EXTRA = [
            ("Population diagram Y [ln(4\u03c0F/h\u03bdA_u g_u)]", "rd_yax"),
            ("Model intensity",                                      "intens"),
            ("Line-center opacity (tau)",                            "tau"),
            ("Molecule / component",                                 "molecule"),
        ]
        _FROM_REG = [
            (_ML.get_text(k), k if k != "lam" else "wavelength")
            for k in ("e_up", "e_low", "a_stein", "g_up", "g_low",
                      "lev_up", "lev_low", "freq", "lam")
        ]
        # Quantum field options — populated from the molecule's schema when
        # pdp has component_data with a recognised molecule_id.
        _QUANTUM: list = []
        if pdp is not None:
            try:
                import iSLAT.Modules.DataTypes.HITRANQuantumSchemas  # noqa: F401 (registers all schemas)
                from iSLAT.Modules.DataTypes.QuantumStateSchema import QuantumStateRegistry
                _comp_data = getattr(pdp, 'component_data', [])
                _mol_ids = [c.get('molecule_id') for c in _comp_data if c.get('molecule_id')]
                if _mol_ids:
                    _schema = QuantumStateRegistry.get_schema(_mol_ids[0])
                    for _f in list(_schema.global_fields) + list(_schema.local_fields):
                        _QUANTUM.append((f"{_f.name} — upper state (quantum)", f"qn_upper:{_f.name}"))
                        _QUANTUM.append((f"{_f.name} — lower state (quantum)", f"qn_lower:{_f.name}"))
            except Exception:
                pass

        PROP_OPTIONS = _EXTRA + _FROM_REG + _QUANTUM
        PROP_LABELS  = [p[0] for p in PROP_OPTIONS]
        PROP_VALUES  = [p[1] for p in PROP_OPTIONS]

        CMAP_OPTIONS = [
            "viridis", "plasma", "inferno", "magma", "cividis",
            "coolwarm", "RdYlBu", "Spectral", "tab10", "tab20",
        ]

        current_mapping = getattr(pdp, '_color_mapping', None) if pdp else None
        initial_prop_idx = 0
        initial_cmap = "viridis"
        initial_vmin = initial_vmax = initial_pmin = initial_pmax = ""
        initial_log_scale = False

        # Properties that naturally span orders of magnitude benefit from
        # a log color scale — pre-tick the checkbox for these.
        _LOG_SUGGESTED = {"intens", "tau", "a_stein"}

        if current_mapping:
            prop_val = current_mapping.get("prop", "e_up")
            if prop_val in PROP_VALUES:
                initial_prop_idx = PROP_VALUES.index(prop_val)
            initial_cmap = current_mapping.get("cmap", "viridis")
            v = current_mapping.get("vmin")
            initial_vmin = str(v) if v is not None else ""
            v = current_mapping.get("vmax")
            initial_vmax = str(v) if v is not None else ""
            v = current_mapping.get("pmin")
            initial_pmin = str(v) if v is not None else ""
            v = current_mapping.get("pmax")
            initial_pmax = str(v) if v is not None else ""
            initial_log_scale = bool(current_mapping.get("log_scale", False))
        else:
            # No existing mapping: suggest log scale for log-distributed props
            initial_prop = PROP_VALUES[initial_prop_idx]
            initial_log_scale = initial_prop in _LOG_SUGGESTED

        win = tk.Toplevel(parent_widget)
        win.title("Population Diagram \u2014 Color By")
        win.resizable(False, False)
        win.grab_set()

        pad = {"padx": 8, "pady": 4}

        tk.Label(win, text="Color by property:").grid(row=0, column=0, sticky="w", **pad)
        prop_var = tk.StringVar(value=PROP_LABELS[initial_prop_idx])
        log_scale_var = tk.BooleanVar(value=initial_log_scale)

        def _on_prop_changed(*_):
            lbl = prop_var.get()
            if lbl in PROP_LABELS:
                pval = PROP_VALUES[PROP_LABELS.index(lbl)]
                if pval in _LOG_SUGGESTED:
                    log_scale_var.set(True)

        prop_combo = ttk.Combobox(win, textvariable=prop_var, values=PROP_LABELS,
                     state="readonly", width=32)
        prop_combo.grid(row=0, column=1, sticky="ew", **pad)
        prop_combo.bind("<<ComboboxSelected>>", _on_prop_changed)

        tk.Label(win, text="Colormap:").grid(row=1, column=0, sticky="w", **pad)
        cmap_var = tk.StringVar(value=initial_cmap)
        ttk.Combobox(win, textvariable=cmap_var, values=CMAP_OPTIONS,
                     state="normal", width=20).grid(row=1, column=1, sticky="ew", **pad)

        tk.Label(win, text="vmin (optional):").grid(row=2, column=0, sticky="w", **pad)
        vmin_var = tk.StringVar(value=initial_vmin)
        tk.Entry(win, textvariable=vmin_var, width=12).grid(row=2, column=1, sticky="w", **pad)

        tk.Label(win, text="vmax (optional):").grid(row=3, column=0, sticky="w", **pad)
        vmax_var = tk.StringVar(value=initial_vmax)
        tk.Entry(win, textvariable=vmax_var, width=12).grid(row=3, column=1, sticky="w", **pad)

        ttk.Separator(win, orient="horizontal").grid(
            row=4, column=0, columnspan=2, sticky="ew", padx=8, pady=2)
        tk.Label(win, text="\u2014 or clip by percentile (overrides vmin/vmax) \u2014",
                 font=("TkDefaultFont", 8), foreground="gray"
                 ).grid(row=5, column=0, columnspan=2, **pad)

        tk.Label(win, text="Min percentile % (optional):").grid(row=6, column=0, sticky="w", **pad)
        pmin_var = tk.StringVar(value=initial_pmin)
        tk.Spinbox(win, from_=0, to=100, increment=1, textvariable=pmin_var,
                   width=8).grid(row=6, column=1, sticky="w", **pad)

        tk.Label(win, text="Max percentile % (optional):").grid(row=7, column=0, sticky="w", **pad)
        pmax_var = tk.StringVar(value=initial_pmax)
        tk.Spinbox(win, from_=0, to=100, increment=1, textvariable=pmax_var,
                   width=8).grid(row=7, column=1, sticky="w", **pad)

        tk.Checkbutton(win, text="Log scale colorbar",
                       variable=log_scale_var).grid(row=8, column=0, columnspan=2,
                                                    sticky="w", **pad)

        btn_frame = tk.Frame(win)
        btn_frame.grid(row=9, column=0, columnspan=2, pady=6)

        def _apply():
            prop_label = prop_var.get()
            if prop_label not in PROP_LABELS:
                return
            prop = PROP_VALUES[PROP_LABELS.index(prop_label)]
            cmap = cmap_var.get().strip() or "viridis"
            try:
                vmin = float(vmin_var.get()) if vmin_var.get().strip() else None
            except ValueError:
                vmin = None
            try:
                vmax = float(vmax_var.get()) if vmax_var.get().strip() else None
            except ValueError:
                vmax = None
            try:
                pmin = float(pmin_var.get()) if pmin_var.get().strip() else None
                if pmin is not None:
                    pmin = max(0.0, min(100.0, pmin))
            except ValueError:
                pmin = None
            try:
                pmax = float(pmax_var.get()) if pmax_var.get().strip() else None
                if pmax is not None:
                    pmax = max(0.0, min(100.0, pmax))
            except ValueError:
                pmax = None

            if pdp is not None:
                pdp.color_by(prop, cmap=cmap, vmin=vmin, vmax=vmax,
                             pmin=pmin, pmax=pmax,
                             log_scale=log_scale_var.get(), regenerate=True)
                draw_idle()
            win.destroy()

        def _cancel():
            win.destroy()

        ttk.Button(btn_frame, text="Apply",  command=_apply).pack(side="left",  padx=4)
        ttk.Button(btn_frame, text="Cancel", command=_cancel).pack(side="left", padx=4)

        win.bind("<Return>", lambda _e: _apply())
        win.bind("<Escape>", lambda _e: _cancel())

    # ------------------------------------------------------------------
    # Axis Settings dialog
    # ------------------------------------------------------------------

    def _open_axis_settings_dialog(
        self,
        pdp: Any,
        parent_widget: Any,
        draw_idle: Callable[[], None],
    ) -> None:
        """Open the *Axis Settings* dialog for *pdp*."""
        try:
            import tkinter as tk
            import tkinter.ttk as ttk
        except ImportError:
            return

        from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine as _ML

        _AXIS_EXTRA = [
            ("Population diagram Y [ln(4\u03c0F/h\u03bdA_u g_u)]", "rd_yax"),
            ("Model intensity",                                      "intens"),
            ("Line-center opacity (tau)",                            "tau"),
            ("Line instrumental FWHM (km/s)",                       "fwhm_instrumental_kms"),
            ("Line convolved FWHM (km/s)",                          "fwhm_convolved_kms"),
        ]
        _FROM_REG_AXIS = [
            (_ML.get_text(k), k if k != "lam" else "wavelength")
            for k in ("e_up", "e_low", "lam", "a_stein", "g_up", "g_low")
        ]
        # Quantum field options — populated from the molecule's schema when
        # pdp has component_data with a recognised molecule_id.
        _QUANTUM_AXIS: list = []
        if pdp is not None:
            try:
                import iSLAT.Modules.DataTypes.HITRANQuantumSchemas  # noqa: F401
                from iSLAT.Modules.DataTypes.QuantumStateSchema import QuantumStateRegistry
                _comp_data = getattr(pdp, 'component_data', [])
                _mol_ids = [c.get('molecule_id') for c in _comp_data if c.get('molecule_id')]
                if _mol_ids:
                    _schema = QuantumStateRegistry.get_schema(_mol_ids[0])
                    for _f in list(_schema.global_fields) + list(_schema.local_fields):
                        _QUANTUM_AXIS.append((f"{_f.name} — upper state (quantum)", f"qn_upper:{_f.name}"))
                        _QUANTUM_AXIS.append((f"{_f.name} — lower state (quantum)", f"qn_lower:{_f.name}"))
            except Exception:
                pass

        AXIS_OPTIONS = [
            (label, "eu" if key == "e_up" else key)
            for label, key in (_AXIS_EXTRA + _FROM_REG_AXIS)
        ] + _QUANTUM_AXIS
        AXIS_LABELS = [o[0] for o in AXIS_OPTIONS]
        AXIS_VALUES = [o[1] for o in AXIS_OPTIONS]

        cur_x     = getattr(pdp, '_x_prop', 'eu')     if pdp else 'eu'
        cur_y     = getattr(pdp, '_y_prop', 'rd_yax') if pdp else 'rd_yax'
        cur_x_log = getattr(pdp, '_x_log', False)     if pdp else False
        cur_y_log = getattr(pdp, '_y_log', False)     if pdp else False

        init_x_idx = AXIS_VALUES.index(cur_x) if cur_x in AXIS_VALUES else 1
        init_y_idx = AXIS_VALUES.index(cur_y) if cur_y in AXIS_VALUES else 0

        win = tk.Toplevel(parent_widget)
        win.title("Population Diagram \u2014 Axis Settings")
        win.resizable(False, False)
        win.grab_set()

        pad = {"padx": 8, "pady": 4}

        tk.Label(win, text="X axis:").grid(row=0, column=0, sticky="w", **pad)
        x_var = tk.StringVar(value=AXIS_LABELS[init_x_idx])
        ttk.Combobox(win, textvariable=x_var, values=AXIS_LABELS,
                     state="readonly", width=38).grid(row=0, column=1, sticky="ew", **pad)

        x_log_var = tk.BooleanVar(value=cur_x_log)
        tk.Checkbutton(win, text="Log scale (X)",
                       variable=x_log_var).grid(row=1, column=0, columnspan=2,
                                                sticky="w", **pad)

        ttk.Separator(win, orient="horizontal").grid(
            row=2, column=0, columnspan=2, sticky="ew", padx=8, pady=2)

        tk.Label(win, text="Y axis:").grid(row=3, column=0, sticky="w", **pad)
        y_var = tk.StringVar(value=AXIS_LABELS[init_y_idx])
        ttk.Combobox(win, textvariable=y_var, values=AXIS_LABELS,
                     state="readonly", width=38).grid(row=3, column=1, sticky="ew", **pad)

        y_log_var = tk.BooleanVar(value=cur_y_log)
        tk.Checkbutton(win, text="Log scale (Y)",
                       variable=y_log_var).grid(row=4, column=0, columnspan=2,
                                                sticky="w", **pad)

        btn_frame = tk.Frame(win)
        btn_frame.grid(row=5, column=0, columnspan=2, pady=6)

        def _apply():
            x_label = x_var.get()
            y_label = y_var.get()
            if x_label not in AXIS_LABELS or y_label not in AXIS_LABELS:
                return
            x_prop = AXIS_VALUES[AXIS_LABELS.index(x_label)]
            y_prop = AXIS_VALUES[AXIS_LABELS.index(y_label)]
            if pdp is not None:
                pdp.set_axes(
                    x_prop=x_prop, y_prop=y_prop,
                    x_log=x_log_var.get(), y_log=y_log_var.get(),
                    regenerate=True,
                )
                draw_idle()
            win.destroy()

        def _reset():
            if pdp is not None:
                pdp.set_axes(regenerate=True)
                draw_idle()
            win.destroy()

        def _cancel():
            win.destroy()

        ttk.Button(btn_frame, text="Apply",         command=_apply).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="Reset Defaults", command=_reset).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="Cancel",        command=_cancel).pack(side="left", padx=4)

        win.bind("<Return>", lambda _e: _apply())
        win.bind("<Escape>", lambda _e: _cancel())