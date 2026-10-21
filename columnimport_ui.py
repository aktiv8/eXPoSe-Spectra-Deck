"""The import dialog for plain column-text files (CSV / ASC / TXT / DAT, SPECS
``.xy``, Excel ``.xlsx``).

Such a file says little about itself, so before it is loaded the user sees what
was found (which column is the energy, which the intensity, binding or kinetic
energy, the unit and the region name), can change any of it with a preview of
the spectrum beside it, and enters what no file of this kind holds: the sample
name, the photon energy and the pass energy. The photon and the pass energy
belong to each file (a header may state them: the SPECS ``.xy`` does); the
boxes under "All files" set them for every file at once, the ones under "This
file" for the selected one. An Excel workbook with several sheets offers the
sheet to read, and a wide CSV written by this program is shown as the
spectra it holds, its names, axes and units taken from its header. The guesses
come from ``columntext`` and are only guesses: nothing is applied that the
dialog does not show.

``items`` are ``(path, columntext.Table, options)``; after ``wait_window`` the
``result`` is ``None`` (cancelled) or ``{path: options}``. Everything the
buttons do is a plain method (``set_row``, ``set_shared``, ``accept``) so the
tests can drive it without a mouse.
"""

from __future__ import annotations

import os
import tkinter as tk
from tkinter import ttk

import columntext
from workbook_ui import _finish

ALL_OTHERS = "All other columns (one spectrum each)"
UNKNOWN = "Unknown"


class ColumnImportDialog(tk.Toplevel):
    def __init__(self, master, app, items, defaults=None):
        super().__init__(master)
        self.app = app
        self.result = None
        d = defaults if isinstance(defaults, dict) else {}
        self.rows = []
        hv0, pe0 = d.get("photon_energy"), d.get("pass_energy")
        for path, table, opts in items:
            o = dict(opts)
            # what the file states stays; the rest starts from what was last
            # entered, and can be changed per file
            if not o.get("photon_energy") and hv0:
                o["photon_energy"] = hv0
                o["anode"] = next((n for n, e in columntext.ANODES.items()
                                   if abs(e - hv0) < 0.05), "")
            if not o.get("pass_energy") and pe0:
                o["pass_energy"] = pe0
            self.rows.append({"path": path, "table": table,
                              "opts": o, "suggested": False})
        self.title("Import spectra from columns of numbers")
        self.transient(master)
        try:
            self.grab_set()
        except tk.TclError:             # not yet viewable (a hidden test root)
            pass
        body = ttk.Frame(self, padding=12)
        body.pack(fill="both", expand=True)
        body.columnconfigure(0, weight=1)
        body.rowconfigure(1, weight=1)
        ttk.Label(
            body, style="Muted.TLabel", wraplength=860, justify="left",
            text="These files hold only numbers, so the program reads them "
                 "the way they look. Check each file: which column is the "
                 "energy, whether it is a binding or a kinetic energy, the "
                 "unit and the name. A name marked ? was suggested from the "
                 "energy range alone; change it if it is wrong. The photon "
                 "energy is needed to turn kinetic energies into binding "
                 "energies; leave it Unknown if you do not know it."
        ).grid(row=0, column=0, sticky="w", pady=(0, 8))

        mid = ttk.Frame(body)
        mid.grid(row=1, column=0, sticky="nsew")
        mid.columnconfigure(0, weight=1)
        cols = ("rows", "energy", "intensity", "axis", "units", "photon",
                "name")
        self.tree = ttk.Treeview(mid, columns=cols, height=9,
                                 selectmode="browse")
        for c, (title, w) in zip(
                ("#0",) + cols, (("File", 210), ("Points", 60),
                                 ("Energy column", 130), ("Intensity", 130),
                                 ("Axis", 100), ("Unit", 70),
                                 ("Photon energy", 110), ("Name", 110))):
            self.tree.heading(c, text=title)
            self.tree.column(c, width=w, stretch=c in ("#0", "name"))
        sb = ttk.Scrollbar(mid, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        sb.grid(row=0, column=1, sticky="ns")
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._on_select())
        self.preview = None
        self._build_preview(mid)

        one = ttk.LabelFrame(body, text="This file", padding=8)
        one.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        self.v_energy = tk.StringVar()
        self.v_inten = tk.StringVar()
        self.v_axis = tk.StringVar()
        self.v_units = tk.StringVar()
        self.v_name = tk.StringVar()
        self.v_sheet = tk.StringVar()
        self.v_hv = tk.StringVar()
        self.v_pe = tk.StringVar()
        for c, (label, var, width) in enumerate((
                ("Energy column", self.v_energy, 22),
                ("Intensity", self.v_inten, 30),
                ("Axis", self.v_axis, 15), ("Unit", self.v_units, 9),
                ("Name", self.v_name, 16))):
            ttk.Label(one, text=label).grid(row=0, column=c, sticky="w",
                                            padx=(0, 6))
            if label == "Name":
                w = ttk.Entry(one, textvariable=var, width=width)
                w.bind("<FocusOut>", lambda e: self._commit())
                w.bind("<Return>", lambda e: self._commit())
                self.name_entry = w
            else:
                w = ttk.Combobox(one, textvariable=var, width=width,
                                 state="readonly")
                w.bind("<<ComboboxSelected>>", lambda e: self._commit())
                setattr(self, "cb_" + label.split()[0].lower(), w)
            w.grid(row=1, column=c, sticky="w", padx=(0, 6))
        self.cb_axis["values"] = [columntext.AXIS_LABELS[a]
                                  for a in columntext.AXES]
        self.cb_unit["values"] = list(columntext.UNITS)
        ttk.Button(one, text="Use this axis and unit for every file",
                   command=self.axis_to_all).grid(row=1, column=5, padx=(8, 0))
        ttk.Label(one, text="Photon energy").grid(row=2, column=0, sticky="w",
                                                  pady=(6, 0))
        self.cb_hv = ttk.Combobox(
            one, textvariable=self.v_hv, width=22,
            values=[UNKNOWN] + [columntext.photon_label(e)
                                for e in columntext.ANODES.values()])
        self.cb_hv.grid(row=3, column=0, sticky="w", padx=(0, 6))
        for ev in ("<<ComboboxSelected>>", "<FocusOut>", "<Return>"):
            self.cb_hv.bind(ev, lambda e: self._commit_energies())
        ttk.Label(one, text="Pass energy (eV)").grid(row=2, column=1,
                                                     sticky="w", pady=(6, 0))
        self.pe_entry = ttk.Entry(one, textvariable=self.v_pe, width=10)
        self.pe_entry.grid(row=3, column=1, sticky="w", padx=(0, 6))
        self.pe_entry.bind("<FocusOut>", lambda e: self._commit_energies())
        self.pe_entry.bind("<Return>", lambda e: self._commit_energies())
        ttk.Label(one, text="Sheet").grid(row=2, column=2, sticky="w",
                                          pady=(6, 0))
        self.cb_sheet = ttk.Combobox(one, textvariable=self.v_sheet, width=18,
                                     state="disabled")
        self.cb_sheet.grid(row=3, column=2, sticky="w", padx=(0, 6))
        self.cb_sheet.bind("<<ComboboxSelected>>",
                           lambda e: self.set_sheet(self.current,
                                                    self.v_sheet.get()))

        every = ttk.LabelFrame(body, text="All files", padding=8)
        every.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        self.v_sample = tk.StringVar(value=d.get("sample", ""))
        hv0 = d.get("photon_energy")
        self.v_photon = tk.StringVar(
            value=columntext.photon_label(hv0) if hv0 else UNKNOWN)
        self.v_pass = tk.StringVar(
            value=f"{d['pass_energy']:g}" if d.get("pass_energy") else "")
        ttk.Label(every, text="Sample name").grid(row=0, column=0, sticky="w")
        ttk.Entry(every, textvariable=self.v_sample, width=28).grid(
            row=1, column=0, sticky="w", padx=(0, 10))
        ttk.Label(every, text="Photon energy (every file)").grid(
            row=0, column=1, sticky="w")
        self.cb_photon = ttk.Combobox(
            every, textvariable=self.v_photon, width=22,
            values=[UNKNOWN] + [columntext.photon_label(e)
                                for e in columntext.ANODES.values()])
        self.cb_photon.grid(row=1, column=1, sticky="w", padx=(0, 10))
        for ev in ("<<ComboboxSelected>>", "<Return>"):
            self.cb_photon.bind(ev, lambda e: self._apply_all())
        ttk.Label(every, text="Pass energy (eV, every file)").grid(
            row=0, column=2, sticky="w")
        pass_entry = ttk.Entry(every, textvariable=self.v_pass, width=10)
        pass_entry.grid(row=1, column=2, sticky="w")
        pass_entry.bind("<Return>", lambda e: self._apply_all())
        ttk.Button(every, text="Apply to every file",
                   command=self._apply_all).grid(row=1, column=3, padx=(10, 0))

        self.notes = ttk.Label(body, style="Muted.TLabel", wraplength=860,
                               justify="left")
        self.notes.grid(row=4, column=0, sticky="w", pady=(8, 0))
        btns = ttk.Frame(body)
        btns.grid(row=5, column=0, sticky="e", pady=(10, 0))
        ttk.Button(btns, text="Import", command=self.accept).pack(
            side="left", padx=(0, 8))
        ttk.Button(btns, text="Cancel", command=self.cancel).pack(side="left")
        self.bind("<Escape>", lambda e: self.cancel())
        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self._suggest_names()
        for i in range(len(self.rows)):
            self.tree.insert("", "end", iid=str(i))
            self._refresh_row(i)
        self._write_notes()
        if self.rows:
            self.tree.selection_set("0")
            self._on_select()
        _finish(self, app, 1200, 660)

    # -- the rows -----------------------------------------------------------------
    def _suggest_names(self):
        """Offer a core level for each file the header did not name (only when
        one common element's strongest line is in the window)."""
        hv = columntext.parse_photon(self.v_photon.get())
        for row in self.rows:
            o, t = row["opts"], row["table"]
            if t.pairs or o["axis"] != "BE" or len(o["intensity_cols"]) != 1:
                continue
            if columntext.region_name(o["names"][0],
                                      t.columns[o["energy_col"]]):
                continue
            e = t.columns[o["energy_col"]]
            name = columntext.suggest_name(min(e), max(e), hv)
            if name:
                o["names"] = [name]
                row["suggested"] = True

    def _name_text(self, row):
        o = row["opts"]
        if row["table"].pairs:
            return f"{len(row['table'].pairs)} spectra"
        if len(o["intensity_cols"]) > 1:
            return f"{len(o['intensity_cols'])} spectra"
        name = o["names"][0] if o["names"] else ""
        shown = columntext.region_name(
            name, row["table"].columns[o["energy_col"]])
        if not shown:
            return "(file name)"
        return shown + (" ?" if row["suggested"] else "")

    def _refresh_row(self, i):
        row = self.rows[i]
        o, t = row["opts"], row["table"]
        inten = ", ".join(t.label(c) for c in o["intensity_cols"])
        axis, units = columntext.AXIS_LABELS[o["axis"]], o["units"]
        if t.pairs:
            inten, axis, units = "(from the file)", "(from the file)", "…"
        text = os.path.basename(row["path"])
        if t.sheet:
            text += f" [{t.sheet}]"
        self.tree.item(str(i), text=text, values=(
            t.n_rows, t.label(o["energy_col"]), inten, axis, units,
            columntext.photon_label(o["photon_energy"])
            if o["photon_energy"] else UNKNOWN, self._name_text(row)))

    def _write_notes(self):
        lines = []
        for row in self.rows:
            t = row["table"]
            if t.trailing:
                lines.append(f"{os.path.basename(row['path'])}: {t.trailing} "
                             f"line(s) after the table (from line "
                             f"{t.trailing_at}) will be ignored.")
            if not columntext.find_energy_col(t)[1]:
                lines.append(f"{os.path.basename(row['path'])}: no column "
                             f"runs steadily one way, so the energy column is "
                             f"a guess.")
            if columntext.is_specs(t.facts):
                lines.append(f"{os.path.basename(row['path'])}: a SPECS file; "
                             "its header gives the excitation and pass energy "
                             "(change them here if they are wrong).")
            others = [n for n, r in t.sheets if n != t.sheet and r]
            if others:
                lines.append(f"{os.path.basename(row['path'])}: other sheets "
                             f"with a table: {', '.join(others)} (choose one "
                             "under Sheet).")
        self.notes.config(text="\n".join(lines[:6])
                          + ("\n…" if len(lines) > 6 else ""))

    # -- editing one file ------------------------------------------------------------
    @property
    def current(self):
        sel = self.tree.selection()
        return int(sel[0]) if sel else None

    def _col_choices(self, t):
        return [f"{c + 1}: {t.label(c)}" for c in range(t.n_cols)]

    def _on_select(self):
        i = self.current
        if i is None:
            return
        row = self.rows[i]
        o, t = row["opts"], row["table"]
        choices = self._col_choices(t)
        self.cb_energy["values"] = choices
        inten = list(choices)
        if t.n_cols > 2:
            inten.append(ALL_OTHERS)
        self.cb_intensity["values"] = inten
        self.v_energy.set(choices[o["energy_col"]])
        others = [c for c in range(t.n_cols) if c != o["energy_col"]]
        if len(o["intensity_cols"]) > 1 and o["intensity_cols"] == others:
            self.v_inten.set(ALL_OTHERS)
        else:
            self.v_inten.set(choices[o["intensity_cols"][0]])
        self.v_axis.set(columntext.AXIS_LABELS[o["axis"]])
        self.v_units.set(o["units"])
        multi = len(o["intensity_cols"]) > 1 or bool(t.pairs)
        self.v_name.set("" if multi else (o["names"][0] if o["names"] else ""))
        self.name_entry.config(state="disabled" if multi else "normal")
        wide = "disabled" if t.pairs else "readonly"
        for cb in (self.cb_energy, self.cb_intensity, self.cb_axis,
                   self.cb_unit):
            cb.config(state=wide)
        self.v_hv.set(columntext.photon_label(o["photon_energy"])
                      if o["photon_energy"] else UNKNOWN)
        self.v_pe.set(f"{o['pass_energy']:g}" if o["pass_energy"] else "")
        names = [n for n, _r in t.sheets]
        self.cb_sheet["values"] = names
        self.v_sheet.set(t.sheet)
        self.cb_sheet.config(state="readonly" if len(names) > 1 else "disabled")
        self._draw_preview()

    def _commit_energies(self):
        i = self.current
        if i is not None:
            self.set_row(i, photon=self.v_hv.get(), pass_energy=self.v_pe.get())

    def _apply_all(self):
        """The photon and pass energy typed under "All files" go to every
        file."""
        for k in range(len(self.rows)):
            self.set_row(k, photon=self.v_photon.get(),
                         pass_energy=self.v_pass.get())

    def set_sheet(self, i, name):
        """Read sheet ``name`` of an Excel file instead (the table, its guesses
        and so its options change; the energies and the sample stay)."""
        if i is None:
            return
        row = self.rows[i]
        if not name or name == row["table"].sheet:
            return
        try:
            table = columntext.load_table(row["path"], name)
        except (OSError, ValueError):
            return
        if table is None:
            return
        keep = {k: row["opts"][k] for k in ("photon_energy", "pass_energy",
                                            "anode", "sample")}
        row["table"] = table
        row["opts"] = dict(columntext.guess_options(table, row["path"]),
                           **{k: v for k, v in keep.items() if v})
        row["suggested"] = False
        self._suggest_names()
        self._refresh_row(i)
        if self.current == i:
            self._on_select()
        self._write_notes()

    def _commit(self):
        """Apply the boxes under the table to the selected file."""
        i = self.current
        if i is None:
            return
        t = self.rows[i]["table"]
        if t.pairs:                         # its header settles all of this
            return
        choices = self._col_choices(t)

        def col(text):
            return choices.index(text) if text in choices else None

        ecol = col(self.v_energy.get())
        inten = self.v_inten.get()
        axis = next((a for a, lab in columntext.AXIS_LABELS.items()
                     if lab == self.v_axis.get()), None)
        self.set_row(i, energy_col=ecol,
                     intensity="all" if inten == ALL_OTHERS else col(inten),
                     axis=axis, units=self.v_units.get(),
                     name=self.v_name.get())

    def set_row(self, i, energy_col=None, intensity=None, axis=None,
                units=None, name=None, photon=None, pass_energy=None):
        """Change what file ``i`` is read with. ``intensity`` is a column
        number, a list of them or ``"all"`` (every column but the energy);
        the energy and intensity columns never coincide. ``photon`` (text:
        "Al Kα (1486.6 eV)", "1253.6", "Unknown") and ``pass_energy`` (text)
        are the file's own."""
        row = self.rows[i]
        o, t = dict(row["opts"]), row["table"]
        if t.pairs:                         # columns, axis, unit, name: its own
            energy_col = intensity = axis = units = name = None
        if photon is not None:
            o["photon_energy"] = columntext.parse_photon(photon)
            o["anode"] = next((n for n, e in columntext.ANODES.items()
                               if o["photon_energy"]
                               and abs(e - o["photon_energy"]) < 0.05), "")
        if pass_energy is not None:
            o["pass_energy"] = columntext.parse_photon(pass_energy)
        if energy_col is not None or intensity is not None:
            old_cols, old_names = o["intensity_cols"], o["names"]
            if energy_col is not None:
                o["energy_col"] = energy_col
            others = [c for c in range(t.n_cols) if c != o["energy_col"]]
            if intensity == "all":
                cols = others
            elif isinstance(intensity, int):
                cols = [intensity]
            elif isinstance(intensity, list):
                cols = list(intensity)
            else:
                cols = [c for c in old_cols if c != o["energy_col"]] \
                    or others[:1]
            cols = [c for c in cols if c != o["energy_col"]] or others[:1]
            if len(cols) > 1:
                names = [t.header[c].strip() if t.header else ""
                         for c in cols]
            else:
                keep = cols == old_cols and len(old_names) == 1
                names = [old_names[0] if keep else ""]
                if not keep:
                    row["suggested"] = False
            o["intensity_cols"], o["names"] = cols, names
        if axis in columntext.AXES:
            o["axis"] = axis
        if units in columntext.UNITS:
            o["units"] = units
        if name is not None and len(o["intensity_cols"]) == 1:
            new = name.strip()
            if new != (o["names"][0] if o["names"] else ""):
                o["names"] = [new]
                row["suggested"] = False
        row["opts"] = columntext.sanitise_options(o, t.n_cols)
        self._refresh_row(i)
        if self.current == i:
            self._on_select()

    def axis_to_all(self):
        i = self.current
        if i is None:
            return
        o = self.rows[i]["opts"]
        for k in range(len(self.rows)):
            if k != i:
                self.set_row(k, axis=o["axis"], units=o["units"])

    # -- the settings every file shares ----------------------------------------------
    def set_shared(self, sample=None, photon=None, pass_energy=None):
        """Set the boxes under "All files"; the photon and the pass energy go
        to every file at once (each can be changed again for one file)."""
        if sample is not None:
            self.v_sample.set(sample)
        if photon is not None:
            self.v_photon.set(photon)
        if pass_energy is not None:
            self.v_pass.set(str(pass_energy))
        if photon is not None or pass_energy is not None:
            self._apply_all()

    def shared(self):
        hv = columntext.parse_photon(self.v_photon.get())
        anode = next((n for n, e in columntext.ANODES.items()
                      if hv and abs(e - hv) < 0.05), "")
        pe = columntext.parse_photon(self.v_pass.get())
        return {"sample": self.v_sample.get().strip(), "photon_energy": hv,
                "anode": anode, "pass_energy": pe}

    # -- preview ----------------------------------------------------------------------
    def _build_preview(self, parent):
        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from matplotlib.figure import Figure
        except Exception:                        # noqa: BLE001 - optional
            return
        pal = self.app.palette
        self.fig = Figure(figsize=(3.4, 2.4), dpi=90, facecolor=pal["bg"])
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().grid(row=0, column=2, padx=(10, 0),
                                         sticky="n")
        self.preview = self.canvas

    def _draw_preview(self):
        if self.preview is None or self.current is None:
            return
        row = self.rows[self.current]
        o, t = row["opts"], row["table"]
        pal = self.app.palette
        ax = self.ax
        ax.clear()
        ax.set_facecolor(pal["plot_bg"])
        if t.pairs:
            for pr in t.pairs[:6]:
                ax.plot(t.columns[pr["energy"]], t.columns[pr["intensity"]],
                        lw=0.8)
        else:
            x = t.columns[o["energy_col"]]
            for c in o["intensity_cols"][:6]:
                ax.plot(x, t.columns[c], lw=0.8)
        ax.tick_params(colors=pal["muted"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(pal["muted"])
        ax.set_xlabel(f"{columntext.AXIS_LABELS[o['axis']]} (eV)",
                      color=pal["plot_fg"], fontsize=8)
        if o["axis"] == "BE" and not t.pairs:
            ax.invert_xaxis()
        self.fig.subplots_adjust(left=0.2, right=0.96, top=0.95, bottom=0.2)
        self.canvas.draw_idle()

    # -- finishing ----------------------------------------------------------------------
    def accept(self):
        self._commit()
        self._commit_energies()
        sample = self.v_sample.get().strip()
        out = {}
        for row in self.rows:
            o = dict(row["opts"])
            if sample:                      # one sample name for every file
                o["sample"] = sample
            out[row["path"]] = columntext.sanitise_options(
                o, row["table"].n_cols)
        self.result = out
        self.destroy()

    def cancel(self):
        self.result = None
        self.destroy()
