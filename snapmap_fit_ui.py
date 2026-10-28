"""The *Fit* tab of the SnapMap window (a mixin of ``snapmap_ui.SnapMapDialog``).

Components of a core level (Si and SiO2 in Si 2p, or whatever the level has),
seeded from the chemical-state table or the line table, are refined on the mean
spectrum of the whole map or of an area; the maps are then drawn from the frozen
shapes (``snapfit``, ``snapmapfig.draw_fit_view``). Every control has a plain
method so tests drive it without a mouse.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

import snapfit
import snapmapfig

SHOW_LABELS = (("components", "Components"), ("fraction", "Fraction"),
               ("residual", "Residual"), ("spectrum", "Spectrum"))


class FitMixin:
    # -- state ----------------------------------------------------------------------
    def _fit_init(self):
        """Variables of the Fit tab (before any widget exists)."""
        self.fit_spec = snapfit.sanitise({})
        self.fit_numer = None              # None = the default (last of two or more)
        self.fit_on = None                 # None = whole map, else an area index
        self.fit_result = None             # the last snapfit.Result of a refinement
        self.v_fwin_lo = tk.StringVar(value="")
        self.v_fwin_hi = tk.StringVar(value="")
        self.v_fbg = tk.StringVar(value="shirley")
        self.v_fbin = tk.StringVar(value="1")
        self.v_fshared = tk.BooleanVar(value=False)
        self.v_fshow = {k: tk.BooleanVar(value=True) for k, _t in SHOW_LABELS}
        self.v_fon = tk.StringVar(value="Whole map")
        self._numer_vars = {}
        self._fit_cache = (None, None)

    def fit_view(self):
        """The fit part of the view, as ``snapmapfig`` stores it."""
        spec = dict(self.fit_spec, window=self._fit_window(),
                    background=self.v_fbg.get(), bin=self._fit_bin(),
                    shared_fwhm=bool(self.v_fshared.get()))
        return {"spec": spec,
                "show": [k for k, _t in SHOW_LABELS if self.v_fshow[k].get()],
                "numerator": self.fit_numer, "on": self.fit_on}

    def _apply_fit(self, fit):
        """Set the Fit tab from a saved view's fit part."""
        f = snapmapfig.sanitise_fit(fit)
        self.fit_spec = f["spec"]
        w = f["spec"]["window"]
        self.v_fwin_lo.set(f"{w[0]:.2f}" if w else "")
        self.v_fwin_hi.set(f"{w[1]:.2f}" if w else "")
        self.v_fbg.set(f["spec"]["background"])
        self.v_fbin.set(str(f["spec"]["bin"]))
        self.v_fshared.set(f["spec"]["shared_fwhm"])
        for k, _t in SHOW_LABELS:
            self.v_fshow[k].set(k in f["show"])
        self.fit_numer = f["numerator"]
        self.fit_on = f["on"]
        self.fit_result = None

    def _fit_window(self):
        try:
            a, b = float(self.v_fwin_lo.get()), float(self.v_fwin_hi.get())
        except ValueError:
            return None
        return [min(a, b), max(a, b)] if a != b else None

    def _fit_bin(self):
        try:
            return int(self.v_fbin.get())
        except ValueError:
            return 1

    # -- the tab ---------------------------------------------------------------------------
    def _build_fit_tab(self, tab):
        cols = ("name", "be", "fwhm", "rng", "free")
        self.fit_tree = ttk.Treeview(tab, columns=cols, show="headings",
                                     height=4, selectmode="browse")
        for c, text, w, a in (("name", "State", 70, "w"), ("be", "BE", 55, "e"),
                              ("fwhm", "FWHM", 50, "e"), ("rng", "±", 38, "e"),
                              ("free", "Free", 50, "w")):
            self.fit_tree.heading(c, text=text)
            self.fit_tree.column(c, width=w, anchor=a)
        self.fit_tree.pack(fill="x")
        self.fit_tree.bind("<Double-Button-1>", lambda e: self._fit_edit_selected())
        for k in range(snapfit.MAX_COMPONENTS):
            self.fit_tree.tag_configure(
                f"c{k}", foreground=snapmapfig.COMPONENT_COLOURS[
                    k % len(snapmapfig.COMPONENT_COLOURS)])
        row = ttk.Frame(tab)
        row.pack(fill="x", pady=(4, 0))
        self.fit_add_btn = ttk.Menubutton(row, text="Add ▾")
        self.fit_add_btn.pack(side="left")
        ttk.Button(row, text="Edit…", command=self._fit_edit_selected).pack(
            side="left", padx=3)
        ttk.Button(row, text="Remove", command=self._fit_remove_selected).pack(
            side="left")
        ttk.Button(row, text="Reseed", command=self.fit_seed).pack(
            side="left", padx=3)

        opts = ttk.Frame(tab)
        opts.pack(fill="x", pady=(6, 0))
        ttk.Label(opts, text="Window").grid(row=0, column=0, sticky="w")
        for c, v in ((1, self.v_fwin_lo), (2, self.v_fwin_hi)):
            e = ttk.Entry(opts, textvariable=v, width=7)
            e.grid(row=0, column=c, padx=2)
            e.bind("<Return>", lambda ev: self._fit_changed())
            e.bind("<FocusOut>", lambda ev: self._fit_changed())
        ttk.Label(opts, text="eV").grid(row=0, column=3, sticky="w")
        ttk.Label(opts, text="Background").grid(row=1, column=0, sticky="w")
        bg = ttk.Combobox(opts, state="readonly", width=9,
                          textvariable=self.v_fbg, values=("shirley", "linear"))
        bg.grid(row=1, column=1, columnspan=2, sticky="w", padx=2, pady=2)
        bg.bind("<<ComboboxSelected>>", lambda e: self._fit_changed())
        ttk.Label(opts, text="Bin").grid(row=2, column=0, sticky="w")
        b = ttk.Combobox(opts, state="readonly", width=4, textvariable=self.v_fbin,
                         values=("1", "2", "3"))
        b.grid(row=2, column=1, sticky="w", padx=2)
        b.bind("<<ComboboxSelected>>", lambda e: self._fit_changed())
        ttk.Checkbutton(opts, text="Same width", variable=self.v_fshared,
                        command=self._fit_changed).grid(
            row=2, column=2, columnspan=2, sticky="w")
        ttk.Label(opts, text="Fit on").grid(row=3, column=0, sticky="w")
        self.fon_cb = ttk.Combobox(opts, state="readonly", width=14,
                                   textvariable=self.v_fon,
                                   values=["Whole map"])
        self.fon_cb.grid(row=3, column=1, columnspan=3, sticky="w", padx=2,
                         pady=2)
        self.fon_cb.bind("<<ComboboxSelected>>", lambda e: self._fit_on_changed())

        show = ttk.Frame(tab)
        show.pack(fill="x", pady=(2, 0))
        for i, (k, text) in enumerate(SHOW_LABELS):
            ttk.Checkbutton(show, text=text, variable=self.v_fshow[k],
                            command=self._fit_changed).grid(
                row=i // 2, column=i % 2, sticky="w", padx=(0, 8))
        self.numer_btn = ttk.Menubutton(tab, text="Fraction of ▾")
        self.numer_btn.pack(anchor="w", pady=(4, 0))
        self.numer_menu = tk.Menu(self.numer_btn, tearoff=0)
        self.numer_btn["menu"] = self.numer_menu
        ttk.Button(tab, text="Fit spectrum", command=self.fit_spectrum).pack(
            anchor="w", pady=(8, 0))

        rcols = ("name", "be", "fwhm", "area", "pct")
        self.fit_res = ttk.Treeview(tab, columns=rcols, show="headings",
                                    height=4, selectmode="none")
        for c, text, w, a in (("name", "State", 62, "w"), ("be", "BE (eV)", 74, "e"),
                              ("fwhm", "FWHM", 62, "e"), ("area", "Area", 56, "e"),
                              ("pct", "%", 34, "e")):
            self.fit_res.heading(c, text=text)
            self.fit_res.column(c, width=w, anchor=a)
        self.fit_res.pack(fill="x", pady=(6, 0))
        self.fit_note = ttk.Label(tab, style="Muted.TLabel", wraplength=250,
                                  justify="left", text="")
        self.fit_note.pack(anchor="w", pady=(6, 0))
        self._refresh_fit_ui()

    # -- actions -----------------------------------------------------------------------------
    def _level_name(self):
        """The core level of the shown region ("Si 2p", "Cu 2p3/2" -> "Cu 2p")."""
        import re
        m = re.match(r"\s*([A-Z][a-z]?)\s*(\d[spdf])", self.region.name or "")
        return f"{m.group(1)} {m.group(2)}" if m else (self.region.name or "")

    def fit_seed(self):
        """Start the components from the tables for this core level (Si 2p:
        Si and SiO2; another level: its main line, on the measured peak)."""
        be, _hv, note = snapmapfig._fit_axis(self.region)
        if be is None:
            self.hint.configure(text=note + ".")
            return False
        data = self.cube.array3d()
        mean = data.reshape(-1, data.shape[2]).mean(axis=0)
        comps = snapfit.seed_components(self._level_name(), be=be, y=mean)
        if not comps:
            self.hint.configure(text=f"No starting components are known for "
                                     f"{self.region.name}: add some with Add.")
        self.fit_spec = snapfit.sanitise({"components": comps})
        if comps:
            w = snapfit.suggest_window(be, self.fit_spec["components"])
            self.v_fwin_lo.set(f"{w[0]:.2f}")
            self.v_fwin_hi.set(f"{w[1]:.2f}")
        self.fit_numer = None
        self.fit_result = None
        self._fit_changed()
        return bool(comps)

    def fit_add_library(self, name):
        """Add a state of this level from the chemical-state table."""
        for e in snapfit.library_states(self._level_name()):
            if e["name"] == name:
                return self._fit_add(snapfit.new_component(
                    e["name"], e["be"], e["fwhm"],
                    doublet=snapfit.doublet_of_level(self._level_name())))
        return None

    def fit_add_custom(self, name="New", be=None, fwhm=1.0):
        """Add a component of your own (at ``be``, default the middle of the
        fit window)."""
        if be is None:
            w = self._fit_window() or [self.cube.energy[-1], self.cube.energy[0]]
            be = (w[0] + w[1]) / 2.0
        return self._fit_add(snapfit.new_component(
            name, be, fwhm, doublet=snapfit.doublet_of_level(self._level_name())))

    def _fit_add(self, comp):
        comps = list(self.fit_spec["components"]) + [comp]
        if len(comps) > snapfit.MAX_COMPONENTS:
            self.hint.configure(
                text=f"At most {snapfit.MAX_COMPONENTS} components.")
            return None
        self.fit_spec = snapfit.sanitise(dict(self.fit_spec, components=comps))
        self.fit_result = None
        self._fit_changed()
        return len(self.fit_spec["components"]) - 1

    def fit_edit(self, i, **fields):
        """Change fields of component ``i`` (name, be, fwhm, be_range,
        free_be, free_fwhm, shape, doublet)."""
        comps = [dict(c) for c in self.fit_spec["components"]]
        if not 0 <= i < len(comps):
            return False
        comps[i].update(fields)
        self.fit_spec = snapfit.sanitise(dict(self.fit_spec, components=comps))
        self.fit_result = None
        self._fit_changed()
        return True

    def fit_remove(self, i):
        comps = [c for k, c in enumerate(self.fit_spec["components"]) if k != i]
        self.fit_spec = snapfit.sanitise(dict(self.fit_spec, components=comps))
        self.fit_numer = None
        self.fit_result = None
        self._fit_changed()

    def _fit_selected(self):
        sel = self.fit_tree.selection()
        return int(sel[0]) if sel else None

    def _fit_remove_selected(self):
        i = self._fit_selected()
        if i is not None:
            self.fit_remove(i)

    def _fit_edit_selected(self):
        i = self._fit_selected()
        if i is not None:
            self._fit_edit_dialog(i)

    def _fit_edit_dialog(self, i):
        from workbook_ui import _finish
        c = self.fit_spec["components"][i]
        win = tk.Toplevel(self)
        win.title(f"Component – {c['name']}")
        win.transient(self)
        vars_ = {
            "name": tk.StringVar(value=c["name"]),
            "be": tk.StringVar(value=f"{c['be']:.3f}"),
            "fwhm": tk.StringVar(value=f"{c['fwhm']:.3f}"),
            "be_range": tk.StringVar(value=f"{c['be_range']:.2f}"),
            "shape": tk.StringVar(value=c["shape"]),
            "split": tk.StringVar(value=f"{c['doublet']['split']:.3f}"
                                  if c["doublet"] else ""),
            "ratio": tk.StringVar(value=f"{c['doublet']['ratio']:.3f}"
                                  if c["doublet"] else ""),
            "free_be": tk.BooleanVar(value=c["free_be"]),
            "free_fwhm": tk.BooleanVar(value=c["free_fwhm"]),
        }
        rows = (("Name", "name"), ("Binding energy (eV)", "be"),
                ("FWHM (eV)", "fwhm"), ("May move ± (eV)", "be_range"),
                ("Shape", "shape"), ("Doublet split (eV)", "split"),
                ("Doublet area ratio", "ratio"))
        for r, (text, key) in enumerate(rows):
            ttk.Label(win, text=text).grid(row=r, column=0, sticky="w",
                                           padx=8, pady=2)
            ttk.Entry(win, textvariable=vars_[key], width=14).grid(
                row=r, column=1, padx=8, pady=2)
        ttk.Checkbutton(win, text="Refine the position",
                        variable=vars_["free_be"]).grid(
            row=len(rows), column=0, columnspan=2, sticky="w", padx=8)
        ttk.Checkbutton(win, text="Refine the width",
                        variable=vars_["free_fwhm"]).grid(
            row=len(rows) + 1, column=0, columnspan=2, sticky="w", padx=8)

        def num(key, default):
            try:
                return float(vars_[key].get())
            except ValueError:
                return default

        def ok():
            split, ratio = num("split", 0.0), num("ratio", 0.0)
            self.fit_edit(
                i, name=vars_["name"].get(), be=num("be", c["be"]),
                fwhm=num("fwhm", c["fwhm"]),
                be_range=num("be_range", c["be_range"]),
                shape=vars_["shape"].get().strip() or c["shape"],
                free_be=vars_["free_be"].get(),
                free_fwhm=vars_["free_fwhm"].get(),
                doublet=({"split": split, "ratio": ratio}
                         if split > 0 and ratio > 0 else None))
            win.destroy()
        bar = ttk.Frame(win, padding=8)
        bar.grid(row=len(rows) + 2, column=0, columnspan=2, sticky="e")
        ttk.Button(bar, text="OK", command=ok).pack(side="right")
        ttk.Button(bar, text="Cancel", command=win.destroy).pack(
            side="right", padx=6)
        _finish(win, self.app, 340, 330)

    # -- the spectrum to fit and the refinement ---------------------------------------------------
    def _fit_spectrum_data(self):
        """``(y, n_pixels, label)``: the mean spectrum of the whole map or of
        the area chosen in *Fit on*."""
        data = self.cube.array3d()
        masks = self.masks()
        on = self.fit_on
        if on is not None and on < len(masks) and masks[on] is not None:
            return (data[masks[on]].mean(axis=0), int(masks[on].sum()),
                    self.rois[on]["label"])
        return (data.reshape(-1, data.shape[2]).mean(axis=0),
                self.cube.nx * self.cube.ny, "whole map")

    def fit_spectrum(self):
        """Refine the free positions and widths on the mean spectrum of the
        whole map (or the area in *Fit on*); the maps are redrawn from the
        refined shapes. Returns the :class:`snapfit.Result` or None."""
        be, hv, note = snapmapfig._fit_axis(self.region)
        if be is None:
            self.hint.configure(text=note + ".")
            return None
        if not self.fit_spec["components"]:
            self.hint.configure(text="Add components first (Add, or Reseed).")
            return None
        y, npx, _label = self._fit_spectrum_data()
        try:
            res = snapfit.fit_mean(be, y, hv, self.fit_view()["spec"], n_pix=npx,
                                   gain=self._gain())
        except ValueError as e:
            self.hint.configure(text=str(e).capitalize() + ".")
            return None
        keep = self.fit_view()["spec"]
        self.fit_spec = snapfit.sanitise(dict(keep, components=res.spec["components"]))
        self.fit_result = res
        self._fit_changed()
        return res

    def _gain(self):
        """The map's noise gain (``snapfit.noise_gain``), found once."""
        key = id(self.cube)
        if getattr(self, "_gain_cache", (None, 1.0))[0] != key:
            self._gain_cache = (key, snapfit.noise_gain(self.cube))
        return self._gain_cache[1]

    def fit_maps(self):
        """The per-pixel partition for the current fit (``snapfit.FitMaps``),
        or None when there is no binding-energy axis or no components."""
        view = self.fit_view()
        key = repr((view["spec"], self.region.name, id(self.cube),
                    self.region.energy[:1]))
        if self._fit_cache[0] == key:
            return self._fit_cache[1]
        maps, _sp, _note = snapmapfig.fit_result(self.region, view)
        self._fit_cache = (key, maps)
        return maps

    def fit_fraction(self):
        """The fraction image of the chosen components, or None."""
        maps = self.fit_maps() if self.fit_view()["spec"]["components"] else None
        if maps is None:
            return None
        names = snapmapfig.sanitise_fit(self.fit_view())["numerator"]
        return snapfit.fraction(maps, names) if names else None

    # -- changes and what is shown -----------------------------------------------------------------
    def _fit_on_changed(self):
        label = self.v_fon.get()
        self.fit_on = None
        for k, r in enumerate(self.rois):
            if r["label"] == label:
                self.fit_on = k
        self._fit_changed()

    def _numer_changed(self):
        self.fit_numer = [n for n, v in self._numer_vars.items() if v.get()]
        self._fit_changed()

    def _fit_changed(self):
        self._fit_cache = (None, None)
        self._refresh_fit_ui()
        if self.mode.get() == "fit":
            self._rebuild(keep_view=True)

    def _refresh_fit_on(self):
        """The *Fit on* choices follow the areas."""
        if not hasattr(self, "fon_cb"):
            return
        names = ["Whole map"] + [r["label"] for r in self.rois]
        self.fon_cb.configure(values=names)
        if self.fit_on is not None and self.fit_on >= len(self.rois):
            self.fit_on = None
        self.v_fon.set(names[(self.fit_on or -1) + 1] if self.fit_on is not None
                       else "Whole map")

    def fit_results_rows(self):
        """``[(name, be, be_err, fwhm, fwhm_err, area, percent)]`` of the
        components on the spectrum chosen in *Fit on* (errors only after a
        refinement that left the components as they are)."""
        spec = self.fit_view()["spec"]
        if not spec["components"]:
            return []
        be, hv, _note = snapmapfig._fit_axis(self.region)
        if be is None:
            return []
        y, npx, _label = self._fit_spectrum_data()
        try:
            res = snapfit.evaluate(be, y, hv, spec, n_pix=npx,
                                   gain=self._gain())
        except ValueError:
            return []
        errs = {}
        if self.fit_result is not None and [
                (c["name"], round(c["be"], 6), round(c["fwhm"], 6))
                for c in self.fit_result.spec["components"]] == [
                (c["name"], round(c["be"], 6), round(c["fwhm"], 6))
                for c in spec["components"]]:
            errs = self.fit_result.errors
        total = float(res.amplitudes.sum()) or 1.0
        self._fit_chi2 = res.chi2red
        return [(c["name"], c["be"], errs.get((c["name"], "be")), c["fwhm"],
                 errs.get((c["name"], "fwhm")), float(a), 100.0 * float(a) / total)
                for c, a in zip(spec["components"], res.amplitudes)]

    def _refresh_fit_ui(self):
        if not hasattr(self, "fit_tree"):
            return
        self.fit_tree.delete(*self.fit_tree.get_children())
        for k, c in enumerate(self.fit_spec["components"]):
            free = "".join(x for x, f in (("BE ", c["free_be"]),
                                          ("W", c["free_fwhm"])) if f) or "–"
            self.fit_tree.insert("", "end", iid=str(k), tags=(f"c{k}",), values=(
                c["name"], f"{c['be']:.2f}", f"{c['fwhm']:.2f}",
                f"{c['be_range']:.1f}", free))
        self.fit_res.delete(*self.fit_res.get_children())
        self._fit_chi2 = None
        for name, be, be_e, fw, fw_e, area, pct in self.fit_results_rows():
            self.fit_res.insert("", "end", values=(
                name, f"{be:.2f}" + (f" ±{be_e:.2f}" if be_e is not None else ""),
                f"{fw:.2f}" + (f" ±{fw_e:.2f}" if fw_e is not None else ""),
                f"{area:.0f}", f"{pct:.0f}"))
        # the choices of the fraction
        self.numer_menu.delete(0, "end")
        self._numer_vars = {}
        chosen = snapmapfig.sanitise_fit(self.fit_view())["numerator"]
        for c in self.fit_spec["components"]:
            v = tk.BooleanVar(value=c["name"] in chosen)
            self._numer_vars[c["name"]] = v
            self.numer_menu.add_checkbutton(label=c["name"], variable=v,
                                            command=self._numer_changed)
        # the Add menu: this level's states from the table, or your own
        menu = tk.Menu(self.fit_add_btn, tearoff=0)
        for e in snapfit.library_states(self._level_name()):
            menu.add_command(label=f"{e['name']}   {e['be']:.2f} eV",
                             command=lambda n=e["name"]: self.fit_add_library(n))
        if menu.index("end") is not None:
            menu.add_separator()
        menu.add_command(label="Your own…",
                         command=lambda: self._fit_edit_new())
        self.fit_add_btn["menu"] = menu
        notes = []
        if self._fit_chi2 is not None:
            g = self._gain()
            notes.append(f"Reduced χ² {self._fit_chi2:.2f}"
                         + (f" (noise: one count is worth about {g:.0f} in "
                            f"the map's values)." if g > 3 or g < 0.3 else "."))
        if self.fit_result is not None:
            notes += [n.capitalize() + "." for n in self.fit_result.notes]
        self.fit_note.configure(text=" ".join(notes))

    def _fit_edit_new(self):
        i = self.fit_add_custom()
        if i is not None:
            self._fit_edit_dialog(i)
