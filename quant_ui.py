"""The desktop "Quantification" tab: live at%/RSF numbers from this
codebase's own recomputation of a CasaXPS fit (``quant.py`` /
``resultspages.py``), for whatever is currently ticked in the tree -- the
same computation and dedup rules (a region "counted once", a preferred-line
exclusion, the RSF fallback tiers) the PDF/deck/HTML browser already show,
just not routed through a report. Embedded as a tab beside Images / Stage
map / CasaXPS quant (``Workspace._build_side_tabs`` / ``_refresh_info``),
shown only when at least one ticked region has a fit.

With CasaXPS's own quantification files beside the data, a sample's
percentages are CasaXPS's (``casamatch``): the survey is its own entry in the
sample list ("<sample> (survey)") and the ticked high-resolution regions
another, never one total; "Use CasaXPS's own numbers" off recomputes
everything from the fits.

Unlike the CasaXPS-quant tab (``casaquant_ui.py``: CasaXPS's own *exported*
numbers, unconditionally shown when a folder had them), this tab has its own
RSF-library choice, independent of the Report generator's spec, so a
substitute can be explored live without touching report settings -- so it
calls ``resultspages.collect`` directly rather than through the memoized
``Workspace._results()``. A sample with no fit of its own (``sample.
casaxps`` only, i.e. ``levels`` empty) already has its home in the
CasaXPS-quant tab and is left out here rather than shown twice.

What the user can change here -- which regions count (the tick in the first
column), dividing the transmission function out, which depth level to show
(or the depth profile, in one of three modes) and writing the table as a CSV
-- lives in ``quantview.ViewState``; the numbers are recomputed from the
report's rows, never edited, and the report itself is not touched. The
ticks are keyed by what a region is (``quantview.entry_key``), so they
survive the tree being re-ticked and are saved, with the transmission choice
and the RSF library, in the workbook state.
"""

from __future__ import annotations

import csv
import tkinter as tk
from tkinter import filedialog, ttk

import quantview
import repeats
import reportspec
import resultspages
import rsf as rsf_lib

RSF_OFF = "off"
RSF_LABELS = {RSF_OFF: "Off (no substitute)", **rsf_lib.LIBRARIES}
RSF_CHOICES = tuple(RSF_LABELS)

PROFILE_LABEL = "Depth profile"
NOTES_MIN, NOTES_MAX, NOTES_WIDTH = 3, 10, 52   # lines; characters per line
CASA_COLUMNS = (0, 4, 5)        # region, CasaXPS %At, at % of composition_cells
TICK_ON, TICK_OFF = "☑", "☐"          # ballot box with / without check


class QuantPanel(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.samples = []
        self.sample = None
        self.view = quantview.ViewState()
        self._row_entry = {}            # tree row id -> entry index (ticks)
        self._level_of_row = None       # level index the rows belong to
        self._choices = []              # [(label, level index or None)]
        self.group = None               # name of the repeat group shown
        self._group_labels = {}         # sample-box text -> group name
        self._groups_dialog = None

        top = ttk.Frame(self)
        top.pack(side="top", fill="x", padx=6, pady=(6, 2))
        ttk.Label(top, text="Sample:").pack(side="left")
        self.sample_var = tk.StringVar()
        self.sample_box = ttk.Combobox(top, textvariable=self.sample_var,
                                       state="readonly", width=26)
        self.sample_box.pack(side="left", padx=(6, 12))
        self.sample_box.bind("<<ComboboxSelected>>",
                             lambda e: self._show_sample())

        ttk.Label(top, text="RSF:").pack(side="left")
        self.rsf_var = tk.StringVar(value=RSF_LABELS[RSF_OFF])
        self.rsf_box = ttk.Combobox(
            top, textvariable=self.rsf_var, state="readonly", width=32,
            values=[RSF_LABELS[k] for k in RSF_CHOICES])
        self.rsf_box.pack(side="left")
        self.rsf_box.bind("<<ComboboxSelected>>", lambda e: self.refresh())

        self.mid = mid = ttk.Frame(self)
        mid.pack(side="top", fill="x", padx=6, pady=(0, 2))
        ttk.Label(mid, text="Show:").pack(side="left")
        self.view_var = tk.StringVar()
        self.view_box = ttk.Combobox(mid, textvariable=self.view_var,
                                     state="readonly", width=26)
        self.view_box.pack(side="left", padx=(6, 12))
        self.view_box.bind("<<ComboboxSelected>>",
                           lambda e: self._view_chosen())
        self.mode_var = tk.StringVar(value=dict(
            quantview.PROFILE_MODES)[self.view.profile_mode])
        self.mode_box = ttk.Combobox(
            mid, textvariable=self.mode_var, state="disabled", width=28,
            values=[label for _k, label in quantview.PROFILE_MODES])
        self.mode_box.pack(side="left")
        self.mode_box.bind("<<ComboboxSelected>>",
                           lambda e: self._mode_chosen())

        self.bar = bar = ttk.Frame(self)
        bar.pack(side="top", fill="x", padx=6, pady=(0, 2))
        self.trans_var = tk.BooleanVar(value=False)
        self.trans_check = ttk.Checkbutton(
            bar, text="Divide out the transmission function",
            variable=self.trans_var, command=self._transmission_toggled,
            state="disabled")
        self.trans_check.pack(side="left")
        self.casa_var = tk.BooleanVar(value=True)
        self.casa_check = ttk.Checkbutton(
            bar, text="Use CasaXPS's own numbers", variable=self.casa_var,
            command=self._casa_toggled, state="disabled")
        self.casa_check.pack(side="left", padx=(12, 0))
        self.export_btn = ttk.Button(bar, text="Export CSV…",
                                     command=self.export_csv)
        self.export_btn.pack(side="right")
        self.groups_btn = ttk.Button(bar, text="Groups…",
                                     command=self.open_groups)
        self.groups_btn.pack(side="right", padx=(0, 6))
        self.reset_btn = ttk.Button(bar, text="Reset ticks",
                                    command=self.reset_ticks)
        self.reset_btn.pack(side="right", padx=(0, 6))

        tree_frame = ttk.Frame(self)
        tree_frame.pack(side="top", fill="both", expand=True, padx=6, pady=6)
        self.tree = ttk.Treeview(tree_frame, show="tree headings", height=10)
        vsb = ttk.Scrollbar(tree_frame, orient="vertical",
                            command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.tree.column("#0", width=34, minwidth=34, stretch=False,
                         anchor="center")
        self.tree.heading("#0", text="Use")
        self.tree.tag_configure("state", foreground="#5b6470")
        self.tree.tag_configure("off", foreground="#8a9099")
        self.tree.tag_configure("changed", font=("TkDefaultFont", 9, "bold"))
        self.tree.bind("<Button-1>", self._on_click)

        self.notes = tk.Text(self, height=NOTES_MIN, wrap="word", relief="flat",
                             state="disabled")
        self.notes.pack(side="top", fill="x", padx=6, pady=(0, 6))

    # -- choices -----------------------------------------------------------
    def rsf_choice(self):
        """The chosen RSF library key ("off" or one of ``rsf.LIBRARIES``)."""
        label = self.rsf_var.get()
        return next((k for k in RSF_CHOICES if RSF_LABELS[k] == label),
                    RSF_OFF)

    def set_transmission(self, on):
        """Set (and show) the transmission choice; a no-op for the numbers
        when no row has a transmission-corrected area."""
        self.view.transmission = bool(on)
        self.trans_var.set(bool(on))
        self._show_sample()

    def set_casa_numbers(self, on):
        """Take the percentages from CasaXPS's quantification files (the
        default, where a folder had them) or recompute them from the fits."""
        self.view.casa_numbers = bool(on)
        self.casa_var.set(bool(on))
        self.refresh()

    def _casa_toggled(self):
        self.set_casa_numbers(self.casa_var.get())

    def _transmission_toggled(self):
        self.view.transmission = bool(self.trans_var.get())
        self._show_sample()

    def _view_chosen(self):
        if self.sample is None:
            return
        label = self.view_var.get()
        li = next((i for t, i in self._choices if t == label), 0)
        self.view.level[self.sample.key] = li
        self._show_sample()

    def _mode_chosen(self):
        label = self.mode_var.get()
        self.view.profile_mode = next(
            (k for k, t in quantview.PROFILE_MODES if t == label), "element")
        self._show_sample()

    def reset_ticks(self):
        """Back to the report's own choice of what counts (in a group view,
        of every member)."""
        if self.group is not None:
            for k in self.app.ann.repeat_groups.get(self.group, {}).get(
                    "samples", []):
                self.view.reset(k)
        else:
            self.view.reset(self.sample.key if self.sample else None)
        self._show_sample()

    def reset_all_ticks(self):
        """Back to the automatic choice for every sample (the Report
        generator's reset)."""
        self.view.reset()
        self._show_sample()

    # -- data --------------------------------------------------------------
    def refresh(self):
        """Reload from the app's ticked regions and this panel's own RSF
        choice; call whenever the ticks, files or annotations may have
        changed (``Workspace._refresh_info``)."""
        app = self.app
        rsf_key = self.rsf_choice()
        rsf_table = app.rsf_entries() if rsf_key != RSF_OFF else None
        results = resultspages.collect(
            app.docs, app._display, lambda p: reportspec.doc_key(p),
            app.casa_quant, ticked=lambda r: id(r) in app.checked,
            rsf_table=rsf_table, rsf_library=rsf_key,
            prefer_csv=bool(app.csv_curves_var.get()),
            casa_numbers=self.view.casa_numbers)
        self.casa_check.configure(
            state="normal" if app.casa_quant else "disabled")
        # a sample with no fit of its own (CasaXPS's own export only) is
        # already shown in the "CasaXPS quant" tab -- not duplicated here
        self.samples = [s for s in results.samples if s.levels]
        self._set_box_values()
        has_t = quantview.transmission_available(self.samples)
        self.trans_check.configure(state="normal" if has_t else "disabled")
        if not has_t and self.view.transmission:
            self.view.transmission = False
            self.trans_var.set(False)
        if not self.samples:
            self.sample_var.set("")
            self.sample = None
            self._show_sample()
            return
        shown = next((lab for lab, n in self._group_labels.items()
                      if n == self.group), None)
        if shown is not None:               # keep showing the group
            self.sample_var.set(shown)
            self._show_sample()
            return
        self.group = None
        key = self.sample.key if self.sample else None
        match = next((s for s in self.samples if s.key == key), None) \
            or self.samples[0]
        self.sample = match
        self.sample_var.set(match.label)
        self._show_sample()

    def _set_box_values(self):
        """The samples, then one entry per saved repeat group that has at
        least one sample here."""
        keys = {s.key for s in self.samples}
        self._group_labels = {}
        for name, g in self.app.ann.repeat_groups.items():
            n = sum(1 for k in g["samples"] if k in keys)
            if n:
                plural = "s" if n != 1 else ""
                self._group_labels[f"Group: {name} ({n} sample{plural})"] = name
        self.sample_box["values"] = ([s.label for s in self.samples]
                                     + list(self._group_labels))

    def _show_group(self, name):
        """The mean, spread and count of a repeat group (``repeats.py``) with
        the user's ticks applied; no tick column, no depth-level choice."""
        self.group = name
        self.sample = None
        self._row_entry = {}
        self._level_of_row = None
        self._choices = []
        self._show_row(False)
        self.view_box["values"] = []
        self.view_var.set("")
        self.view_box.configure(state="disabled")
        self.mode_box.configure(state="disabled")
        g = quantview.repeat_groups(
            self.samples, {name: self.app.ann.repeat_groups[name]},
            self.view)[0]
        header, rows = quantview.repeat_table(g)
        self.tree.configure(show="headings")
        self._set_columns(header)
        self._fill([(k, c) for k, c, _p in rows])
        for iid, (_k, _c, partial) in zip(self.tree.get_children(), rows):
            if partial:
                self.tree.item(iid, tags=tuple(self.tree.item(iid, "tags"))
                               + ("changed",))
        self._set_notes(list(g.notes) + [
            "Each value is that sample's own composition with your ticks "
            "applied. A bold row was found in fewer samples than the group "
            "has; a mean over one sample has no spread."])

    def open_groups(self):
        """The repeat-group editor (one window at a time)."""
        d = self._groups_dialog
        try:
            if d is not None and d.winfo_exists():
                d.lift()
                return d
        except tk.TclError:
            pass
        self._groups_dialog = RepeatGroupsDialog(self, self.app)
        return self._groups_dialog

    def _show_sample(self):
        label = self.sample_var.get()
        gname = self._group_labels.get(label)
        if gname is not None:
            self._show_group(gname)
            return
        self.group = None
        self.sample = next((s for s in self.samples if s.label == label),
                           None)
        self._row_entry = {}
        self._level_of_row = None
        if self.sample is None:
            self._choices = []
            self._show_row(False)
            self.view_box["values"] = []
            self.view_var.set("")
            self.view_box.configure(state="disabled")
            self.mode_box.configure(state="disabled")
            self.tree.configure(show="headings")
            self._set_columns(())
            self._fill([])
            self._set_notes([])
            return
        s = self.sample
        shown = self.view.shown(s)
        self._choices = ([(PROFILE_LABEL, None)] if s.is_profile else []) + [
            (quantview.level_label(lv, i), i) for i, lv in enumerate(s.levels)]
        self.view_box["values"] = [t for t, _i in self._choices]
        self.view_box.configure(
            state="readonly" if len(self._choices) > 1 else "disabled")
        self._show_row(len(self._choices) > 1)
        self.view_var.set(next(t for t, i in self._choices if i == shown))
        self.mode_box.configure(state="readonly" if shown is None
                                else "disabled")
        self.mode_var.set(dict(quantview.PROFILE_MODES)[
            self.view.profile_mode])
        if shown is None:
            header, rows = quantview.profile_cells(s, self.view)
            self.tree.configure(show="headings")        # no tick column
            self._set_columns(header)
            self._fill([(None, r) for r in rows])
        else:
            self._level_of_row = shown
            eff = quantview.effective(s, shown, self.view)
            changed = quantview.changed(s, shown, self.view)
            self.tree.configure(show="tree headings")
            stats = self.app.fit_stats()
            header = resultspages.composition_header(eff, stats)
            cells = resultspages.composition_cells(eff, stats)
            if header is resultspages.CASA_HEADER:
                # CasaXPS's own numbers: the region and the two percentages;
                # background, RSF, area and fit quality are the recomputed
                # view's business and would push these off the column
                keep = CASA_COLUMNS
                header = tuple(header[i] for i in keep)
                cells = [(k, [c[i] for i in keep]) for k, c in cells]
            self._set_columns(header)
            self._fill(cells, ticks=(eff.include, changed))
        self._set_notes(list(s.notes) + quantview.transmission_notes(
            s, self.view))

    def _show_row(self, on):
        """The "Show:" row (depth level / profile, profile mode) is only there
        when the sample has more than one level to choose between."""
        shown = bool(self.mid.winfo_manager())
        if on and not shown:
            self.mid.pack(side="top", fill="x", padx=6, pady=(0, 2),
                          before=self.bar)
        elif not on and shown:
            self.mid.pack_forget()

    # -- ticks ---------------------------------------------------------------
    def _on_click(self, event):
        if self.tree.identify_region(event.x, event.y) != "tree":
            return
        if self.tree.identify_column(event.x) != "#0":
            return
        iid = self.tree.identify_row(event.y)
        if iid in self._row_entry:
            self.toggle_row(iid)
            return "break"

    def toggle_row(self, iid):
        """Flip the tick of a region row (a state row has none)."""
        ei = self._row_entry.get(iid)
        if ei is None or self.sample is None or self._level_of_row is None:
            return
        quantview.toggle_entry(self.sample, self._level_of_row, ei, self.view)
        self._show_sample()

    # -- output --------------------------------------------------------------
    def export_csv(self, path=None):
        """Write the table of every sample here (the ticks and the
        transmission choice applied) as a CSV; asks for the file unless
        ``path`` is given. Returns the path written, or None."""
        if not self.samples:
            return None
        in_group = self.group is not None
        if in_group:
            rows = repeats.csv_rows(quantview.repeat_groups(
                self.samples, self.app.ann.repeat_groups, self.view))
        if path is None:
            path = filedialog.asksaveasfilename(
                parent=self, defaultextension=".csv",
                initialfile=("quantification_repeat_groups.csv" if in_group
                             else "quantification.csv"),
                filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return None
        if not in_group:
            rows = quantview.csv_table(self.samples, self.view)
        # UTF-8 with a BOM so Excel reads the Greek/degree signs
        with open(path, "w", encoding="utf-8-sig", newline="") as fh:
            csv.writer(fh).writerows(rows)
        status = getattr(self.app, "status", None)
        if status is not None:
            status.config(text=f"Saved {path}")
        return path

    # -- table ---------------------------------------------------------------
    def _set_columns(self, header):
        cols = list(range(len(header)))
        self.tree["columns"] = cols
        for i, text in enumerate(header):
            self.tree.heading(i, text=text)
            self.tree.column(i, anchor="w" if i == 0 else "center",
                             width=140 if i == 0 else 100, stretch=True)

    def _fill(self, rows, ticks=None):
        self.tree.delete(*self.tree.get_children())
        self._row_entry = {}
        ei = -1
        for kind, cells in rows:
            tags, text = (("state",) if kind == "state" else ()), ""
            if ticks is not None and kind == "region":
                ei += 1
                inc, changed = ticks
                text = TICK_ON if inc[ei] else TICK_OFF
                if not inc[ei]:
                    tags += ("off",)
                if changed[ei]:
                    tags += ("changed",)
            iid = self.tree.insert("", "end", text=text, values=cells,
                                   tags=tags)
            if ticks is not None and kind == "region":
                self._row_entry[iid] = ei

    def _set_notes(self, notes):
        lines = sum(max(1, -(-len(n) // NOTES_WIDTH)) for n in notes)
        self.notes.configure(height=max(NOTES_MIN, min(NOTES_MAX, lines)))
        self.notes.configure(state="normal")
        self.notes.delete("1.0", "end")
        self.notes.insert("end", "\n".join(notes))
        self.notes.configure(state="disabled")


class RepeatGroupsDialog(tk.Toplevel):
    """Create and edit repeat groups (``repeats.py``): which samples are
    repeats of one another. Every change is saved to the workbook's
    annotations at once and the Quantification tab redrawn. The behaviour is
    in plain methods (``new_group``, ``select``, ``set_name``, ``set_kind``,
    ``toggle_sample``, ``suggest``, ``delete_group``) that the widgets call,
    so tests need no mouse.

    A group with no sample is kept only in this window (the annotations store
    none); a sample key that is no longer among the quantified samples stays
    in its group, listed as missing, until it is unticked."""

    def __init__(self, panel, app):
        super().__init__(panel)
        self.panel, self.app = panel, app
        self.title("Repeat groups")
        self.drafts = {name: {"kind": g["kind"],
                              "samples": list(g["samples"])}
                       for name, g in app.ann.repeat_groups.items()}
        self.current = next(iter(self.drafts), None)
        self._build()
        self._fill_groups()
        import workbook_ui
        workbook_ui._finish(self, app, 760, 460)
        workbook_ui.keep_in_front(self, panel.winfo_toplevel())

    # -- what can be a member -------------------------------------------------
    def _reels_samples(self):
        """The samples with a stored REELS π / (σ+π) ratio (``reelstable``): a
        REELS-only sample has no quantification, so the Quantification tab's
        own list never holds it."""
        import reelstable
        try:
            return reelstable.suggest_samples(self.app._reels_table())
        except AttributeError:                      # an app with no REELS
            return []

    def candidates(self):
        """``[(key, label, usable)]``: the quantified samples, then the keys
        a group names that are no longer there. A depth profile is not a
        repeat, so it is listed but cannot be ticked."""
        out = [(s.key, s.label, bool(s.levels) and not s.is_profile)
               for s in self.panel.samples]
        have = {k for k, _l, _u in out}
        for s in self._reels_samples():             # REELS spectra, which
            if s.key not in have:                   # may have no fit at all
                have.add(s.key)
                out.append((s.key, f"{s.label} (REELS)", True))
        for g in self.drafts.values():
            for k in g["samples"]:
                if k not in have:
                    have.add(k)
                    out.append((k, f"{k} (missing)", True))
        return out

    # -- behaviour --------------------------------------------------------------
    def _commit(self, name):
        g = self.drafts.get(name)
        if g and g["samples"]:
            self.app.ann.set_repeat_group(name, g["kind"], g["samples"])
        else:
            self.app.ann.remove_repeat_group(name)
        self.app._ann_changed(relabel=False)
        self.panel.refresh()

    def new_group(self, name=None):
        base = (name or "Group").strip() or "Group"
        name, n = base, 1
        while name in self.drafts:
            n += 1
            name = f"{base} {n}"
        self.drafts[name] = {"kind": repeats.DEFAULT_KIND, "samples": []}
        self.current = name
        self._fill_groups()
        return name

    def select(self, name):
        if name in self.drafts:
            self.current = name
            self._fill_groups()

    def set_name(self, new):
        """Rename the current group; False when the name is empty or taken."""
        new = (new or "").strip()
        old = self.current
        if old is None or not new or (new != old and new in self.drafts):
            return False
        if new != old:
            self.drafts = {(new if k == old else k): v
                           for k, v in self.drafts.items()}
            self.app.ann.remove_repeat_group(old)
            self.current = new
            self._commit(new)
            self._fill_groups()
        return True

    def set_kind(self, kind):
        g = self.drafts.get(self.current)
        if g is None or kind not in repeats.KINDS:
            return False
        g["kind"] = kind
        self._commit(self.current)
        return True

    def toggle_sample(self, key):
        """Add or remove a sample; False for one that cannot be a member."""
        g = self.drafts.get(self.current)
        usable = {k: u for k, _l, u in self.candidates()}
        if g is None or key not in usable:
            return False
        if key in g["samples"]:
            g["samples"].remove(key)
        elif usable[key]:
            g["samples"].append(key)
        else:
            return False
        self._commit(self.current)
        self._fill_members()
        return True

    def suggest(self):
        """Add the groups the sample names suggest; none that exist already
        and no sample list is changed. Returns the names added."""
        added = []
        have = {s.key for s in self.panel.samples}
        samples = list(self.panel.samples) + [
            s for s in self._reels_samples() if s.key not in have]
        for name, g in repeats.suggest(samples).items():
            if name not in self.drafts:
                self.drafts[name] = g
                self._commit(name)
                added.append(name)
        if added:
            self.current = self.current or added[0]
            self._fill_groups()
        return added

    def delete_group(self, name=None):
        name = name or self.current
        if name not in self.drafts:
            return
        del self.drafts[name]
        self.app.ann.remove_repeat_group(name)
        self.current = next(iter(self.drafts), None)
        self.app._ann_changed(relabel=False)
        self.panel.refresh()
        self._fill_groups()

    # -- widgets ------------------------------------------------------------------
    def _build(self):
        left = ttk.Frame(self)
        left.pack(side="left", fill="y", padx=(8, 4), pady=8)
        ttk.Label(left, text="Groups").pack(anchor="w")
        self.listbox = tk.Listbox(left, width=24, exportselection=False,
                                  height=10)
        self.listbox.pack(fill="y", expand=True, pady=(2, 6))
        self.listbox.bind("<<ListboxSelect>>", self._picked)
        row = ttk.Frame(left)
        row.pack(fill="x")
        ttk.Button(row, text="New", command=self.new_group).pack(side="left")
        ttk.Button(row, text="Delete", command=self.delete_group).pack(
            side="left", padx=(6, 0))
        ttk.Button(left, text="Suggest from names",
                   command=self.suggest).pack(fill="x", pady=(8, 0))

        right = ttk.Frame(self)
        right.pack(side="left", fill="both", expand=True, padx=(4, 8), pady=8)
        form = ttk.Frame(right)
        form.pack(fill="x")
        ttk.Label(form, text="Name").grid(row=0, column=0, sticky="w")
        self.name_var = tk.StringVar()
        name = ttk.Entry(form, textvariable=self.name_var)
        name.grid(row=0, column=1, sticky="we", padx=6, pady=2)
        name.bind("<Return>", lambda e: self.set_name(self.name_var.get()))
        name.bind("<FocusOut>", lambda e: self.set_name(self.name_var.get()))
        ttk.Label(form, text="Kind").grid(row=1, column=0, sticky="w")
        self.kind_var = tk.StringVar()
        self.kind_box = ttk.Combobox(
            form, textvariable=self.kind_var, state="readonly",
            values=list(repeats.KINDS.values()))
        self.kind_box.grid(row=1, column=1, sticky="we", padx=6, pady=2)
        self.kind_box.bind("<<ComboboxSelected>>", lambda e: self.set_kind(
            next((k for k, t in repeats.KINDS.items()
                  if t == self.kind_var.get()), None)))
        form.columnconfigure(1, weight=1)
        ttk.Label(right, text="Samples in this group").pack(
            anchor="w", pady=(8, 2))
        frame = ttk.Frame(right)
        frame.pack(fill="both", expand=True)
        self.members = ttk.Treeview(frame, show="tree", selectmode="none")
        vsb = ttk.Scrollbar(frame, orient="vertical",
                            command=self.members.yview)
        self.members.configure(yscrollcommand=vsb.set)
        self.members.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.members.tag_configure("off", foreground="#8a9099")
        self.members.bind("<Button-1>", self._member_click)
        self._member_key = {}
        ttk.Label(right, foreground="#5b6470", wraplength=480,
                  text="Renaming a sample takes it out of its group. "
                       "Depth profiles are not repeats.").pack(
            anchor="w", pady=(6, 0))
        ttk.Button(right, text="Close", command=self.destroy).pack(
            anchor="e", pady=(6, 0))

    def _picked(self, _event=None):
        sel = self.listbox.curselection()
        if sel:
            self.select(self.listbox.get(sel[0]))

    def _fill_groups(self):
        self.listbox.delete(0, "end")
        for name in self.drafts:
            self.listbox.insert("end", name)
        if self.current in self.drafts:
            i = list(self.drafts).index(self.current)
            self.listbox.selection_set(i)
        g = self.drafts.get(self.current)
        self.name_var.set(self.current or "")
        self.kind_var.set(repeats.KINDS[g["kind"]] if g else "")
        self._fill_members()

    def _fill_members(self):
        self.members.delete(*self.members.get_children())
        self._member_key = {}
        g = self.drafts.get(self.current)
        if g is None:
            return
        for key, label, usable in self.candidates():
            on = key in g["samples"]
            text = (TICK_ON if on else TICK_OFF) + " " + label
            if not usable:
                text += "  (depth profile: not a repeat)"
            iid = self.members.insert("", "end", text=text,
                                      tags=() if usable else ("off",))
            self._member_key[iid] = key

    def _member_click(self, event):
        iid = self.members.identify_row(event.y)
        if iid in self._member_key:
            self.toggle_sample(self._member_key[iid])
        return "break"
