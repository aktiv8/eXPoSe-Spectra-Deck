"""Plain numbers in columns (``.csv`` / ``.asc`` / ``.txt`` / ``.dat``): the
files journals such as *Surface Science Spectra* ask for, and what most
programs export.

The parsing, the column guesses and the options are in ``columntext`` (Tk-free,
shared with the Annotations that store how a file was read); this reader turns
them into regions. ``load(path, options)`` takes the options the import dialog
produced, else reads the file the way ``columntext.guess_options`` says. One
region is made per intensity column. A file carries no photon energy, pass
energy or sample, so those are what the user enters (or they stay unset).
"""

from __future__ import annotations

import os

import columntext

from .base import Region, SpectrumFile, read_bytes

PATTERNS = tuple("*" + e for e in columntext.EXTS)
HEAD = 4096


def sniff(head: bytes, ext: str) -> bool:
    """A file of numbers in columns, one of which runs steadily one way (the
    energy axis), with at least ``columntext.MIN_ROWS`` lines. Strict, so a
    results table or a settings file in a folder is not taken for a spectrum."""
    if ext not in columntext.EXTS or b"\x00" in head[:512]:
        return False
    lines = columntext.decode(head).splitlines()
    if len(head) >= HEAD and lines:
        lines = lines[:-1]                  # the line the read cut in two
    table = columntext.read_table("\n".join(lines))
    return table is not None and columntext.find_energy_col(table)[1]


def _median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else None


class ColumnTextFile(SpectrumFile):
    format_name = "Column text (CSV / ASC / TXT)"

    def __init__(self):
        super().__init__()
        self.import_options = {}        # what the file was read with (stored)
        self.table = None               # the columntext.Table, for the dialog

    def load(self, path: str, options=None):
        self.path = path
        table = columntext.read_table(columntext.decode(read_bytes(path)))
        if table is None:
            raise ValueError(
                f"no table of numbers found (it needs at least "
                f"{columntext.MIN_ROWS} lines of numbers in columns)")
        self.table = table
        merged = columntext.guess_options(table, path)
        if isinstance(options, dict):
            merged.update({k: v for k, v in options.items() if k in merged})
        opts = columntext.sanitise_options(merged, table.n_cols)
        self.import_options = opts
        if table.trailing:
            self.warnings.append(
                f"{table.trailing} line(s) after the table, from line "
                f"{table.trailing_at}, were ignored (not numbers in "
                f"{table.n_cols} columns).")
        if not columntext.find_energy_col(table)[1]:
            self.warnings.append("The energy column is not evenly spaced or "
                                 "does not run one way; check the columns.")
        if table.preamble:
            self.instrument["Comments"] = " | ".join(table.preamble[:5])
        stem = os.path.splitext(os.path.basename(path))[0]
        for k, col in enumerate(opts["intensity_cols"]):
            self._add_region(k, col, opts, table, stem)
        return self._finish()

    def _add_region(self, k, col, opts, table, stem):
        energy = list(table.columns[opts["energy_col"]])
        counts = list(table.columns[col])
        hv, axis = opts["photon_energy"], opts["axis"]
        notes = []
        if axis == "KE" and hv:
            energy = [hv - e for e in energy]
            label = "Binding Energy"
            notes.append(f"Kinetic energies converted to binding energy with "
                         f"hν = {hv:g} eV.")
        elif axis == "KE":
            label = "Kinetic Energy"
            self.warnings.append(
                f"{os.path.basename(self.path)}: kinetic energies, and no "
                f"photon energy to convert them: the kinetic axis is kept.")
        else:
            label = "Binding Energy"
        descending = label == "Binding Energy"          # BE high to low
        if energy[0] != energy[-1] and (energy[0] < energy[-1]) == descending:
            energy, counts = energy[::-1], counts[::-1]
            notes.append("Stored in the opposite order to the file "
                         "(binding energy runs high to low).")
        step = _median([abs(b - a) for a, b in zip(energy, energy[1:])
                        if b != a])
        names = opts["names"]
        raw = names[k] if k < len(names) else ""
        name = columntext.region_name(raw, energy)
        if not name:
            # never guessed from the energies here: the import dialog may
            # suggest a core level, the user decides
            name = stem if len(opts["intensity_cols"]) == 1 \
                else f"{stem} {k + 1}"
        reg = Region(
            name=name, index=k, offset=k + 1, technique="XPS",
            energy=energy, counts=counts, energy_label=label,
            energy_units="eV", count_label="Intensity",
            count_units=opts["units"], decodable=True,
            note=" ".join(notes), sample=opts["sample"], photon_energy=hv,
            pass_energy=opts["pass_energy"],
            step=round(step, 6) if step else None, anode=opts["anode"])
        reg.extra["import_axis"] = axis
        self.regions.append(reg)
