"""Plain numbers in columns (``.csv`` / ``.asc`` / ``.txt`` / ``.dat``, the
SPECS ``.xy`` export and Excel ``.xlsx`` workbooks): the files journals such as
*Surface Science Spectra* ask for, and what most programs export.

The parsing, the column guesses and the options are in ``columntext`` (Tk-free,
shared with the Annotations that store how a file was read); this reader turns
them into regions. ``load(path, options)`` takes the options the import dialog
produced, else reads the file the way ``columntext.guess_options`` says. One
region is made per intensity column. A file carries no photon energy, pass
energy or sample, so those are what the user enters (or they stay unset),
except where a header states them (the SPECS ``.xy`` export gives the
excitation and the pass energy, the dwell time, the number of scans, the
analyser, its lens and scan mode and the work function; those are used as
stated). The wide CSV this program writes comes back one spectrum per
energy / intensity pair, with its names, axes and units.
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
    if ext in columntext.XLSX_EXTS:
        # a workbook; whether a sheet holds a spectrum is found on loading
        return head[:2] == b"PK" and b"xl/" in head
    if ext not in columntext.EXTS or b"\x00" in head[:512]:
        return False
    lines = columntext.decode(head).splitlines()
    if len(head) >= HEAD and lines:
        lines = lines[:-1]                  # the line the read cut in two
    table = columntext.read_any("\n".join(lines))
    if table is None or not columntext.find_energy_col(table)[1]:
        return False
    # a bare ".xy" is also an X-ray diffraction pattern: only the SPECS header
    # (excitation and pass energy) says it is a photoelectron spectrum
    return ext != ".xy" or columntext.is_specs(table.facts)


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
        want = options.get("sheet") if isinstance(options, dict) else None
        table = columntext.load_table(path, want or None)
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
        specs = columntext.is_specs(table.facts)
        if table.trailing:
            more = (" (a SPECS file can hold several blocks of data; only "
                    "the longest was read)" if specs else "")
            self.warnings.append(
                f"{table.trailing} line(s) after the table, from line "
                f"{table.trailing_at}, were ignored (not numbers in "
                f"{table.n_cols} columns){more}.")
        if table.sheet and len(table.sheets) > 1:
            others = [f"{n} ({r} rows)" for n, r in table.sheets
                      if n != table.sheet and r]
            if others:
                self.warnings.append(
                    f"Sheet '{table.sheet}' was read; the other sheet(s) "
                    f"with a table were not: {', '.join(others)}.")
        if not columntext.find_energy_col(table)[1]:
            self.warnings.append("The energy column is not evenly spaced or "
                                 "does not run one way; check the columns.")
        if specs:
            self._specs_facts(table)
        elif table.preamble:
            self.instrument["Comments"] = " | ".join(table.preamble[:5])
        stem = os.path.splitext(os.path.basename(path))[0]
        if table.pairs:
            for k, pr in enumerate(table.pairs):
                self._add_region(k, pr["energy"], pr["intensity"], opts,
                                 table, stem, pair=pr)
        else:
            for k, col in enumerate(opts["intensity_cols"]):
                self._add_region(k, opts["energy_col"], col, opts, table, stem)
        return self._finish()

    def _specs_facts(self, table):
        """What a SPECS ``.xy`` header states about the instrument."""
        f = table.facts
        got = {"Instrument": f.get("analyzer", ""),
               "Work function (eV)": f.get("eff. workfunction", ""),
               "Comments": f.get("comment", "")}
        self.instrument.update({k: v for k, v in got.items() if v})
        n = columntext.fact_float(f, "values/curve")
        if n is not None and int(n) != table.n_rows:
            self.warnings.append(
                f"The header says {int(n)} values per curve but {table.n_rows}"
                " were read.")

    _ANALYSER_MODES = {"fixedanalyzertransmission": "Constant analyser energy "
                       "(CAE)", "fixedretardingratio": "Constant retard ratio "
                       "(CRR)"}

    def _add_region(self, k, ecol, icol, opts, table, stem, pair=None):
        energy = list(table.columns[ecol])
        counts = list(table.columns[icol])
        hv = opts["photon_energy"]
        axis = pair["axis"] if pair else opts["axis"]
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
        raw = pair["name"] if pair else (names[k] if k < len(names) else "")
        name = columntext.region_name(raw, energy)
        if not name:
            # never guessed from the energies here: the import dialog may
            # suggest a core level, the user decides
            many = len(table.pairs) if pair else len(opts["intensity_cols"])
            name = stem if many == 1 else f"{stem} {k + 1}"
        f = table.facts if columntext.is_specs(table.facts) else {}
        scans = columntext.fact_float(f, "number of scans")
        dwell = columntext.fact_float(f, "dwell time")
        reg = Region(
            name=name, index=k, offset=k + 1, technique="XPS",
            energy=energy, counts=counts, energy_label=label,
            energy_units="eV", count_label="Intensity",
            count_units=pair["units"] if pair else opts["units"],
            decodable=True,
            note=" ".join(notes),
            sample=opts["sample"] or (pair["sample"] if pair else ""),
            photon_energy=hv, pass_energy=opts["pass_energy"],
            step=round(step, 6) if step else None, anode=opts["anode"],
            dwell=dwell, lens_mode=f.get("analyzer lens", ""))
        if scans and scans >= 1 and scans == int(scans):
            reg.extra["n_scans"] = int(scans)
        mode = self._ANALYSER_MODES.get(f.get("scan mode", "").lower())
        if mode:
            reg.extra["analyser_mode"] = mode
        if f:
            reg.extra["info"] = dict(f)
        reg.extra["import_axis"] = axis
        self.regions.append(reg)
