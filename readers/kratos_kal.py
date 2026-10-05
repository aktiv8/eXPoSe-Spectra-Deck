"""Kratos Vision ``.kal`` reader (ASCII "dump_dataset" text).

A ``.kal`` file holds one or more objects, each introduced by
``Object name = ...`` and made of numbered ``id Name = value`` fields
(``3 Spectrum scan start``, ``4 step size``, ``12 Ordinate values = {...}``,
``42 Pass energy``, ``3113/3114`` chemical symbol / transition, the
``Transmission Function (ke,t)`` lists, ...). The abscissa is kinetic energy;
binding energy = anode energy - KE.
"""

from __future__ import annotations

import re

import kratosmap

from .base import (Region, SpectrumFile, analyser_mode_name, canon_region_name,
                   guess_region_name, read_bytes)

# X-ray line energies (eV) by the anode named in Kratos flags
_ANODES = {"AL": 1486.6, "MG": 1253.6, "AG": 2984.2, "ZR": 2042.4,
           "TI": 4510.9, "CR": 5414.7}
_FIELD_RE = re.compile(r"^\s*(\d+)\s+(.+?)\s*=\s*(.*)$")


def sniff(head: bytes, ext: str) -> bool:
    return head.lstrip().startswith(b"Dataset filename") and b"Object name" in head


def _num(s, default=None):
    m = re.search(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s or "")
    return float(m.group(0)) if m else default


def _list(value):
    body = value.strip().lstrip("{").rstrip("}")
    return [float(t) for t in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", body)]


def _date(s):
    m = re.match(r"^(\d{2})/(\d{2})/(\d{2})\s+(\d{2}:\d{2}:\d{2})", s or "")
    if not m:
        return (s or "").strip()
    yy = int(m.group(1))
    return f"{1900 + yy if yy >= 70 else 2000 + yy}-{m.group(2)}-{m.group(3)} {m.group(4)}"


class KratosKalFile(SpectrumFile):
    format_name = "Kratos Vision (.kal)"

    def load(self, path: str):
        self.path = path
        text = read_bytes(path).decode("latin-1").replace("\r", "")
        objects, cur, pending = [], None, None
        for line in text.split("\n"):
            if pending is not None:                    # multi-line { ... } list
                pending[1].append(line)
                if "}" in line:
                    cur[pending[0]] = " ".join(pending[1])
                    pending = None
                continue
            if line.startswith("Dataset filename"):
                self._dataset = line.split("=", 1)[1].strip()
                continue
            if line.startswith("Object name"):
                cur = {"_object": line.split("=", 1)[1].strip()}
                objects.append(cur)
                continue
            m = _FIELD_RE.match(line)
            if not m or cur is None:
                continue
            key, val = m.group(2).strip(), m.group(3).strip()
            if val == "{":                             # container: ignore
                continue
            if val.startswith("{") and "}" not in val:
                pending = (key, [val])
                continue
            cur[key] = val
        return self._from_objects(objects)

    def _from_objects(self, objects):
        """Turn the dumped objects (``{field name: value text}``, each with
        ``_object``) into regions. Shared with the binary ``.dset`` reader,
        which rebuilds the same text from the Vision2 file itself."""
        if not objects:
            raise ValueError("no 'Object name' entries found")
        first = objects[0]
        self._no_hv, self._no_hv_maps = [], []
        self._position = ""              # the last "Sample Position" object
        for i, o in enumerate(objects):
            if o.get("Stage Position Name"):
                self._position = o["Stage Position Name"].strip()
            self._add_object(i, o)
        for names, one, many in ((self._no_hv, "spectrum", "spectra"),
                                 (self._no_hv_maps, "map", "maps")):
            if names:                    # one line for the file, not one a region
                shown = ", ".join(names[:4]) + (f" … ({len(names)} in all)"
                                                if len(names) > 4 else "")
                self.warnings.append(
                    f"X-ray energy unknown for {len(names)} "
                    f"{one if len(names) == 1 else many} ({shown}): "
                    "kinetic-energy axis.")
        anode = self._anode(first)
        self.instrument = {k: v for k, v in {
            "Instrument": "Kratos (Vision)",
            "Acquisition software": "Kratos Vision",
            "X-ray source": (f"{anode[0]} ({anode[1]:g} eV)" if anode else ""),
            "Lens mode": first.get("HSA Lens Mode", "").replace(
                "F_HSA_LENS_", "").title(),
        }.items() if v}
        return self._finish()

    @staticmethod
    def _anode(o):
        for key in ("Xray Reference Energy", "Xray Gun Anode"):
            v = o.get(key, "")
            for sym, e in _ANODES.items():
                if re.search(rf"(_|^){sym}(_|$)", v.upper()):
                    return f"{sym.title()} anode", e
        return None

    @staticmethod
    def _aperture(o):
        """'Slot' from the aperture-size descriptor, plus the iris position
        when it says something else."""
        size = o.get("Descriptor for aperture size used in acquisition",
                     "").strip()
        iris = o.get("Descriptor for iris position used in acquisition",
                     "").strip()
        if iris and iris.lower() != size.lower():
            return f"{size}, iris {iris}" if size else f"iris {iris}"
        return size

    @staticmethod
    def _neutraliser(o):
        """The neutraliser as the file states it. Units are not given for the
        filament current, bias and balance, so none are written."""
        state = o.get("Neutraliser Switch State", "")
        if not state or "OFF" in state.upper():
            return ""
        words = state.replace("F_NEUTRALISER_", "").replace("_", " ").lower()
        parts = [f"{label} {o[key].strip()}" for label, key in (
            ("filament current", "Charge Neutraliser Filament Current"),
            ("bias", "Charge Neutraliser Filament Bias"),
            ("balance", "Charge Neutraliser Charge Balance")) if o.get(key)]
        return words + (": " + ", ".join(parts) if parts else "")

    def _add_object(self, idx, o):
        if kratosmap.is_map(o):
            self._add_map(idx, o)
            return
        vals = _list(o.get("Ordinate values", ""))
        start = _num(o.get("Spectrum scan start"))
        step = _num(o.get("Spectrum scan step size"))
        n = len(vals)
        if not vals or start is None or step is None:
            if "SPECTRUM" in o.get("Scan type", "").upper():
                self.warnings.append(
                    f"Object {o['_object']}: no spectrum data.")
            return                       # snapshots, counters, positions, ...
        native = [start + i * step for i in range(n)]
        anode = self._anode(o)
        hv = anode[1] if anode else None
        label = o.get("Abscissa label", "Kinetic Energy")
        if "kinetic" in label.lower() and hv:
            energy, e_label = [hv - ke for ke in native], "Binding Energy"
        else:
            energy, e_label = native, label
            # "refer to none" is the instrument's own statement that the scan
            # has no X-ray reference (a transmission test, say), so the
            # kinetic-energy axis is what was recorded, not something missing
            if ("kinetic" in label.lower() and "REFER_TO_NONE"
                    not in o.get("Xray Reference Energy", "").upper()):
                self._no_hv.append(o["_object"])
        name = guess_region_name(canon_region_name(
            f"{o.get('Chemical symbol or formula', '')} "
            f"{o.get('Transition or charge state', '')}".strip()
            or o["_object"]), energy)
        reg = Region(
            name=name, index=idx, offset=idx, technique="XPS",
            energy=energy, counts=vals, energy_label=e_label, energy_units="eV",
            count_label="Intensity",
            count_units=o.get("Ordinate units", "counts").lower(),
            decodable=True, sample=o.get("Acquisition name", o["_object"]),
            photon_energy=hv, pass_energy=_num(o.get("Pass energy")),
            dwell=_num(o.get("Dwell time")), step=abs(step) or None,
            lens_mode=o.get("HSA Lens Mode", "").replace("F_HSA_LENS_",
                                                          "").title(),
            date=_date(o.get("Date Acquired", "")),
            anode=anode[0] if anode else "")
        self._decorate(reg, o)
        reg.extra["fields"] = o
        self.regions.append(reg)

    def _decorate(self, reg, o):
        """What a spectrum and a map share: source power, analyser mode,
        aperture, neutraliser, transmission function, sweeps."""
        # older dumps name the gun current / voltage plainly, newer ones (NICPU
        # electronics) hold them under "NICPU X-ray Gun ..." as '0.012 A' /
        # '12000 V'; without the second pair no power was ever found
        cur = _num(o.get("Xray Gun current")
                   or o.get("NICPU X-ray Gun Emission Current"))
        volt = _num(o.get("Xray Gun voltage")
                    or o.get("NICPU X-ray Gun Anode HT Voltage"))
        if cur and volt:
            reg.conditions["X-ray Power"] = f"{cur * volt:g} W"
        if volt:
            reg.conditions["Anode voltage (kV)"] = f"{volt / 1000:g}"
        if cur:
            reg.conditions["Emission current (mA)"] = f"{cur * 1000:.3g}"
        mode = analyser_mode_name(
            o.get("Analyser Scan Mode", "").replace("F_", "", 1))
        if mode:
            reg.extra["analyser_mode"] = mode
        reg.aperture = self._aperture(o)
        neut = self._neutraliser(o)
        if neut:
            reg.extra["neutraliser"] = neut
        tk_, tv = (_list(o.get("Transmission Function Kinetic Energy", "")),
                   _list(o.get("Transmission Function Value", "")))
        if tk_ and len(tk_) == len(tv):
            reg.tf_ke, reg.tf_values = tk_, tv
        sw = o.get("# Sweeps completed")
        if sw:
            reg.conditions["Sweeps"] = sw
            if sw.strip().isdigit() and int(sw) > 0:
                reg.extra["n_scans"] = int(sw)

    def _add_map(self, idx, o):
        """A stigmatic imaging map: one single-energy image. It has no
        spectrum, so the region is not plottable (``decodable`` False, like a
        file with no data) and the pixels ride along as ``extra["cube"]`` (see
        ``kratosmap``); the viewer opens from the tree, as for a SnapMap."""
        name = o["_object"]
        vals = _list(o.get("Ordinate values", ""))
        try:
            cube = kratosmap.build_cube(o, vals)
        except ValueError as exc:
            self.warnings.append(f"Map {name}: not read ({exc}).")
            return
        anode = self._anode(o)
        hv = anode[1] if anode else None
        ke = cube.energy[0]
        if hv:
            cube.energy = [round(hv - ke, 6)]
        elif "REFER_TO_NONE" not in o.get("Xray Reference Energy", "").upper():
            self._no_hv_maps.append(name)
        lens = (o.get("MHSA Lens Mode", "").replace("F_MHSA_", "")
                .replace("_MAGN", "_MAGNIFICATION").replace("_", " ")
                .capitalize())
        z = _num(o.get("Stage Z Position"))
        reg = Region(
            name=canon_region_name(
                f"{o.get('Chemical symbol or formula', '')} "
                f"{o.get('Transition or charge state', '')}".strip() or name),
            index=idx, offset=idx, technique="XPS imaging",
            energy_label="Binding Energy" if hv else "Kinetic Energy",
            energy_units="eV", count_units="counts", decodable=False,
            note=f"Stigmatic image ({cube.nx} x {cube.ny} pixels) at one "
                 f"energy: open it to see the map. {kratosmap.SCALE_NOTE}",
            sample=self._position or "Maps",
            photon_energy=hv, pass_energy=_num(o.get("Pass energy")),
            dwell=_num(o.get("Dwell time")), lens_mode=lens,
            date=_date(o.get("Date Acquired", "")),
            anode=anode[0] if anode else "")
        self._decorate(reg, o)
        if cube.stage_x_mm is not None:
            reg.pos_x, reg.pos_y = cube.stage_x_mm, cube.stage_y_mm
        reg.extra.update({
            "cube": cube, "map_ke": ke, "acq_mode": "Stigmatic map",
            "position_name": self._position,
            "stage_z_um": kratosmap._six(z * 1e6) if z is not None else None,
            "fields": {k: v for k, v in o.items() if k != "Ordinate values"}})
        self.regions.append(reg)
