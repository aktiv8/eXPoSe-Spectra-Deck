"""NeXus ``.nxs`` reader for analyser (NXxps / NXmpes) spectra, built on the
NEXAFS reader so one ``.nxs`` file may hold both.

**What it was checked against, and what it was not.** No XPS ``.nxs`` file
from an instrument or beamline is available here (the 200 on this machine are
all Diamond B07 NEXAFS scans), so this reader is checked on the files this
application writes (``nexus_export.py``: a round trip of every field it
writes, including real Kratos, Avantage and CasaXPS data) and follows the
NeXus application definitions' own names. A file from another converter
(``pynxtools-xps`` and the like) is read through the same names (class, not
group name: ``NXelectronanalyzer``, ``NXenergydispersion``, ``NXbeam`` ...),
but whether a given beamline's layout fits has not been seen, so a spectrum
that cannot be read is left out with a note instead of being guessed.

One ``NXentry`` gives one ``Region``, taken from its ``NXdata`` (``@signal``,
the first of ``@axes`` or an ``energy`` dataset). Only a one-dimensional
signal is read (an angle- or position-resolved cube is left out with a note).
The energy axis follows its ``@type``: ``binding`` is kept, ``kinetic`` is
turned into a binding energy when the photon energy is known (as every reader
here does) and otherwise stays kinetic; an axis with no type keeps its own
name, never assumed. When the entry carries the ``metadata`` note this
application writes (the whole ``region_metadata`` as JSON), what VAMAS can
restore is restored from it (``vamasmeta.apply_to_region``) and the rest is
kept as preserved metadata. ``NXfit`` groups are not read back: the fit is
not rebuilt from NeXus (a note says how many were left).
"""

from __future__ import annotations

import json

from . import nexus_nexafs as nx
from .base import Region, canon_region_name

#: NXxps ``energy_scan_mode`` -> the analyser-mode words used here
_SCAN_MODES = {"fixed_analyzer_transmission": "Constant analyser energy (CAE)",
               "fixed_retarding_ratio": "Constant retard ratio (CRR)"}
_DEFINITIONS = ("nxxps", "nxmpes")
_SKIP_GROUPS = ("metadata",)


def _children(group):
    h5py = nx._h5py()
    return [(n, g) for n, g in group.items() if isinstance(g, h5py.Group)]


def _find(group, nx_class, depth=4):
    """The first descendant of ``group`` whose ``NX_class`` is ``nx_class``
    (breadth first, so the nearest wins), or None."""
    level = [group]
    for _ in range(depth):
        nxt = []
        for g in level:
            for n, c in _children(g):
                if nx._attr(c, "NX_class") == nx_class:
                    return c
                if not n.startswith(_SKIP_GROUPS):
                    nxt.append(c)
        level = nxt
    return None


def _dataset(group, name, depth=3):
    """Dataset ``name`` in ``group`` or, failing that, the nearest group below
    it (an analyser's ``pass_energy`` is in its ``NXenergydispersion``)."""
    h5py = nx._h5py()
    if group is None:
        return None
    level = [group]
    for _ in range(depth):
        nxt = []
        for g in level:
            ds = g.get(name)
            if isinstance(ds, h5py.Dataset):
                return ds
            nxt += [c for n, c in _children(g) if not n.startswith("data")]
        level = nxt
    return None


def _num(group, name):
    ds = _dataset(group, name)
    return None if ds is None else nx._number(ds[()])


def _str(group, name):
    ds = _dataset(group, name)
    return "" if ds is None else nx._text(ds[()])


def find_data(entry):
    """The ``NXdata`` holding the spectrum: ``data`` or the first with a
    ``@signal``."""
    cand = entry.get("data")
    if cand is not None and nx._attr(cand, "NX_class") == "NXdata":
        return cand
    for _n, c in _children(entry):
        if nx._attr(c, "NX_class") == "NXdata" and nx._attr(c, "signal"):
            return c
    return None


def is_xps_entry(entry):
    """An entry that declares NXxps / NXmpes, or has an analyser, and has a
    spectrum to read."""
    definition = nx._string(entry, "definition").lower()
    ins = entry.get("instrument")
    analyser = ins is not None and _find(ins, "NXelectronanalyzer") is not None
    return (definition.startswith(_DEFINITIONS) or analyser) \
        and find_data(entry) is not None


class NexusFile(nx.NexusNexafsFile):
    """A ``.nxs`` holding NEXAFS scans, analyser spectra, or both."""

    format_name = "NeXus (.nxs)"
    _what = "NXxps spectrum or photon-energy scan"

    def _accepts(self, entry):
        return is_xps_entry(entry)

    def _after_entries(self, entries, used):
        kinds = set()
        for _n, g in used:
            kinds.add("scan" if nx._is_scan(g) else "xps")
        self.format_name = {"scan": "NeXus NEXAFS (.nxs)",
                            "xps": "NeXus NXxps (.nxs)"}.get(
            next(iter(kinds)) if len(kinds) == 1 else "", "NeXus (.nxs)")
        left = len(entries) - len(used)
        if left:
            self.warnings.append(f"{left} NeXus entr{'y' if left == 1 else 'ies'}"
                                 " that are neither an analyser spectrum nor "
                                 "a NEXAFS scan were left out.")

    # -- one analyser spectrum ----------------------------------------------
    def _add_other(self, name, entry, stem):
        data = find_data(entry)
        sig_name = nx._attr(data, "signal") or "data"
        if sig_name not in data:
            self.warnings.append(f"{name}: its NXdata names no signal.")
            return
        sig = data[sig_name]
        if sig.ndim != 1:
            self.warnings.append(
                f"{name}: the signal has {sig.ndim} dimensions (angle- or "
                "position-resolved): only a one-dimensional spectrum is read.")
            return
        axis_name = nx._attr(data, "axes")
        if axis_name not in data:
            axis_name = next((k for k in ("energy", "binding_energy",
                                          "kinetic_energy") if k in data), "")
        if not axis_name:
            self.warnings.append(f"{name}: no energy axis in its NXdata.")
            return
        axis = data[axis_name]
        ys = [float(v) for v in sig[()]]
        xs = [float(v) for v in axis[()]]
        if axis.ndim != 1 or len(xs) != len(ys):
            self.warnings.append(
                f"{name}: the energy axis has {len(xs)} values for "
                f"{len(ys)} intensities: left out.")
            return

        ins = entry.get("instrument")
        ana = _find(ins, "NXelectronanalyzer") if ins is not None else None
        beam = _find(ins, "NXbeam") if ins is not None else None
        src = _find(ins, "NXsource") if ins is not None else None
        smp = _find(entry, "NXsample")
        md = self._note(entry)
        region_md = md.get("region") if isinstance(md.get("region"), dict) \
            else {}

        hv = _num(beam, "incident_energy")
        kind = (nx._attr(axis, "type") or axis_name).lower()
        label = "Energy"
        if "binding" in kind:
            label = "Binding Energy"
        elif "kinetic" in kind:
            if hv:
                xs, label = [hv - x for x in xs], "Binding Energy"
            else:
                label = "Kinetic Energy"
        units = nx._attr(axis, "units") or "eV"

        # a note of ours knows the sample even when it is empty; otherwise
        # the file name stands in (a NeXus entry has no sample-less form)
        sample = _str(smp, "name") or (
            str(region_md.get("Sample", "")) if region_md else stem)
        title = nx._string(entry, "title")
        rname = region_md.get("Region") or (
            title[len(sample):].strip() if sample and title.startswith(sample)
            and title[len(sample):].strip() else title or name)
        wall0, utc0 = nx.parse_time(nx._string(entry, "start_time"))
        wall1, utc1 = nx.parse_time(nx._string(entry, "end_time"))

        r = Region(
            name=canon_region_name(rname), index=len(self.regions),
            offset=len(self.regions), technique="XPS",
            energy=xs, counts=ys, energy_label=label, energy_units=units,
            count_label=nx._attr(sig, "long_name") or "Intensity",
            count_units=nx._attr(sig, "units") or "arb.", decodable=True,
            sample=sample, photon_energy=hv or None,
            pass_energy=_num(ana, "pass_energy"),
            dwell=_num(ana, "dwell_time"), step=_num(ana, "energy_step"),
            lens_mode=_str(ana, "lens_mode"), anode=_str(src, "name"),
            date=wall0)
        scans = _num(ana, "number_of_scans")
        if scans:
            r.extra["n_scans"] = int(scans)
        mode = _SCAN_MODES.get(_str(ana, "energy_scan_mode"), "")
        if mode:
            r.extra["analyser_mode"] = mode
        if (utc0 or wall0) and (utc1 or wall1):    # a run window
            r.extra["t_start"] = utc0 or wall0
            r.extra["t_end"] = utc1 or wall1
            if utc0:
                r.extra["tz"] = "UTC"
        power = _num(src, "power")
        if power:
            r.conditions["X-ray Power"] = f"{power:g} W"
        volts = _num(src, "voltage")
        if volts:
            r.conditions["Anode voltage (kV)"] = f"{volts / 1000.0:g}"
        emission = _num(src, "emission_current")
        if emission:
            r.conditions["Emission current (mA)"] = f"{emission:g}"
        spot = _num(beam, "extent")
        if spot:
            r.conditions["X-ray spot (µm)"] = f"{spot:g}"
        if region_md:               # what VAMAS can hold, as it was ...
            import vamasmeta
            vamasmeta.apply_to_region(r, region_md)
            self._run_window(r, region_md)
            self._dwell_total(r, region_md)
        self._sample_fields(r, smp)      # ... then the exact numbers win
        self._referencing(r, entry, label)
        tf = _find(ana, "NXdata") if ana is not None else None
        if tf is not None and "kinetic_energy" in tf \
                and "relative_intensity" in tf:
            pairs = sorted(zip((float(v) for v in tf["kinetic_energy"][()]),
                               (float(v) for v in
                                tf["relative_intensity"][()])))
            r.tf_ke = [p[0] for p in pairs]     # ascending, as the
            r.tf_values = [p[1] for p in pairs]  # interpolation needs
        r.extra["nexus"] = {"entry": name,
                            "program": nx._string(entry, "program_name"),
                            "definition": nx._string(entry, "definition")}
        if md.get("entered_by_user"):
            r.extra["nexus"]["entered_by_user"] = md["entered_by_user"]
        fits = sum(1 for _n, g in _children(entry)
                   if nx._attr(g, "NX_class") == "NXfit")
        if fits:
            self.warnings.append(
                f"{name}: {fits} NXfit group{'s' if fits != 1 else ''} not "
                "read (a fit is not rebuilt from NeXus).")
        if not self.instrument:
            self.instrument = self._analyser_instrument(entry, ins, ana, md)
        self.regions.append(r)

    @staticmethod
    def _note(entry):
        """The ``metadata`` note this application writes, as a dict ({} for
        anything else or text that is not JSON)."""
        note = entry.get("metadata")
        raw = nx._read(note, "data") if note is not None else None
        if raw is None:
            return {}
        try:
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            doc = json.loads(str(raw))
        except (ValueError, UnicodeDecodeError):
            return {}
        return doc if isinstance(doc, dict) else {}

    @staticmethod
    def _dwell_total(r, region_md):
        """A Kratos ``.experiment`` stores the dwell already summed over the
        sweeps; NeXus holds it per sweep. The note's own "Dwell (s)" says
        which the source used (the same time either way, so only what the
        metadata shows changes)."""
        scans = r.extra.get("n_scans") or 0
        try:
            shown = float(region_md.get("Dwell (s)"))
        except (TypeError, ValueError):
            return
        if (scans > 1 and r.dwell
                and abs(shown - r.dwell * scans) <= 0.005 * shown
                and abs(shown - r.dwell) > 0.005 * shown):
            r.dwell = r.dwell * scans
            r.extra["dwell_total"] = True

    @staticmethod
    def _run_window(r, region_md):
        """A run window (start and finish) the metadata note holds, which the
        entry's own ``start_time`` cannot say when only the start is known."""
        import timing
        t0 = timing.parse_ts(str(region_md.get("Run started", ""))[:19])
        if t0 is None:
            return
        t1 = timing.parse_ts(str(region_md.get("Run finished", ""))[:19])
        r.extra["t_start"] = t0.strftime("%Y-%m-%d %H:%M:%S")
        if t1 is not None:
            r.extra["t_end"] = t1.strftime("%Y-%m-%d %H:%M:%S")
        if str(region_md["Run started"]).rstrip().endswith("UTC"):
            r.extra["tz"] = "UTC"

    @staticmethod
    def _sample_fields(r, smp):
        if smp is None:
            return
        pos = _dataset(smp, "stage_position")
        if pos is not None and pos.shape and pos.shape[0] >= 2:
            x, y = (float(v) for v in pos[()][:2])
            r.pos_x, r.pos_y = x, y
        level = _num(smp, "etch_level")
        if level is not None:
            r.etch_level = int(level)
        t = _num(smp, "etch_time")
        if t is not None:
            r.etch_time = t

    @staticmethod
    def _referencing(r, entry, label):
        """``energy_referencing``: a charge correction the file states.
        ``applied`` means the axis already has it (nothing is added again);
        otherwise it is the correction to add, as ``calibration_shift``."""
        cal = _find(entry, "NXcalibration")
        offset = _num(cal, "offset") if cal is not None else None
        if not offset:
            return
        applied = _dataset(cal, "applied")
        applied = bool(applied[()]) if applied is not None else False
        measured, assigned = _num(cal, "measured"), _num(cal, "assigned")
        r.extra["nexus_referencing"] = {"offset": offset, "applied": applied}
        if applied or label != "Binding Energy":
            return
        r.calibration_shift = offset
        if measured is not None and assigned is not None:
            r.extra["casa_calib"] = {"measured": measured,
                                     "assigned": assigned, "shift": offset,
                                     "inherited": False}

    @staticmethod
    def _analyser_instrument(entry, ins, ana, md):
        user = _find(entry, "NXuser")
        dev = _find(ins, "NXfabrication") if ins is not None else None
        out = {
            "Instrument": _str(dev, "model"),
            "Operator": _str(user, "name"),
            "Institution": _str(user, "affiliation"),
            "Work function (eV)": ("" if _num(ana, "work_function") is None
                                   else f"{_num(ana, 'work_function'):g}"),
        }
        exp = md.get("experiment")
        if isinstance(exp, dict):
            for k, v in exp.items():        # what the file once said about the job
                if v and not out.get(k):
                    out[k] = str(v)
        return {k: v for k, v in out.items() if v}
