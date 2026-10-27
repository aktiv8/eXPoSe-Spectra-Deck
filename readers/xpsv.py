"""XPSView / HarwellXPS sessions (``.xpsv``): gzip-compressed JSON.

Two schemas have been seen. Schema 1 is the bare experiment model of XPSView
(``{"model_version": 1, "samples": [...], "spectra": [...], ...}``); schema 2,
written by "XPSView HarwellXPS", wraps it: ``{"schema": 2, "experiment": {...},
"project": ..., "mapping": ..., "content": ..., "calibration": ...,
"provenance": ...}``. Only the experiment is read: the spectra with their
acquisition metadata, the samples and positions, the camera pictures and the
CasaXPS fits. The report-authoring parts (research questions, figures, analyses,
review states, the project record) are not spectra and are left alone.

The layout was worked out from three real sessions and checked against the
instrument files they were made from (``XPS_XPSV_CORPUS``, tests/test_xpsv.py):
every count of every level equals the original's. What differs is the axis.

* A spectrum is ``ke0`` / ``dke`` / ``n`` plus the photon energy ``hv`` and one
  row of raw counts per level (``y``); BE = hv - KE.
* **A CasaXPS charge correction is baked into the stored axis.** When a fit is
  pulled back from CasaXPS the spectrum's ``ke0`` is lowered by the correction
  (``meta.casa_shift``, the sum of the ``Calib`` lines) and ``calib_shift`` stays
  0. The raw axis is ``KE_raw = ke0 + casa_shift``; this reader restores it and
  keeps the correction as ``Region.calibration_shift``, the way a CasaXPS VAMAS
  file is read, so a fit's positions (still in CasaXPS's frame) line up.
* The stored ``hv`` is not always the photon energy: XPSView lets the user edit
  it, and three spectra of one session carry 1484.69 / 1485.69 eV beside
  1486.69 for the same monochromatic source. The value the other spectra of the
  same source agree on is used and the spectrum is named in a warning.
* The Casa lines of the block are kept verbatim in
  ``meta.vamas_block.comment``; ``casafit.parse`` reads the fit from them, so
  nothing is reconstructed. A spectrum with components but no such lines has its
  fit left out (and says so).
* Sample names are not unique (a spot sequence, its repeats and a depth profile
  of one sample are all "ATH PR034"); a distinguishing suffix is added.
"""

from __future__ import annotations

import base64
import collections
import gzip
import json
import re
import zlib

import casafit

from .base import (Region, ImageBlob, SpectrumFile, canon_region_name,
                   clean_text, read_bytes)
from . import vamas

PATTERNS = ("*.xpsv",)
KNOWN_SCHEMAS = (1, 2)
_TECHNIQUES = {"XPS", "UPS", "AES", "ISS", "REELS"}


def _inflate_head(head: bytes, limit: int = 2000) -> str:
    """The first characters of the JSON inside a (possibly truncated) gzip
    stream, or of plain JSON text; "" when it is neither."""
    if head[:2] == b"\x1f\x8b":
        try:
            raw = zlib.decompressobj(31).decompress(head, limit)
        except zlib.error:
            return ""
    else:
        raw = head[:limit]
    return raw.decode("utf-8", "replace")


def sniff(head: bytes, ext: str) -> bool:
    """A gzip (or plain) JSON file that opens like an XPSView session."""
    text = _inflate_head(head).lstrip("﻿ \r\n\t")
    if not text.startswith("{"):
        return False
    if re.match(r'\{\s*"schema"\s*:\s*\d+\s*,\s*"experiment"\s*:', text):
        return True
    return ext == ".xpsv" or ('"model_version"' in text[:600]
                              and '"samples"' in text[:3000])


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f == f and abs(f) != float("inf") else None


def _unwrap(doc):
    """(experiment dict, schema number, wrapper dict or {})."""
    if isinstance(doc, dict) and isinstance(doc.get("experiment"), dict):
        return doc["experiment"], doc.get("schema"), doc
    return doc, 1, {}


def read_session(path: str):
    """(experiment, schema, wrapper) of an ``.xpsv`` file."""
    raw = read_bytes(path)
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return _unwrap(json.loads(raw.decode("utf-8-sig")))


def _watts(text):
    """'225.00W' -> '225 W' (the form the other readers use)."""
    v = _num(re.sub(r"[^0-9.+\-eE]", "", str(text or "")))
    return f"{v:g} W" if v else ""


def _neutraliser(acq):
    if acq.get("neutraliser_on") is None:
        return ""
    if not acq["neutraliser_on"]:
        return "off"
    n = acq.get("neutraliser") or {}
    parts = []
    for key, label, unit in (("filament_current_A", "filament current", "A"),
                             ("charge_balance_V", "charge balance", "V"),
                             ("filament_bias_V", "bias", "V")):
        if _num(n.get(key)) is not None:
            parts.append(f"{label} {n[key]:g} {unit}")
    return "on" + (f" ({', '.join(parts)})" if parts else "")


def _stage(acq, sample):
    """(x_mm, y_mm, z_mm) of an analysis: the spectrum's own stage record,
    else the sample's position; a record in nm (Avantage) is turned into mm."""
    for src in (acq.get("stage"), (sample or {}).get("position")):
        if not isinstance(src, dict):
            continue
        x, y, z = (_num(src.get(k)) for k in ("x", "y", "z"))
        if x is None or y is None:
            continue
        scale = {"nm": 1e-6, "um": 1e-3, "µm": 1e-3, "m": 1e3}.get(
            str(src.get("unit", "mm")).lower(), 1.0)
        return x * scale, y * scale, (None if z is None else z * scale)
    return None


def _anode(acq):
    a, src = clean_text(acq.get("anode")), clean_text(acq.get("source"))
    if not a:
        return ""
    if "mono" in src.lower() and "non" not in src.lower():
        return f"{a} mono"
    return a


def _modal_photon_energies(spectra):
    """{(technique, anode, source): the photon energy most of its spectra
    store}, for the sources of binding-energy spectra."""
    seen = collections.defaultdict(collections.Counter)
    for sp in spectra:
        if (sp.get("technique") or "XPS") not in ("XPS", "UPS"):
            continue
        acq = sp.get("acq") or {}
        hv = _num(sp.get("hv"))
        if hv is not None and acq.get("anode"):
            seen[(sp.get("technique"), acq.get("anode"),
                  acq.get("source"))][round(hv, 4)] += 1
    return {k: c.most_common(1)[0][0] for k, c in seen.items()}


class XpsvFile(SpectrumFile):
    format_name = "XPSView session (.xpsv)"

    def __init__(self):
        super().__init__()
        self._etch = set()

    def load(self, path: str):
        self.path = path
        try:
            exp, schema, _wrap = read_session(path)
        except (OSError, ValueError, EOFError) as exc:
            raise ValueError(f"Not an XPSView session: {exc}") from exc
        if not isinstance(exp, dict) or "spectra" not in exp:
            raise ValueError("Not an XPSView session: no spectra in it.")
        if schema not in KNOWN_SCHEMAS:
            self.warnings.append(
                f"Session schema {schema} is newer than the ones this reader "
                "was written for (1 and 2); it was read as schema 2.")
        samples = {s.get("id"): s for s in exp.get("samples") or []}
        names = self._sample_names(exp, samples)
        self._instrument(exp)
        modal = _modal_photon_energies(exp.get("spectra") or [])
        for sp in exp.get("spectra") or []:
            self._add_spectrum(sp, samples.get(sp.get("sample_id")),
                               names, modal, exp)
        self._images(exp)
        self._ion_gun()
        self._profile()
        self.warnings.append(
            "Only the spectra, fits, positions and pictures of the session "
            "were read; its report text, figures, analyses and project record "
            "are not spectra and were left out.")
        return self._finish()

    # -- names ----------------------------------------------------------
    @staticmethod
    def _sample_names(exp, samples):
        """{sample id: unique label}. A name used by several samples gets the
        first of: its group, its location, its run, its first start time, then
        a number, that makes the labels differ."""
        first_spec = {}
        for sp in exp.get("spectra") or []:
            first_spec.setdefault(sp.get("sample_id"), sp)
        base = {sid: clean_text(s.get("display_name") or s.get("name"))
                or "(unnamed)" for sid, s in samples.items()}
        clash = collections.defaultdict(list)
        for sid, b in base.items():
            clash[b].append(sid)

        def run_label(sid, name):
            run = ((first_spec.get(sid) or {}).get("meta") or {}).get("run")
            lab = clean_text(run.get("label") if isinstance(run, dict) else "")
            toks = []
            for tok in re.split(r"\s+[-–—]\s+", lab):
                tok = re.sub(r"^[\s,;:]+|[\s,;:]+$", "", tok.replace(name, ""))
                if tok:
                    toks.append(tok)
            return " - ".join(toks)

        def start(sid):
            acq = (first_spec.get(sid) or {}).get("acq") or {}
            return clean_text(acq.get("acq_start"))[:16]

        out = {}
        for b, sids in clash.items():
            if len(sids) == 1:
                out[sids[0]] = b
                continue
            chosen = None
            for pick in (lambda s: clean_text(samples[s].get("group")),
                         lambda s: clean_text(
                             (samples[s].get("meta") or {}).get("location")),
                         lambda s: run_label(s, b),
                         start):
                tags = [pick(s) for s in sids]
                if len(set(tags)) == len(sids) and all(tags):
                    chosen = tags
                    break
            if chosen is None:        # a partial label plus a running number
                tags = []
                seen = collections.Counter()
                for s in sids:
                    t = (clean_text(samples[s].get("group"))
                         or run_label(s, b) or start(s))
                    seen[t] += 1
                    tags.append(t if seen[t] == 1 and t else
                                f"{t} #{seen[t]}".strip())
                chosen = [t or f"#{i + 1}" for i, t in enumerate(tags)]
            for s, t in zip(sids, chosen):
                out[s] = f"{b} - {t}"
        # the labels must be unique whatever the data said
        used = collections.Counter(out.values())
        count = collections.Counter()
        for sid, lab in list(out.items()):
            if used[lab] > 1:
                count[lab] += 1
                out[sid] = f"{lab} #{count[lab]}"
        return out

    def _instrument(self, exp):
        meta = exp.get("meta") or {}
        ins = {"Instrument": clean_text(exp.get("instrument")),
               "Operator": clean_text(exp.get("operator")),
               "Project": clean_text(exp.get("project")),
               "Experiment": clean_text(exp.get("name")),
               "Acquisition computer": clean_text(meta.get("machine"))}
        for sp in exp.get("spectra") or []:
            src = clean_text((sp.get("acq") or {}).get("source"))
            if src:
                ins["X-ray source"] = src
                break
        if clean_text(exp.get("notes")):
            ins["Comments"] = exp["notes"].strip()
        self.instrument = {k: v for k, v in ins.items() if v}

    # -- spectra --------------------------------------------------------
    def _add_spectrum(self, sp, sample, names, modal, exp):
        import numpy as np
        tech = (sp.get("technique") or "XPS").upper()
        tech = tech if tech in _TECHNIQUES else "XPS"
        rows = sp.get("y") or []
        n = int(sp.get("n") or 0)
        if not rows or not n:
            self.warnings.append(
                f"{sp.get('name')}: no data in the session; left out.")
            return
        acq = sp.get("acq") or {}
        meta = sp.get("meta") or {}
        sid = sp.get("sample_id")
        sample_label = names.get(sid, clean_text((sample or {}).get("name")))
        label = sp.get("display_name") or sp.get("name") or ""
        kinetic = tech in ("ISS", "REELS") or (
            sp.get("energy_scale") == "KE" and _num(sp.get("hv")) is None)

        hv_stored = _num(sp.get("hv"))
        hv = hv_stored
        shift = _num(meta.get("casa_shift")) or 0.0
        if not kinetic and hv_stored is not None and acq.get("anode"):
            key = (sp.get("technique"), acq.get("anode"), acq.get("source"))
            m = modal.get(key)
            if m is not None and abs(m - hv_stored) > 1e-3:
                hv = m
                self.warnings.append(
                    f"{label} ({sample_label}): XPSView stored a photon energy "
                    f"of {hv_stored:g} eV but the other spectra of the same "
                    f"source say {m:g} eV; {m:g} eV was used.")
        user_shift = _num(sp.get("calib_shift")) or 0.0
        if abs(user_shift) > 1e-9 and not kinetic:
            self.warnings.append(
                f"{label} ({sample_label}): XPSView holds a binding-energy "
                f"shift of {user_shift:+g} eV set by hand; it is not applied "
                "(set a shift on the spectrum here if you want it).")

        ke0, dke = _num(sp.get("ke0")), _num(sp.get("dke"))
        if ke0 is None or not dke:
            self.warnings.append(f"{label}: no energy axis; left out.")
            return
        ke = ke0 + dke * np.arange(n)
        if kinetic:
            energy = [round(float(v), 6) for v in ke]
            e_label = "Kinetic Energy"
            photon = None
        else:
            if hv is None:
                energy, e_label, photon = ([round(float(v) + shift, 6)
                                            for v in ke], "Kinetic Energy",
                                           None)
                self.warnings.append(
                    f"{label} ({sample_label}): no photon energy in the "
                    "session; the kinetic-energy axis is kept.")
            else:
                energy = [round(float(hv - (k + shift)), 6) for k in ke]
                e_label, photon = "Binding Energy", hv
        step = abs(round(dke, 6))
        levels = sp.get("levels") or [{}]
        if len(levels) != len(rows):
            levels = [(levels[i] if i < len(levels) else {})
                      for i in range(len(rows))]
        cps = bool(sp.get("y_is_cps"))
        pos = _stage(acq, sample)
        comment = list((meta.get("vamas_block") or {}).get("comment") or [])
        has_fit_lines = any(str(l).strip().startswith("CASA ")
                            for l in comment)
        if (sp.get("components") and not has_fit_lines
                and len(rows) == 1):
            self.warnings.append(
                f"{label} ({sample_label}): the session has fitted components "
                "but not the CasaXPS lines they came from; the fit was left "
                "out.")
        multi_fit = (sp.get("components") and len(rows) > 1)
        if multi_fit:
            self.warnings.append(
                f"{label} ({sample_label}): fits of a multi-level spectrum "
                "are not read from a session.")

        for li, row in enumerate(rows):
            y = [None if v is None else float(v) for v in row]
            if any(v is None for v in y):
                y = [float("nan") if v is None else v for v in y]
            lv = levels[li] if isinstance(levels[li], dict) else {}
            r = Region(
                name=canon_region_name(label), index=len(self.regions),
                offset=len(self.regions), technique=tech,
                energy=energy, counts=y,
                energy_label=e_label, energy_units="eV",
                count_label="Intensity",
                count_units="counts/s" if cps else "counts",
                decodable=True, sample=sample_label,
                photon_energy=photon,
                pass_energy=_num(sp.get("pass_energy")) or None,
                dwell=_num(sp.get("dwell")) or None, step=step or None,
                lens_mode=clean_text(acq.get("lens_mode")),
                aperture=clean_text(acq.get("aperture")),
                anode="" if kinetic else _anode(acq))
            r.extra["xpsv"] = {"spectrum": sp.get("id"), "sample": sid,
                               "level": li, "name": sp.get("name"),
                               "casa_shift_baked": shift,
                               "hv_stored": hv_stored}
            if abs(user_shift) > 1e-9:
                r.extra["xpsv"]["user_shift"] = user_shift
            if sp.get("notes"):
                r.extra["xpsv"]["notes"] = sp["notes"]
            self._acquisition(r, sp, acq, lv, exp, pos)
            et = _num(lv.get("etch_time"))
            if et is not None and len(rows) > 1:
                r.etch_time = et
            if lv.get("iteration") is not None:
                r.extra["iteration"] = lv["iteration"]
            if has_fit_lines and len(rows) == 1:
                self._fit(r, comment)
            if tech in ("XPS", "UPS") and not kinetic:
                self._calibration(r, comment, shift, lv)
            self._transmission(r, sp)
            self.regions.append(r)

    def _calibration(self, r, comment, shift, level):
        """Charge correction: CasaXPS's own Calib lines when the block has
        them (so the fit frame is right), else the stored sum."""
        before = r.calibration_shift
        vamas.VamasFile._own_calibration(r, comment)
        if r.calibration_shift == before and abs(shift) > 1e-9:
            r.calibration_shift = shift
            r.extra["casa_calib"] = {"measured": 0.0, "assigned": shift,
                                     "shift": shift, "inherited": False}
        got = r.calibration_shift
        if abs(got - shift) > 1e-3 and abs(shift) > 1e-9:
            self.warnings.append(
                f"{r.name} ({r.sample}): the Calib lines add up to {got:+.3f} "
                f"eV but the session stored a shift of {shift:+.3f} eV; the "
                "axis was restored with the stored shift.")
            r.calibration_shift = shift
        if abs(shift) > 1e-9:
            r.extra["casa_calib"]["shift"] = r.calibration_shift

    def _fit(self, r, comment):
        fit = casafit.parse(comment)
        if fit is not None:
            r.fit = fit
        r.extra["comment_lines"] = [str(l) for l in comment
                                    if not str(l).startswith("===")]

    @staticmethod
    def _transmission(r, sp):
        tf = sp.get("tf") or []
        pairs = [(float(a), float(b)) for a, b in
                 (p for p in tf if isinstance(p, (list, tuple)) and len(p) == 2
                  and _num(p[0]) is not None and _num(p[1]) is not None)]
        if len(pairs) >= 2:
            pairs.sort()
            r.tf_ke = [a for a, _ in pairs]
            r.tf_values = [b for _, b in pairs]

    def _acquisition(self, r, sp, acq, level, exp, pos):
        vendor = (exp.get("vendor") or "").lower()
        values = level.get("values") if isinstance(level.get("values"),
                                                   dict) else {}
        r.extra["n_scans"] = int(sp.get("sweeps") or acq.get("sweeps") or 1)
        # XPSView keeps the dwell of one sweep and the sum beside it
        # (acq.dwell_total_s); ESCApe, like our Kratos reader, holds the sum
        total = _num(acq.get("dwell_total_s"))
        if vendor == "kratos" and total:
            r.dwell = total
            r.extra["dwell_total"] = True
        mode = clean_text(acq.get("analyser_mode"))
        if mode:
            from .base import analyser_mode_name
            r.extra["analyser_mode"] = analyser_mode_name(mode)
        power = _watts(acq.get("xray_power")
                       or (f"{acq['xray_power_W']}W"
                           if _num(acq.get("xray_power_W")) else ""))
        if power:
            r.conditions["X-ray Power"] = power
        if _num(acq.get("emission_current_mA")):
            r.conditions["Emission current (mA)"] = \
                f"{acq['emission_current_mA']:g}"
        if _num(acq.get("anode_voltage_V")):
            r.conditions["Anode voltage (kV)"] = \
                f"{acq['anode_voltage_V'] / 1000:g}"
        spot = acq.get("spot_um")
        if isinstance(spot, (list, tuple)) and spot and _num(spot[0]):
            w, h = _num(spot[0]), _num(spot[-1])
            r.conditions["X-ray spot (µm)"] = (
                f"{w:g}" if h is None or abs(w - h) < 1e-9
                else f"{w:g} × {h:g}")
        neut = _neutraliser(acq)
        if neut:
            r.extra["neutraliser"] = neut
        start = (clean_text(values.get("acquired"))
                 or clean_text(acq.get("acq_start")))
        if start:
            r.date = start
            r.extra["t_start"] = start
        end = clean_text(acq.get("acq_end"))
        if end and not values.get("acquired"):
            r.extra["t_end"] = end
        if pos:
            r.pos_x, r.pos_y = pos[0], pos[1]
            if pos[2] is not None:
                r.extra["pos_z"] = pos[2]
        run = (sp.get("meta") or {}).get("run")
        scan = clean_text(run.get("label")) if isinstance(run, dict) else ""
        if scan:
            r.extra["scan"] = scan
        if clean_text(acq.get("etch_ion_beam")):
            self._etch.add(clean_text(acq["etch_ion_beam"]))
        bt = (sp.get("meta") or {}).get("block_type")
        rows = {}
        if bt:
            why = clean_text((sp.get("meta") or {}).get("block_type_reason"))
            rows["XPSView block type"] = bt + (f" ({why})" if why else "")
        for key, label in (("analysis_area", "Analysis area"),
                           ("flood_gun", "Flood gun"),
                           ("detector_mode", "Detector mode"),
                           ("system_number", "System number"),
                           ("calibration_energy_scale",
                            "Energy-scale calibration")):
            if clean_text(acq.get(key)):
                rows[label] = clean_text(acq[key])
        for key, label, fmt in (("work_function_eV", "Work function (eV)",
                                 "{:.2f}"),
                                ("source_energy_eV",
                                 "Primary beam energy (eV)", "{:g}"),
                                ("flood_gun_current_uA",
                                 "Flood gun current (µA)", "{:g}"),
                                ("flood_gun_energy_eV",
                                 "Flood gun energy (eV)", "{:g}")):
            if _num(acq.get(key)) is not None:
                rows[label] = fmt.format(acq[key])
        if acq.get("transmission_corrected") is not None:
            rows["Transmission corrected"] = (
                "yes" if acq["transmission_corrected"] else "no")
        if rows:
            r.extra["acq_metadata"] = rows

    # -- pictures, ion gun, depth profile ---------------------------------
    def _images(self, exp):
        samples = {s.get("id"): s for s in exp.get("samples") or []}
        for im in exp.get("images") or []:
            try:
                data = base64.b64decode(im.get("data") or "")
            except (ValueError, TypeError):
                data = b""
            if not data:
                continue
            fmt = "png" if "png" in str(im.get("mime", "")).lower() else "jpeg"
            note = ("Picture from an XPSView session; XPSView's own "
                    "calibration of it is not used.")
            sid = im.get("sample_id")
            self.images.append(ImageBlob(
                clean_text(im.get("name")) or "Picture", len(self.images),
                data, True, note, fmt=fmt,
                sample=clean_text((samples.get(sid) or {}).get("name")),
                taken=clean_text(im.get("acquired"))))

    def _ion_gun(self):
        """The etch beam, when every sputtered spectrum names the same one."""
        if len(self._etch) == 1:
            self.instrument["Ion gun / sputtering"] = next(iter(self._etch))

    def _profile(self):
        """A depth profile when some spectrum has etch times for its levels:
        every region then gets a level number (as a VAMAS profile does)."""
        times = sorted({r.etch_time for r in self.regions
                        if r.etch_time is not None})
        if not times:
            return
        levels, names = {}, []
        for r in self.regions:
            if r.name not in names:
                names.append(r.name)
            k = (r.sample, r.name)
            levels[k] = levels.get(k, -1) + 1
            r.etch_level = levels[k]
        self.depth_profile = {
            "is_profile": True, "n_levels": max(levels.values()) + 1,
            "regions_per_level": len(names),
            "etch_per_level": (times[1] - times[0]) if len(times) > 1 else 0.0,
            "total_etch_time": times[-1], "cumulative": times,
            "etch_source": "etch times stored in the XPSView session"}
