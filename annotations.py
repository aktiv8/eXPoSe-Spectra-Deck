"""User edits kept *beside* the data: display names, notes, edited metadata,
binding-energy shifts, peak markers and the calibration statement.

The readers' output is never modified. The app applies these annotations when
it draws, exports and reports, and stores them in the workbook, so they
survive reloading the (unchanged) instrument files. Keys use the workbook's
file id (``f1``...) plus the sample and region names as read from the file:

* sample key  ``"<fid>|<sample>"``
* region key  ``"<fid>|<sample>|<region name>"`` (all levels of that region)
* metadata key ``"<fid>|<position in parser.regions>"`` (one spectrum)

Tk-free and JSON-serialisable.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field

import columntext
import elements
import nexafs
import nexus_settings
import repeats
import sputter as sputter_mod

VERSION = 1


def sample_key(fid, sample):
    return f"{fid}|{sample}"


def region_key(fid, sample, name):
    return f"{fid}|{sample}|{name}"


def meta_key(fid, pos):
    return f"{fid}|{pos}"


def _clean_reels(d):
    """A valid REELS construction ``{elastic, p1, p2, gap, base, slope}``
    (gap, base and slope may be absent) or None."""
    import math
    if not isinstance(d, dict):
        return None
    try:
        out = {"elastic": float(d["elastic"]),
               "p1": [float(d["p1"][0]), float(d["p1"][1])],
               "p2": [float(d["p2"][0]), float(d["p2"][1])]}
    except (KeyError, TypeError, ValueError, IndexError):
        return None
    for k in ("gap", "base", "slope"):
        try:
            v = float(d[k])
            out[k] = v if math.isfinite(v) else None
        except (KeyError, TypeError, ValueError):
            out[k] = None
    if not all(math.isfinite(v) for v in (out["elastic"], *out["p1"],
                                          *out["p2"])):
        return None
    return out


def _clean_plasmon(d):
    """A valid REELS plasmon result (``reels.plasmon`` plus the elastic peak
    it was measured from) or None."""
    import math
    if not isinstance(d, dict):
        return None
    try:
        out = {"elastic": float(d["elastic"]), "ep": float(d["ep"]),
               "window": [float(d["window"][0]), float(d["window"][1])],
               "method": str(d["method"]),
               "background": str(d["background"])}
    except (KeyError, TypeError, ValueError, IndexError):
        return None
    if out["method"] not in ("gauss", "d2") or \
            out["background"] not in ("Linear", "Shirley", "None"):
        return None
    if not all(math.isfinite(v) for v in (out["elastic"], out["ep"],
                                          *out["window"])):
        return None
    for k in ("fwhm", "frac", "smooth"):
        try:
            v = float(d[k])
            out[k] = v if math.isfinite(v) else None
        except (KeyError, TypeError, ValueError):
            out[k] = None
    bg = []
    for p in d.get("bg") or []:
        try:
            a, b = float(p[0]), float(p[1])
        except (TypeError, ValueError, IndexError):
            continue
        if math.isfinite(a) and math.isfinite(b):
            bg.append([a, b])
    out["bg"] = bg
    return out


def _finite(v):
    import math
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _clean_loss_fit(d):
    """A valid REELS loss-region fit (``reels.decompose`` plus the elastic peak
    it was measured from, and optionally the ``geometry`` and ``layers`` of a
    layer estimate) or None."""
    if not isinstance(d, dict):
        return None
    elastic = _finite(d.get("elastic"))
    try:
        win = [float(d["window"][0]), float(d["window"][1])]
    except (KeyError, TypeError, ValueError, IndexError):
        return None
    if elastic is None or _finite(win[0]) is None or _finite(win[1]) is None             or d.get("background") not in ("Linear", "Shirley", "None"):
        return None
    comps = []
    for c in d.get("components") or []:
        try:
            item = {"key": str(c["key"]), "name": str(c["name"]),
                    "loss": float(c["loss"]), "fwhm": float(c["fwhm"]),
                    "area": float(c["area"])}
        except (KeyError, TypeError, ValueError):
            continue
        if all(_finite(item[k]) is not None for k in ("loss", "fwhm", "area"))                 and item["fwhm"] > 0 and item["area"] >= 0:
            comps.append(item)
    if not comps:
        return None
    out = {"elastic": elastic, "window": win, "background": d["background"],
           "tol": _finite(d.get("tol")) or 0.0, "components": comps,
           "r": _finite(d.get("r")), "rms": _finite(d.get("rms")), "bg": []}
    for p in d.get("bg") or []:
        try:
            a, b = float(p[0]), float(p[1])
        except (TypeError, ValueError, IndexError):
            continue
        if _finite(a) is not None and _finite(b) is not None:
            out["bg"].append([a, b])
    g = d.get("geometry")
    if isinstance(g, dict):
        geo = {k: _finite(g.get(k)) for k in
               ("lam", "a_in", "a_out", "d", "r_graphite", "p")}
        if all(geo[k] is not None for k in ("lam", "a_in", "a_out", "d",
                                            "r_graphite", "p")):
            out["geometry"] = geo
            ly = d.get("layers")
            if isinstance(ly, dict):
                out["layers"] = {k: _finite(ly.get(k)) for k in
                                 ("n", "at_least", "id", "n_graphite", "L")}
    return out


def _clean_window(w):
    try:
        lo, hi = float(w[0]), float(w[1])
    except (TypeError, ValueError, IndexError):
        return None
    return [lo, hi] if _finite(lo) is not None and _finite(hi) is not None \
        and hi > lo else None


def _clean_pi_ratio(d):
    """A valid REELS π / (σ+π) ratio result (``reels.pi_ratio`` plus the
    elastic peak it was measured from) or None. The ratio may be None (a σ+π
    area that is not positive) and a π area may be negative: both are kept as
    measured."""
    if not isinstance(d, dict):
        return None
    pw, sw = _clean_window(d.get("pi_window")), _clean_window(d.get("sp_window"))
    elastic, i_pi, i_sp = (_finite(d.get(k)) for k in
                           ("elastic", "i_pi", "i_sp"))
    try:
        avg = int(d["avg"])
    except (KeyError, TypeError, ValueError):
        return None
    if None in (pw, sw, elastic, i_pi, i_sp) or avg < 1:
        return None
    out = {"pi_window": pw, "sp_window": sw, "avg": avg, "elastic": elastic,
           "i_pi": i_pi, "i_sp": i_sp, "ratio": _finite(d.get("ratio"))}
    for k in ("pos_pi", "pos_sp", "dip_pi", "dip_sp"):
        out[k] = _finite(d.get(k))
    for k in ("edge_pi", "edge_sp"):
        out[k] = bool(d.get(k))
    for key in ("bg_pi", "bg_sp"):
        pts = []
        for p in d.get(key) or []:
            try:
                a, b = float(p[0]), float(p[1])
            except (TypeError, ValueError, IndexError):
                continue
            if _finite(a) is not None and _finite(b) is not None:
                pts.append([a, b])
        out[key] = pts
    return out


def _clean_reels_cal(d):
    """A valid REELS energy calibration ``{elastic, nominal, shift, source}``
    or None: the elastic peak as measured, the primary energy it was moved to
    and the shift (kinetic energy added to every point); ``source`` says where
    the primary energy came from ("file" or "entered")."""
    if not isinstance(d, dict):
        return None
    elastic, nominal, shift = (_finite(d.get(k)) for k in
                               ("elastic", "nominal", "shift"))
    if None in (elastic, nominal, shift) or nominal <= 0:
        return None
    return {"elastic": elastic, "nominal": nominal, "shift": shift,
            "source": "file" if d.get("source") == "file" else "entered"}


def _clean_reels_ref(d):
    """A valid reference ratio (HOPG or diamond) or None: a ``ratio``, a
    ``label`` and ``kind`` ("spectrum" = measured here, "entered" = typed),
    with the windows, end-point channels and elastic-peak energy it was measured
    with when it was measured here."""
    if not isinstance(d, dict):
        return None
    ratio = _finite(d.get("ratio"))
    if ratio is None:
        return None
    out = {"ratio": ratio, "label": str(d.get("label") or "reference"),
           "kind": "entered" if d.get("kind") == "entered" else "spectrum",
           "source": str(d.get("source") or ""),
           "pi_window": _clean_window(d.get("pi_window")),
           "sp_window": _clean_window(d.get("sp_window")),
           "e0": _finite(d.get("e0"))}
    try:
        out["avg"] = int(d["avg"]) if d.get("avg") is not None else None
    except (TypeError, ValueError):
        out["avg"] = None
    return out


@dataclass
class Annotations:
    sample_names: dict = field(default_factory=dict)     # sample key -> name
    region_names: dict = field(default_factory=dict)     # region key -> name
    sample_notes: dict = field(default_factory=dict)     # sample key -> text
    region_notes: dict = field(default_factory=dict)     # region key -> text
    md_edits: dict = field(default_factory=dict)         # meta key -> {k: v}
    shifts: dict = field(default_factory=dict)           # scope key -> eV
    calibration: list = field(default_factory=list)      # log of calibrations
    calibration_statement: str = ""                      # optional override
    markers: dict = field(default_factory=dict)          # region key -> [...]
    experiment_notes: str = ""
    sputter: dict = field(default_factory=dict)          # sample key -> settings
    reels: dict = field(default_factory=dict)            # region key -> band gap
    plasmon: dict = field(default_factory=dict)          # region key -> REELS plasmon energy
    loss_fit: dict = field(default_factory=dict)         # region key -> REELS loss-region fit (reels.decompose)
    pi_ratio: dict = field(default_factory=dict)         # region key -> REELS pi / (sigma+pi) ratio (reels.pi_ratio)
    reels_refs: dict = field(default_factory=dict)       # "hopg" / "diamond" -> the reference ratios of this workbook
    reels_cal: dict = field(default_factory=dict)        # region key -> REELS energy calibration (reels.calibration)
    instrument: dict = field(default_factory=dict)       # file id -> NeXus settings
    iss: dict = field(default_factory=dict)              # file id -> {ion, e0, theta}
    imports: dict = field(default_factory=dict)          # file id -> how a column-text file is read
    nexafs_ring: bool = False                            # NEXAFS scaled to the mean ring current
    nexafs_edge: dict = field(default_factory=dict)      # region key -> edge normalisation (nexafs.py)
    repeat_groups: dict = field(default_factory=dict)    # name -> {kind, samples} (repeats.py)
    extra: dict = field(default_factory=dict)            # unknown keys kept

    # -- samples a reader split into several (Kratos runs) ------------------------
    ADOPTED = ("region_names", "sample_notes", "region_notes", "shifts",
               "markers", "sputter", "reels", "plasmon", "loss_fit",
               "pi_ratio", "reels_cal", "nexafs_edge")

    def adopt_aliases(self, fid, aliases):
        """``aliases`` maps a sample the reader made out of a plainer one
        (``SpectrumFile.sample_aliases``) to that plain name. Whatever was
        saved for the plain sample (notes, BE shifts, markers, sputter
        settings ...) is copied to the new one where it has nothing of its own,
        so edits made before a file's runs were told apart still apply. The
        plain sample's rename is not copied (the new name says which run it
        is). Returns how many entries were copied."""
        n = 0
        for new, plain in (aliases or {}).items():
            if new == plain:
                continue
            ps, ns = sample_key(fid, plain), sample_key(fid, new)
            for table in self.ADOPTED:
                d = getattr(self, table)
                for key in list(d):
                    if key == ps:
                        target = ns
                    elif key.startswith(ps + "|"):
                        target = ns + key[len(ps):]
                    else:
                        continue
                    if target not in d:
                        d[target] = copy.deepcopy(d[key])
                        n += 1
        return n

    # -- names -----------------------------------------------------------------
    def sample_label(self, fid, sample):
        return self.sample_names.get(sample_key(fid, sample)) or sample

    def region_label(self, fid, sample, name):
        return self.region_names.get(region_key(fid, sample, name)) or name

    def set_name(self, table, key, new, original):
        """Rename (or, with an empty / unchanged name, reset)."""
        new = (new or "").strip()
        d = getattr(self, table)
        if not new or new == original:
            d.pop(key, None)
        else:
            d[key] = new

    # -- notes -----------------------------------------------------------------
    def set_note(self, table, key, text):
        text = (text or "").strip()
        d = getattr(self, table)
        if text:
            d[key] = text
        else:
            d.pop(key, None)

    def notes_for(self, fid, sample, name):
        """Sample note and region note, joined for display."""
        parts = [self.sample_notes.get(sample_key(fid, sample), ""),
                 self.region_notes.get(region_key(fid, sample, name), "")]
        return " / ".join(p for p in parts if p)

    # -- metadata edits ---------------------------------------------------------
    def edited_keys(self, fid, pos):
        return set(self.md_edits.get(meta_key(fid, pos), {}))

    def set_meta(self, fid, pos, key, value):
        d = self.md_edits.setdefault(meta_key(fid, pos), {})
        d[key] = str(value)

    def reset_meta(self, fid, pos, key=None):
        k = meta_key(fid, pos)
        if key is None:
            self.md_edits.pop(k, None)
        elif k in self.md_edits:
            self.md_edits[k].pop(key, None)
            if not self.md_edits[k]:
                del self.md_edits[k]

    def apply_metadata(self, fid, pos, region, md):
        """``md`` (a region's metadata dict) with names, notes, the BE shift
        and any edits applied. Returns a new dict; key order is kept."""
        out = dict(md)
        if "Sample" in out:
            out["Sample"] = self.sample_label(fid, region.sample)
        if "Region" in out:
            out["Region"] = self.region_label(fid, region.sample, region.name)
        shift = self.shift_for(fid, region.sample, region.name,
                               getattr(region, "calibration_shift", 0.0))
        if shift:
            out["BE shift (eV)"] = f"{shift:+.3f}"
        if region.etch_level is not None:
            sset = self.sputter.get(sample_key(fid, region.sample))
            if sset:
                for k, v in sputter_mod.metadata_rows(
                        sset, region.etch_time).items():
                    out[k] = v
        if getattr(region, "is_nexafs", False):
            edge = self.edge_for(fid, region.sample, region.name)
            res = nexafs.process(region, self.nexafs_ring, edge)
            if res is not None and res.mean is not None:
                out["Normalisation"] = (
                    f"scaled to the mean ring current ({res.mean:.1f}; the "
                    "unit is not recorded)")
            if res is not None and res.edge is not None:
                out["Edge normalisation"] = nexafs.edge_text(edge)
                if "step" in res.edge:
                    out["Edge step"] = (f"{res.edge['step']:.6g} "
                                        f"{region.count_units}".strip())
        iss = self.iss_for(fid) if getattr(region, "is_iss", False) else {}
        if iss.get("ion"):
            out["ISS ion"] = iss["ion"]
        if iss.get("e0"):
            out["ISS beam energy (eV)"] = f"{iss['e0']:g}"
            out.pop("ISS calibrated beam energy (eV)", None)  # theirs is *the* energy
        if iss.get("theta"):
            out["ISS scattering angle (°)"] = f"{iss['theta']:g}"
        rl = self.reels.get(region_key(fid, region.sample, region.name))
        if rl and rl.get("gap") is not None:
            out["REELS band gap (eV)"] = f"{rl['gap']:.2f}"
            out["REELS elastic peak (eV KE)"] = f"{rl['elastic']:.2f}"
        pl = self.plasmon.get(region_key(fid, region.sample, region.name))
        if pl:
            import reels
            out["REELS plasmon energy (eV)"] = f"{pl['ep']:.2f}"
            out["REELS plasmon method"] = reels.method_text(pl)
            out.setdefault("REELS elastic peak (eV KE)",
                           f"{pl['elastic']:.2f}")
        rc = self.reels_cal.get(region_key(fid, region.sample, region.name))
        if rc:
            import reels
            out["REELS energy calibration"] = reels.calibration_text(rc)
        lf = self.loss_fit.get(region_key(fid, region.sample, region.name))
        if lf:
            import reels
            out.update(reels.loss_fit_metadata(lf))
        pr = self.pi_ratio.get(region_key(fid, region.sample, region.name))
        if pr:
            import reels
            out.update(reels.ratio_metadata(pr, self.reels_refs))
        note = self.notes_for(fid, region.sample, region.name)
        if note:
            out["Notes"] = note
        for k, v in self.md_edits.get(meta_key(fid, pos), {}).items():
            out[k] = v
        return out

    # -- binding-energy shifts ---------------------------------------------------
    def shift_for(self, fid, sample, name, default=0.0):
        """The BE shift for a region: the most specific scope wins (region,
        then sample, then the whole file). ``default`` is returned when the
        user set none: pass the region's ``calibration_shift`` so the charge
        correction recorded in the file applies until the user's own
        calibration replaces it (``None`` tells the two cases apart)."""
        for key in (region_key(fid, sample, name), sample_key(fid, sample),
                    str(fid)):
            if key in self.shifts:
                return float(self.shifts[key])
        return default

    def set_shift(self, scope_key, shift, entry=None):
        if abs(shift) < 1e-12:
            self.shifts.pop(scope_key, None)
        else:
            self.shifts[scope_key] = float(shift)
        if entry:
            self.calibration.append(dict(entry, scope=scope_key,
                                         shift=float(shift)))

    # -- markers -----------------------------------------------------------------
    def markers_for(self, fid, sample, name):
        return list(self.markers.get(region_key(fid, sample, name), []))

    def add_marker(self, fid, sample, name, be, label, kin=False, src=""):
        """A peak marker. ``be`` is a binding energy, or with ``kin=True`` a
        kinetic energy (ISS peaks): it then does not move with an energy
        calibration. ``src="casa"`` marks a name taken from the file's own
        CasaXPS regions, which needs no alternative lines beside it."""
        lst = self.markers.setdefault(region_key(fid, sample, name), [])
        if not any(abs(m["be"] - be) < 1e-6 and m["label"] == label
                   for m in lst):
            m = {"be": float(be), "label": label}
            if kin:
                m["kin"] = True
            if src == "casa":
                m["src"] = "casa"
            lst.append(m)

    # -- REELS band gaps (per spectrum) ---------------------------------------
    def reels_for(self, fid, sample, name):
        """The stored REELS construction of a spectrum, or None."""
        d = self.reels.get(region_key(fid, sample, name))
        return dict(d) if d else None

    def set_reels(self, fid, sample, name, data):
        key = region_key(fid, sample, name)
        clean = _clean_reels(data)
        if clean:
            self.reels[key] = clean
        else:
            self.reels.pop(key, None)

    # -- REELS plasmon energies (per spectrum) --------------------------------
    def plasmon_for(self, fid, sample, name):
        """The stored plasmon result of a spectrum, or None."""
        d = self.plasmon.get(region_key(fid, sample, name))
        return copy.deepcopy(d) if d else None

    def set_plasmon(self, fid, sample, name, data):
        key = region_key(fid, sample, name)
        clean = _clean_plasmon(data)
        if clean:
            self.plasmon[key] = clean
        else:
            self.plasmon.pop(key, None)

    # -- REELS loss-region fits and layer estimates (per spectrum) ------------
    def loss_fit_for(self, fid, sample, name):
        """The stored loss-region fit of a spectrum, or None."""
        d = self.loss_fit.get(region_key(fid, sample, name))
        return copy.deepcopy(d) if d else None

    def set_loss_fit(self, fid, sample, name, data):
        key = region_key(fid, sample, name)
        clean = _clean_loss_fit(data)
        if clean:
            self.loss_fit[key] = clean
        else:
            self.loss_fit.pop(key, None)

    # -- REELS pi / (sigma+pi) ratios (per spectrum) and their references -----
    def pi_ratio_for(self, fid, sample, name):
        """The stored π / (σ+π) ratio result of a spectrum, or None."""
        d = self.pi_ratio.get(region_key(fid, sample, name))
        return copy.deepcopy(d) if d else None

    def set_pi_ratio(self, fid, sample, name, data):
        key = region_key(fid, sample, name)
        clean = _clean_pi_ratio(data)
        if clean:
            self.pi_ratio[key] = clean
        else:
            self.pi_ratio.pop(key, None)

    def set_reels_ref(self, which, data):
        """Store (or, with None, remove) the workbook's HOPG or diamond
        reference ratio. ``which`` is "hopg" or "diamond"."""
        clean = _clean_reels_ref(data)
        if which in ("hopg", "diamond") and clean:
            self.reels_refs[which] = clean
        else:
            self.reels_refs.pop(which, None)

    # -- REELS energy calibration (per spectrum) -------------------------------
    # Results of the REELS dialog keep the elastic peak they were measured from,
    # in the frame the spectrum was drawn in. Moving the energy scale moves that
    # frame, so a stored elastic peak moves with it: the loss coordinates of the
    # results (a band-gap pick, a window) are unchanged.
    REELS_RESULTS = ("reels", "plasmon", "loss_fit", "pi_ratio")

    def reels_cal_for(self, fid, sample, name):
        """The stored REELS energy calibration of a spectrum, or None."""
        d = self.reels_cal.get(region_key(fid, sample, name))
        return dict(d) if d else None

    def reels_shift_for(self, fid, sample, name):
        """The kinetic-energy shift (eV) applied to a REELS spectrum (0 when it
        is not calibrated)."""
        d = self.reels_cal.get(region_key(fid, sample, name))
        return d["shift"] if d else 0.0

    def set_reels_cal(self, fid, sample, name, data):
        """Store (or, with None, remove) a spectrum's calibration and move the
        elastic peak of every REELS result stored for it by the change of
        shift. Returns the change (eV)."""
        key = region_key(fid, sample, name)
        old = self.reels_cal.get(key)
        clean = _clean_reels_cal(data)
        delta = (clean["shift"] if clean else 0.0) - (
            old["shift"] if old else 0.0)
        if clean:
            self.reels_cal[key] = clean
        else:
            self.reels_cal.pop(key, None)
        if delta:
            for table in self.REELS_RESULTS:
                d = getattr(self, table).get(key)
                if d and d.get("elastic") is not None:
                    d["elastic"] = d["elastic"] + delta
        return delta

    def remove_marker(self, fid, sample, name, be, label):
        key = region_key(fid, sample, name)
        lst = [m for m in self.markers.get(key, [])
               if not (m["label"] == label and abs(m["be"] - be) < 1e-6)]
        if lst:
            self.markers[key] = lst
        else:
            self.markers.pop(key, None)

    def clear_markers(self, fid, sample, name):
        self.markers.pop(region_key(fid, sample, name), None)

    # -- sputter settings (per sample) ---------------------------------------------
    def sputter_for(self, fid, sample):
        """The sample's sputter settings (complete dict), or None."""
        s = self.sputter.get(sample_key(fid, sample))
        return sputter_mod.sanitise(s) if s else None

    def set_sputter(self, fid, sample, settings):
        key = sample_key(fid, sample)
        if settings is None or sputter_mod.is_empty(settings):
            self.sputter.pop(key, None)
        else:
            self.sputter[key] = sputter_mod.sanitise(settings)

    # -- NEXAFS edge normalisation (per spectrum) -----------------------------------
    def edge_for(self, fid, sample, name):
        """The saved edge normalisation of a NEXAFS spectrum (``nexafs``), or
        {} when the user set none."""
        return nexafs.sanitise_edge(
            self.nexafs_edge.get(region_key(fid, sample, name)))

    def set_edge(self, fid, sample, name, params):
        """Save (or, with {} / invalid parameters, remove) the edge
        normalisation of one spectrum."""
        key, clean = region_key(fid, sample, name), nexafs.sanitise_edge(params)
        if clean:
            self.nexafs_edge[key] = clean
        else:
            self.nexafs_edge.pop(key, None)

    # -- repeat groups (samples that repeat one another; see repeats.py) ------------
    def set_repeat_group(self, name, kind, sample_keys):
        """Create or replace a group; an empty list of samples removes it."""
        name = (name or "").strip()
        if not name:
            return
        g = repeats.sanitise_group({"kind": kind, "samples": sample_keys})
        if g:
            self.repeat_groups[name] = g
        else:
            self.repeat_groups.pop(name, None)

    def remove_repeat_group(self, name):
        self.repeat_groups.pop(name, None)

    # -- instrument settings for NeXus (per file) ------------------------------------
    def instrument_for(self, fid):
        """The file's instrument settings (only what is set), or {}."""
        return nexus_settings.sanitise(self.instrument.get(str(fid)))

    def set_instrument(self, fid, settings):
        s = nexus_settings.sanitise(settings)
        if s:
            self.instrument[str(fid)] = s
        else:
            self.instrument.pop(str(fid), None)

    # -- ion scattering settings the user confirmed (per file) ---------------------
    def iss_for(self, fid):
        """The ion, beam energy (eV) and scattering angle (degrees) saved for
        a file (only what is set), or {}."""
        return elements.sanitise_iss(self.iss.get(str(fid)))

    def set_iss(self, fid, settings):
        s = elements.sanitise_iss(settings)
        if s:
            self.iss[str(fid)] = s
        else:
            self.iss.pop(str(fid), None)

    # -- how a plain column-text file is read (per file) ---------------------------
    def import_for(self, fid):
        """The options the file was imported with (``columntext``), or {}."""
        opts = self.imports.get(str(fid))
        return columntext.sanitise_options(opts) if opts else {}

    def set_import(self, fid, options):
        if options:
            self.imports[str(fid)] = columntext.sanitise_options(options)
        else:
            self.imports.pop(str(fid), None)

    # -- state ---------------------------------------------------------------------
    def is_empty(self):
        return not (self.sample_names or self.region_names
                    or self.sample_notes or self.region_notes or self.md_edits
                    or self.shifts or self.calibration or self.markers
                    or self.calibration_statement or self.experiment_notes
                    or self.sputter or self.reels or self.plasmon
                    or self.loss_fit or self.pi_ratio or self.reels_refs
                    or self.reels_cal
                    or self.instrument
                    or self.iss or self.imports or self.nexafs_ring
                    or self.nexafs_edge or self.repeat_groups)

    def copy(self):
        return copy.deepcopy(self)

    def to_json(self):
        d = {"version": VERSION}
        for name in ("sample_names", "region_names", "sample_notes",
                     "region_notes", "md_edits", "shifts", "calibration",
                     "calibration_statement", "markers", "experiment_notes",
                     "sputter", "reels", "plasmon", "loss_fit", "pi_ratio",
                     "reels_refs", "reels_cal", "instrument", "iss",
                     "imports"):
            d[name] = copy.deepcopy(getattr(self, name))
        if self.nexafs_ring:
            d["nexafs_ring"] = True
        if self.nexafs_edge:
            d["nexafs_edge"] = copy.deepcopy(self.nexafs_edge)
        if self.repeat_groups:
            d["repeat_groups"] = copy.deepcopy(self.repeat_groups)
        d.update(self.extra)
        return d

    @classmethod
    def from_json(cls, data):
        """Tolerant loader: wrong types are ignored, unknown keys kept."""
        a = cls()
        if not isinstance(data, dict):
            return a
        known = {"version", "sample_names", "region_names", "sample_notes",
                 "region_notes", "md_edits", "shifts", "calibration",
                 "calibration_statement", "markers", "experiment_notes",
                 "sputter", "reels", "plasmon", "loss_fit", "pi_ratio",
                 "reels_refs", "reels_cal", "instrument", "iss", "imports",
                 "nexafs_ring", "nexafs_edge", "repeat_groups"}
        for name in ("sample_names", "region_names", "sample_notes",
                     "region_notes"):
            v = data.get(name)
            if isinstance(v, dict):
                setattr(a, name, {str(k): str(x) for k, x in v.items()})
        v = data.get("md_edits")
        if isinstance(v, dict):
            a.md_edits = {str(k): {str(i): str(j) for i, j in e.items()}
                          for k, e in v.items() if isinstance(e, dict)}
        v = data.get("shifts")
        if isinstance(v, dict):
            for k, x in v.items():
                try:
                    a.shifts[str(k)] = float(x)
                except (TypeError, ValueError):
                    pass
        v = data.get("calibration")
        if isinstance(v, list):
            a.calibration = [e for e in v if isinstance(e, dict)]
        v = data.get("markers")
        if isinstance(v, dict):
            for k, lst in v.items():
                if not isinstance(lst, list):
                    continue
                good = []
                for m in lst:
                    try:
                        item = {"be": float(m["be"]), "label": str(m["label"])}
                        if m.get("kin") is True:
                            item["kin"] = True
                        if m.get("src") == "casa":
                            item["src"] = "casa"
                        good.append(item)
                    except (KeyError, TypeError, ValueError):
                        pass
                if good:
                    a.markers[str(k)] = good
        v = data.get("reels")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = _clean_reels(x)
                if clean:
                    a.reels[str(k)] = clean
        v = data.get("plasmon")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = _clean_plasmon(x)
                if clean:
                    a.plasmon[str(k)] = clean
        v = data.get("loss_fit")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = _clean_loss_fit(x)
                if clean:
                    a.loss_fit[str(k)] = clean
        v = data.get("reels_cal")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = _clean_reels_cal(x)
                if clean:
                    a.reels_cal[str(k)] = clean
        v = data.get("pi_ratio")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = _clean_pi_ratio(x)
                if clean:
                    a.pi_ratio[str(k)] = clean
        v = data.get("reels_refs")
        if isinstance(v, dict):
            for k in ("hopg", "diamond"):
                clean = _clean_reels_ref(v.get(k))
                if clean:
                    a.reels_refs[k] = clean
        v = data.get("sputter")
        if isinstance(v, dict):
            a.sputter = {str(k): sputter_mod.sanitise(x) for k, x in v.items()
                         if isinstance(x, dict)
                         and not sputter_mod.is_empty(x)}
        if data.get("nexafs_ring") is True:
            a.nexafs_ring = True
        v = data.get("nexafs_edge")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = nexafs.sanitise_edge(x)
                if clean:
                    a.nexafs_edge[str(k)] = clean
        a.repeat_groups = repeats.sanitise_groups(data.get("repeat_groups"))
        v = data.get("instrument")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = nexus_settings.sanitise(x)
                if clean:
                    a.instrument[str(k)] = clean
        v = data.get("iss")
        if isinstance(v, dict):
            for k, x in v.items():
                clean = elements.sanitise_iss(x)
                if clean:
                    a.iss[str(k)] = clean
        v = data.get("imports")
        if isinstance(v, dict):
            for k, x in v.items():
                if isinstance(x, dict) and x:
                    a.imports[str(k)] = columntext.sanitise_options(x)
        for name in ("calibration_statement", "experiment_notes"):
            if isinstance(data.get(name), str):
                setattr(a, name, data[name])
        a.extra = {k: v for k, v in data.items() if k not in known}
        return a
