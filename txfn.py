"""Thermo Avantage transmission function (Tk-free, pure Python).

Avantage stores four polynomial coefficients per DataSpace and, for the
K-Alpha / NEXSA / HYPULSE instruments, has already divided the intensities by
the transmission function (``DS_ANPROPID_TXFN_REMOVED``). CasaXPS's import of
an Avantage ``.avg`` puts the function back as a ``Transmission`` column and
rescales the counts so that dividing by that column gives the stored values.

* The polynomial is Thermo's *Avantage Data Processing Reference Manual* v5,
  p. 35: ``log10[T / (PE * XSF)] = a + b x + c x^2 + d x^3`` with
  ``x = log10(KE / PE)`` (coefficients in increasing powers).
* The scale ``T = PE * 10**poly / 100`` and the reference energy
  ``REF_KE`` = 1118.0 eV are **not in the manual**: they were read off
  CasaXPS's own import of C 1s, Ga 2p, Ga 3d, GaLM2, O 1s (pass 50) and a
  survey (pass 150) of one K-Alpha+ file. T agrees with CasaXPS's column to
  1e-6 and the counts, ``stored * T(KE) / T(REF_KE)``, to 5e-6 (rounding of
  the printed values) on every point.

**Which function does a CasaXPS VAMAS column hold?** (``classify``; checked on
the three CasaXPS exports of one Nexsa G2 NaCl run, Standard lens, pass
energies 50-200 eV, in ``D:\\Temp\\metadataapp\\datatest``; the counts are
identical in all of them, only the second variable differs.)

* ``avantage``: Avantage's polynomial above, as stored in the dump
  (``nacl-avg-transcurve``; equals ``transmission()`` to 5e-6).
* ``casa_formula``: CasaXPS's "New VG Transmission Function: C0 .. C3"
  (``nacl-quant-recalculated-TF``): the same form with new coefficients and
  **no** ``/ 100``, ``T = PE * 10**poly(log10(KE / PE))`` (1e-6).
* ``casa_file``: a "Transmission file: <measured .vms>" line
  (``nacl-quant-rel-TF-ESCALAB250Xi``, which also prints C0 .. C3): the column
  follows **neither** formula and is not a function of KE / PE alone, so it is
  read from the measured file, which is not part of the export. It cannot be
  rebuilt here; the column is used as written.
"""
from __future__ import annotations

import math
import re

REF_KE = 1118.0         # eV; CasaXPS divides T by its value at this energy
MATCH_TOL = 1e-4        # relative; the printed values have 6 digits

_NEW_VG_RE = re.compile(r"New\s+VG\s+Transmission\s+Function\s*:\s*(.+)", re.I)
_COEFF_RE = re.compile(r"\bC(\d+)\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)")
_FILE_RE = re.compile(r"^\s*Transmission\s+file\s*:\s*(.+?)\s*$", re.I)


def transmission(ke, pass_energy, coeffs):
    """T(KE) for one kinetic energy (eV), the pass energy (eV) and the
    polynomial coefficients; None when it cannot be evaluated."""
    if not (ke and ke > 0 and pass_energy and pass_energy > 0 and coeffs):
        return None
    x = math.log10(ke / pass_energy)
    y = sum(c * x ** i for i, c in enumerate(coeffs))
    try:
        return pass_energy * 10.0 ** y / 100.0
    except OverflowError:
        return None


def casa_form(region):
    """``(counts, transmission)`` as CasaXPS's Avantage import writes them, or
    None when the region cannot be put in that form: only data whose
    transmission function the source says was removed, with the coefficients,
    a pass energy and a kinetic-energy axis."""
    tx = region.extra.get("txfn") or {}
    coeffs = tx.get("coeffs")
    ke = region.kinetic_energy
    if not tx.get("removed") or not coeffs or not ke or not region.counts:
        return None
    pe = region.pass_energy
    ref = transmission(REF_KE, pe, coeffs)
    t = [transmission(e, pe, coeffs) for e in ke]
    if not ref or any(v is None for v in t):
        return None
    return [c * v / ref for c, v in zip(region.counts, t)], t


def casa_formula(ke, pass_energy, coeffs):
    """CasaXPS's "New VG Transmission Function": ``PE * 10**poly(x)`` with
    ``x = log10(KE / PE)``; None when it cannot be evaluated."""
    if not (ke and ke > 0 and pass_energy and pass_energy > 0 and coeffs):
        return None
    x = math.log10(ke / pass_energy)
    try:
        return pass_energy * 10.0 ** sum(c * x ** i
                                         for i, c in enumerate(coeffs))
    except OverflowError:
        return None


def parse_casa_lines(lines):
    """What CasaXPS writes about a transmission function in a block comment:
    ``{"coeffs": [C0 ..] or None, "file": basename}``. Only the file's name is
    kept (the path is the other machine's)."""
    coeffs, name = None, ""
    for line in lines or ():
        m = _NEW_VG_RE.search(line)
        if m:
            found = {int(i): float(v) for i, v in _COEFF_RE.findall(m.group(1))}
            if found and sorted(found) == list(range(len(found))):
                coeffs = [found[i] for i in range(len(found))]
            continue
        m = _FILE_RE.match(line)
        if m:
            name = re.split(r"[\\/]", m.group(1).strip())[-1].strip()
    return {"coeffs": coeffs, "file": name}


def _matches(column, model):
    if any(v is None for v in model):
        return False
    return all(abs(c - v) <= MATCH_TOL * abs(v) for c, v in zip(column, model))


def classify(ke, pass_energy, column, avantage_coeffs=None, casa=None):
    """Which function a ``Transmission`` column holds, by its numbers (a
    comment that no longer describes the column cannot mislabel it).

    ``{"kind": "avantage" | "casa_formula" | "casa_file", "coeffs": [...],
    "file": name, "vs_avantage": [lo, hi]}``, or None when the column matches
    no stated function and no file is named (nothing is then claimed).
    ``vs_avantage`` is the column over Avantage's function, relative to its
    median across the region, when the coefficients are known."""
    casa = casa or {}
    if not column or not ke or len(ke) != len(column) or not (
            pass_energy and pass_energy > 0):
        return None
    out = {"coeffs": [], "file": casa.get("file", "")}
    av = ([transmission(e, pass_energy, avantage_coeffs) for e in ke]
          if avantage_coeffs else None)
    if av and _matches(column, av):
        out.update(kind="avantage", coeffs=list(avantage_coeffs))
        return out
    cf = casa.get("coeffs")
    if cf and _matches(column, [casa_formula(e, pass_energy, cf) for e in ke]):
        out.update(kind="casa_formula", coeffs=list(cf))
    elif out["file"]:
        out["kind"] = "casa_file"
    else:
        return None
    if av and all(v for v in av):
        ratio = sorted(c / v for c, v in zip(column, av))
        mid = ratio[len(ratio) // 2]
        if mid:
            out["vs_avantage"] = [ratio[0] / mid - 1.0, ratio[-1] / mid - 1.0]
    return out


def describe(src):
    """The metadata statement for a ``classify`` result."""
    kind = src.get("kind")
    span = src.get("vs_avantage")
    worst = max(abs(v) for v in span) if span else 0.0
    # relative to the middle of the region; a spread under 0.5 % is not worth
    # a sentence (a narrow high-resolution scan always looks alike)
    differs = (" (across this region its shape departs from Avantage's by up "
               "to %.0f %%)" % (100 * worst)) if worst >= 0.005 else ""
    if kind == "avantage":
        return "Avantage polynomial (the instrument's function, as stored)"
    if kind == "casa_formula":
        return ("CasaXPS recalculated VG function, PE x 10^poly(log10(KE/PE)) "
                "(C0..C3 " + ", ".join(f"{c:g}" for c in src["coeffs"]) + ")"
                + differs)
    if kind == "casa_file":
        return ("Measured transmission file %s%s%s" % (
            src["file"], ", replacing Avantage's function" if span else "",
            differs))
    return ""
