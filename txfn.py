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
"""
from __future__ import annotations

import math

REF_KE = 1118.0         # eV; CasaXPS divides T by its value at this energy


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
