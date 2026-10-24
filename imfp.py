"""Inelastic mean free path (IMFP) proxies for the Scofield RSF fallback
(``quant.py``'s ``"scofield_tpp2m"``/``"scofield_ke06"`` tiers), Tk-free.

The standard XPS quantification relation is n_A proportional to I_A / (sigma_A
* lambda_A(KE) * T_A): a pure Scofield cross section (``rsf.py``'s
"scofield" library, C 1s = 1.0, a theoretical cross section only) has no
lambda term at all, so comparing two lines whose kinetic energies differ a
lot (round 6's own motivating case: O 1s ~= 950 eV KE vs. Pt 4f ~= 1420 eV
KE) with cross sections alone ignores a real, non-negligible factor (the
two functions here differ by ~35% for that KE pair). Both functions below
give a *relative* factor to multiply onto a Scofield cross section, never a
number to be used or shown on its own; onto that library only -- never onto
the Kratos Axis F1s library or a file's own recorded RSF, both empirical
(Wagner-style) values already calibrated against a real instrument, with
their own implicit lambda(E)*T(E) dependence baked in by being measured --
multiplying either of those by an independent term here would double-count
it, not correct it.

**Provenance, checked directly against KherveFitting's own source**
(``D:\\Software\\KherveFitting\\_internal\\libraries\\Peak_Functions.py``,
class ``AtomicConcentrations``, plus a standalone ``TPP-2M.py`` script) --
neither of the two functions here is a proven, actively-exercised feature
the way ``rsf.py``'s RSF libraries were:

- ``imfp_nm``: KherveFitting's ``calculate_imfp_tpp2m`` is a correctly
  formed TPP-2M equation (S. Tanuma, C.J. Powell, D.R. Penn, *Surf.
  Interface Anal.* 21, 165-176, 1993 -- the formula shape here matches the
  published equation, independently checked, and the values it gives
  (e.g. ~2.37 nm at 1000 eV) sit in the expected 0.5-4 nm range for
  50-2000 eV) using fixed "average matrix" constants (cited to Briggs &
  Grant, *Surface Analysis by XPS and AES*, 2nd ed., Wiley 2003, p.84-85 --
  not independently checked past that citation) rather than a real
  material's own density/molecular weight/valence-electron count/band gap,
  which sidesteps the fact that eXPoSe has no way to know an arbitrary
  sample's own values for those (KherveFitting's *other* TPP-2M variant,
  ``calculate_imfp_tpp2m_WITHOUT_VALUES_BUT_GOOD``, takes real per-material
  values as parameters, but nothing in KherveFitting ever calls it and no
  material database exists there to supply them; it also uses ``log2``
  where the published formula needs natural ``log`` -- a real
  transcription bug in code this module does not port). Exhaustive search
  of the whole KherveFitting install found **no call site for either
  function outside their own definitions, and no menu item or dialog
  anywhere that reaches them** -- both are dead code there; KherveFitting's
  own actual atomic-percent math is plain user-typed RSF/TXFN, no lambda
  term anywhere.
- ``ke_power_factor``: not derived from any formula at all -- it is
  literally what KherveFitting's own ``calculate_imfp_tpp2m`` docstring
  references as *Thermo Avantage's* own convention (a stray, otherwise
  unused ``imfp2 = imfp * 26.2`` line in that same function, commented
  "Avantage add a scaling factor so that corrected area matches the one
  obtained with KE^0.6 -- To compare KherveFitting with Avantage we will
  also apply this factor"). Included as a second, simpler, faster
  approximation for comparison, not claimed to be as physically complete
  as TPP-2M.

Only *ratios* between lines matter for atomic percent (the same reason the
Scofield table's own absolute scale is arbitrary -- C 1s = 1.0 is a
convention, not a physical constant), so neither function needs, or is
given, a proportionality constant: they return exactly what the formula/
convention states, nothing normalised or rescaled.
"""

import math

# Briggs & Grant's "average matrix" constants (see the module docstring)
N_V = 4.684
RHO = 6.767          # g/cm^3
MOLAR_MASS = 137.51  # g/mol
E_GAP = 0.0          # eV

MIN_KE, MAX_KE = 50.0, 2000.0   # TPP-2M's own stated validity range (eV)

KE_POWER_EXPONENT = 0.6   # Thermo Avantage's own convention -- see the
                          # module docstring's "ke_power_factor" paragraph


def _tpp2m(e, n_v, rho, molar_mass, e_gap, ep2_coeff=None):
    """The TPP-2M equation (Tanuma, Powell, Penn 1993) in nm for a kinetic
    energy ``e`` (eV); None when the denominator is not positive (very low
    energies). ``ep2_coeff`` is the plasmon constant in Ep^2 = k U (the
    published 28.8^2 = 829.44 unless given)."""
    u = n_v * rho / molar_mass
    ep2 = (829.44 if ep2_coeff is None else ep2_coeff) * u
    beta = (-0.10 + 0.944 / math.sqrt(ep2 + e_gap ** 2)
            + 0.069 * rho ** 0.1)
    gamma = 0.191 * rho ** -0.5
    c = 1.97 - 0.91 * u
    d = 53.4 - 20.8 * u
    den = ep2 * (beta * math.log(gamma * e) - c / e + d / e ** 2)
    if not den > 0:
        return None
    return e / den / 10.0


def imfp_nm(kinetic_energy):
    """The TPP-2M inelastic mean free path, in nm, for an electron of this
    kinetic energy (eV), in the fixed "average matrix" -- None when
    ``kinetic_energy`` is missing or outside TPP-2M's own stated
    50-2000 eV validity range (never an unreliable extrapolation)."""
    if not kinetic_energy or not (MIN_KE <= kinetic_energy <= MAX_KE):
        return None
    return _tpp2m(float(kinetic_energy), N_V, RHO, MOLAR_MASS, E_GAP)


# Avantage's ALTHERMO1 quantification divides each peak by a TPP-2M energy
# term in place of the Scofield libraries' KE^0.6. Fitted to a Peak Table
# Avantage made (an S 2p fit, Na sulfite and sulfate: Area (N) = Area (P) /
# (SF x term)), the term is ``REF^0.6 x lambda(KE) / lambda(REF)`` -- the
# TPP-2M path relative to its value at a reference energy, scaled so that it
# equals KE^0.6 there -- with ``REF`` = 1118 eV (the reference energy of
# Avantage's transmission function, ``txfn.REF_KE``): the four peaks agree to
# 0.04 % (the table's own rounding), the best reference energy lies in
# 1115-1119 eV, and the average-matrix and CasaXPS-matrix paths give the same
# result (the ratio lambda(KE) / lambda(REF) hardly depends on the matrix).
AVANTAGE_REF_KE = 1118.0


def avantage_energy_term(kinetic_energy):
    """The energy term of Avantage's ALTHERMO1 quantification (see above),
    or None outside TPP-2M's 50-2000 eV range."""
    lam, ref = imfp_nm(kinetic_energy), imfp_nm(AVANTAGE_REF_KE)
    if not lam or not ref:
        return None
    return AVANTAGE_REF_KE ** 0.6 * lam / ref


# CasaXPS's own "TPP2M" escape-depth type, read from CasaXPS.exe (Ghidra,
# function 0x1403cc280, case 5): one fixed material, diamond-like carbon, and
# no energy limits. Ep^2 = 829.4 U, U = Nv rho / M.
CASA_N_V, CASA_RHO, CASA_MOLAR_MASS, CASA_E_GAP = 4.0, 3.513, 12.01, 5.4
CASA_EP2_COEFF = 829.4


def casa_imfp_nm(kinetic_energy):
    """CasaXPS's TPP-2M mean free path in nm (its fixed diamond-like matrix,
    no validity range, as the program has none); None when the energy is not
    positive or the equation gives no positive value."""
    if not kinetic_energy or kinetic_energy <= 0:
        return None
    return _tpp2m(float(kinetic_energy), CASA_N_V, CASA_RHO,
                  CASA_MOLAR_MASS, CASA_E_GAP, CASA_EP2_COEFF)


# CasaXPS "ESCAPE DEPTH TYPE" codes (the integer stored in the VAMAS block's
# additional parameters) and the labels CasaXPS shows ("Escape D.C.: ...").
ESCAPE_NONE, ESCAPE_KE_EXP, ESCAPE_EAL, ESCAPE_IMFP = 0, 1, 2, 3
ESCAPE_CH_FILM, ESCAPE_TPP2M, ESCAPE_LAYER = 4, 5, 6
ESCAPE_NAMES = {0: "None", 1: "KE exponent", 2: "EAL", 3: "IMFP",
                4: "CH film", 5: "TPP-2M", 6: "Layer"}


def _eal(ke):                       # CasaXPS's effective attenuation length
    return 0.65 + 0.007 * ke ** 0.93


def _casa_imfp_line(ke):            # the program's "IMFP" type
    return 0.73 + 0.0095 * ke ** 0.872


def casa_escape_factor(kind, ke, ref=1.0, p=0.0):
    """The factor CasaXPS multiplies an intensity by for escape-depth
    correction, lambda(ref) / lambda(KE), read from CasaXPS.exe. ``kind`` is
    the ESCAPE DEPTH TYPE code, ``ke`` the kinetic energy of the point,
    ``ref`` the reference kinetic energy (1 eV unless a reference binding
    energy is set) and ``p`` the block's MFP Exponent field: the exponent for
    type 1 (stored negative, -0.6 for KE^0.6), the overlayer thickness for CH
    film and Layer (in the EAL's unit, nm). Quantification divides the area
    by 1 / factor, i.e. its "MFP" is lambda(KE) / lambda(ref). 1.0 for None,
    an unknown type or an exponent of zero, as in the program."""
    if kind == ESCAPE_KE_EXP:
        return (ke / ref) ** p if abs(p) > 1e-200 else 1.0
    if kind == ESCAPE_EAL:
        l, lr = _eal(ke), _eal(ref)
        return lr / l if l > 1e-200 else l
    if kind == ESCAPE_IMFP:
        l, lr = _casa_imfp_line(ke), _casa_imfp_line(ref)
        return lr / l if l > 1e-200 else l
    if kind == ESCAPE_TPP2M:
        l = _tpp2m(ke, CASA_N_V, CASA_RHO, CASA_MOLAR_MASS, CASA_E_GAP,
                   CASA_EP2_COEFF)
        lr = _tpp2m(ref, CASA_N_V, CASA_RHO, CASA_MOLAR_MASS, CASA_E_GAP,
                    CASA_EP2_COEFF)
        if l is None or lr is None:
            return 1.0
        return lr / l if l > 1e-200 else l
    if kind in (ESCAPE_CH_FILM, ESCAPE_LAYER):
        scale = 4.0 ** 0.38 if kind == ESCAPE_CH_FILM else 1.0
        l, lr = _eal(ke), _eal(ref)
        if not l > 1e-200:
            return l
        if kind == ESCAPE_CH_FILM:
            return (lr / l) * math.exp(p * scale / l) / math.exp(
                p * scale / lr)
        if not p > 1e-200:
            return l * (1 - math.exp(-p / l))
        return lr * (1 - math.exp(-p / lr)) / (l * (1 - math.exp(-p / l)))
    return 1.0


def ke_power_factor(kinetic_energy):
    """A relative mean-free-path proxy, ``kinetic_energy ** 0.6`` -- not a
    physical length (no nm/Angstrom unit, no material assumed at all), see
    the module docstring's "ke_power_factor" paragraph. None when
    ``kinetic_energy`` is not a positive number."""
    if not kinetic_energy or kinetic_energy <= 0:
        return None
    return float(kinetic_energy) ** KE_POWER_EXPONENT
