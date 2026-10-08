"""Plain-language wording for Kratos instrument constants (Tk-free).

Kratos files name settings with internal constants (``F_HSA_LENS_HYBRID``);
the Harwell / CasaXPS VAMAS export copies them without the ``F_``
(``HSA_LENS_HYBRID``). Every reader passes the lens settings through
``friendly`` so the methods text, the metadata and the exports read the same
whichever file route the data came by.

To change a wording, edit ``LENS_MODES`` / ``MAGNIFICATIONS``. Keys are the
constant without ``F_`` and in upper case. A constant that is not listed is
tidied (prefix dropped, underscores to spaces, capitalised), never guessed
or left blank; text that is not a Kratos constant is returned as written.
"""

from __future__ import annotations

import re

# ``HSA Lens Mode`` (field 3049): every value seen in 7 342 Kratos .kal files
LENS_MODES = {
    "HSA_LENS_HYBRID": "Hybrid (electrostatic and magnetic immersion lens)",
    "HSA_LENS_MAGNETIC": "Magnetic (magnetic immersion lens)",
    "HSA_LENS_ELECTROSTATIC": "Electrostatic (magnetic lens off)",
    "HSA_LENS_ISS": "ISS (magnetic lens off)",
}

# ``MHSA Lens Mode`` (field 3050, stigmatic imaging maps): the three fields
# of view of the imaging lens, FoV1 the widest
MAGNIFICATIONS = {
    "MHSA_LOW_MAGN": "Low magnification (FoV1)",
    "MHSA_MEDIUM_MAGN": "Medium magnification (FoV2)",
    "MHSA_HIGH_MAGN": "High magnification (FoV3)",
}

_CONSTANT = re.compile(r"^(?:F_)?(M?HSA_[A-Z0-9_]+)$", re.I)
_PREFIXES = ("MHSA_", "HSA_LENS_", "HSA_")


def is_iss_lens(text) -> bool:
    """True when a lens-mode text is the Kratos ISS lens: the constant (with
    or without ``F_``), its wording in ``LENS_MODES`` (whatever it has been
    changed to) or a bare "ISS". Reports and the readers decide "this is an
    ion scattering spectrum" with this, never by comparing wording."""
    s = (text or "").strip()
    if not s:
        return False
    return (s.upper() == "ISS" or s == LENS_MODES["HSA_LENS_ISS"]
            or friendly(s) == LENS_MODES["HSA_LENS_ISS"])


def friendly(text) -> str:
    """The wording for a Kratos lens constant; anything else comes back
    unchanged (stripped), "" for empty."""
    s = (text or "").strip()
    m = _CONSTANT.match(s)
    if not m:
        return s
    key = m.group(1).upper()
    table = MAGNIFICATIONS if key.startswith("MHSA_") else LENS_MODES
    if key in table:
        return table[key]
    for p in _PREFIXES:
        if key.startswith(p):
            key = key[len(p):]
            break
    return key.replace("_MAGN", "_MAGNIFICATION").replace("_", " ").capitalize()
