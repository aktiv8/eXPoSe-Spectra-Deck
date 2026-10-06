"""Plain numbers in columns: a spectrum saved as CSV / ASC / TXT / DAT.

Tk-free and free of the readers, so the Annotations (which store how a file
was read) and the tests can use it alone. ``readers.column_text`` turns what
is found here into regions.

A file of this kind says little about itself, so nothing is assumed beyond
what the numbers show:

* the table is the longest run of lines that are all numbers with the same
  number of columns (a header line above it, junk below it, are set aside and
  reported, never merged in);
* the energy column is the one that runs steadily one way (an axis), wherever
  it stands: ``intensity, energy`` files exist;
* the header, when there is one, names the axis (``BE_`` / ``KE_`` /
  "Binding Energy"), the region (``BE_Cl2p``) and the unit (``CPS``, ``c/s``);
* what no header says (photon energy, pass energy, sample) is left unset until
  the user supplies it.

``Options`` are plain dicts (``guess_options`` / ``sanitise_options``) so they
are stored in the workbook and the file reads the same way after reopening.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

EXTS = (".csv", ".asc", ".txt", ".dat", ".tsv")
MIN_ROWS = 8                    # fewer numeric lines than this is not a spectrum
DELIMITERS = ("\t", ";", ",", None)     # None = any white space

AXES = ("BE", "KE")
AXIS_LABELS = {"BE": "Binding energy", "KE": "Kinetic energy"}
UNITS = ("counts/s", "counts", "a.u.")
ANODES = {"Al Kα": 1486.6, "Mg Kα": 1253.6}
REGULAR = 0.9                   # share of equal steps that makes an axis
MAX_SUGGEST_SPAN = 40.0         # eV: a wider window names too many lines
SURVEY_ALIASES = {"gen", "general", "overview", "widescan", "wide scan"}


# -- reading the table ----------------------------------------------------------
def decode(raw: bytes) -> str:
    """UTF-8 (with or without a byte-order mark), else Windows-1252."""
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        return raw.decode("cp1252", errors="replace")


def _num(tok, decimal_comma=False):
    t = tok.strip()
    if not t or "_" in t:           # float() accepts "1_0"; a header does not
        return None
    if decimal_comma:
        t = t.replace(",", ".")
    try:
        v = float(t)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _split(line, delim):
    return line.split(delim) if delim else line.split()


def _row(line, delim):
    """The numbers of a line, or None when any cell is not one."""
    toks = _split(line, delim)
    if len(toks) < 2:
        return None
    vals = [_num(t, decimal_comma=delim != ",") for t in toks]
    return None if any(v is None for v in vals) else vals


@dataclass
class Table:
    delimiter: object                     # "\t" ";" "," or None (white space)
    columns: list                         # column-major floats
    header: list = field(default_factory=list)      # one label per column
    preamble: list = field(default_factory=list)    # text lines above the table
    first_line: int = 1                   # 1-based line number of the first row
    trailing: int = 0                     # non-blank lines after the table
    trailing_at: int = 0                  # ... starting at this line number

    @property
    def n_cols(self):
        return len(self.columns)

    @property
    def n_rows(self):
        return len(self.columns[0]) if self.columns else 0

    def label(self, c):
        """A column's name for a list: its header, else ``Column n``."""
        h = self.header[c].strip() if c < len(self.header) else ""
        return h or f"Column {c + 1}"


def _longest_run(lines, delim):
    """``(start, end, rows)`` of the longest run of numeric lines with one
    column count (blank lines inside it are skipped); ``None`` when empty."""
    best, cur_start, cur_n, cur_rows, last = None, None, None, [], None
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        vals = _row(line, delim)
        if vals is not None and (cur_n is None or len(vals) == cur_n):
            if cur_n is None:
                cur_start, cur_n, cur_rows = i, len(vals), []
            cur_rows.append(vals)
            last = i
            continue
        if cur_rows and (best is None or len(cur_rows) > len(best[2])):
            best = (cur_start, last, cur_rows)
        if vals is not None:                  # a new run of another width
            cur_start, cur_n, cur_rows, last = i, len(vals), [vals], i
        else:
            cur_start, cur_n, cur_rows = None, None, []
    if cur_rows and (best is None or len(cur_rows) > len(best[2])):
        best = (cur_start, last, cur_rows)
    return best


def read_table(text: str):
    """The numeric table in ``text`` (see the module notes), or ``None`` when
    there are fewer than ``MIN_ROWS`` numeric lines of one width."""
    lines = text.splitlines()
    best = None
    for delim in DELIMITERS:
        run = _longest_run(lines, delim)
        if run and (best is None or len(run[2]) > len(best[1][2])):
            best = (delim, run)
    if best is None or len(best[1][2]) < MIN_ROWS:
        return None
    delim, (start, end, rows) = best
    cols = [list(c) for c in zip(*rows)]
    header, preamble = [], []
    above = [ln for ln in lines[:start] if ln.strip()]
    if above:
        toks = [t.strip() for t in _split(above[-1], delim)]
        if len(toks) == len(cols) and not any(
                _num(t, delim != ",") is not None for t in toks):
            header = toks
            above = above[:-1]
        preamble = [ln.strip() for ln in above][-20:]
    rest = [(i + 1, ln) for i, ln in enumerate(lines[end + 1:], end + 1)
            if ln.strip()]
    return Table(delim, cols, header, preamble, start + 1, len(rest),
                 rest[0][0] if rest else 0)


# -- which column is the energy axis --------------------------------------------
def axis_quality(vals):
    """``(monotonic, regular)``: the share of steps that go the main way, and
    the share of those that equal the median step (2 %)."""
    d = [b - a for a, b in zip(vals, vals[1:])]
    if not d:
        return 0.0, 0.0
    pos = sum(1 for x in d if x > 0)
    neg = sum(1 for x in d if x < 0)
    sign = 1 if pos >= neg else -1
    mono = max(pos, neg) / len(d)
    steps = sorted(abs(x) for x in d if x * sign > 0)
    if not steps:
        return mono, 0.0
    m = steps[len(steps) // 2]
    tol = 0.02 * m + 1e-12
    reg = sum(1 for x in d if abs(x - sign * m) <= tol) / len(d)
    return mono, reg


_AXIS_HDR = re.compile(r"^\s*(be|binding|ke|kinetic)(?![a-z])", re.I)
_ENERGY_HDR = re.compile(r"^\s*(be|binding|ke|kinetic|energy|e)(?![a-z])", re.I)
_PER_S = re.compile(r"(?<![a-z])(c\s*/\s*s|cps|(counts?|cts)\s*(/|per)\s*s)"
                    r"(?![a-z])", re.I)
_COUNTS = re.compile(r"(?<![a-z])(counts?|cts)(?![a-z])", re.I)
_ARB = re.compile(r"a\.\s?u\.?|arb", re.I)
_NAME = re.compile(r"^\s*(?:be|ke|cps|counts?|cts|intensity|int|c/s)\s*[_:]\s*"
                   r"(.+?)\s*$", re.I)


def header_info(label):
    """What a column header says: ``axis`` ("BE" / "KE" / None), ``energy``
    (it reads as an energy), ``units`` (one of ``UNITS`` or None) and ``name``
    (the region, as in ``BE_Cl2p`` or ``CPS_C 1s``; "" when absent)."""
    text = (label or "").strip()
    m = _AXIS_HDR.match(text)
    axis = None
    if m:
        axis = "BE" if m.group(1).lower() in ("be", "binding") else "KE"
    units = None
    if _PER_S.search(text):
        units = "counts/s"
    elif _COUNTS.search(text):
        units = "counts"
    elif _ARB.search(text):
        units = "a.u."
    nm = _NAME.match(text)
    return {"axis": axis, "energy": bool(_ENERGY_HDR.match(text)),
            "units": units, "name": nm.group(1) if nm else ""}


def find_energy_col(table):
    """``(column, is_axis)``: the column that runs steadily one way (preferring
    one whose header reads as an energy), and whether any column did."""
    best, best_score, ok = 0, -1.0, False
    for c, vals in enumerate(table.columns):
        mono, reg = axis_quality(vals)
        info = header_info(table.header[c]) if table.header else {}
        score = 0.5 * mono + 0.5 * reg
        if info.get("energy"):
            score += 0.3
        elif info.get("units"):
            score -= 0.3
        if score > best_score:
            best, best_score = c, score
            ok = mono >= 0.99 and reg >= REGULAR
    return best, ok


# -- naming ---------------------------------------------------------------------
def region_name(name, energy=None):
    """The region name as the rest of the app spells it (``C1s`` → ``C 1s``;
    survey aliases such as ``Gen`` → ``Survey``; a blank name on a very wide
    axis → ``Survey``)."""
    from readers.base import canon_region_name, is_survey_span
    name = (name or "").strip()
    if name.lower() in SURVEY_ALIASES:
        return "Survey"
    if not name:
        return "Survey" if energy and is_survey_span(energy) else ""
    return canon_region_name(name)


def suggest_name(be_lo, be_hi, hv=None):
    """A core-level name for a binding-energy window when exactly one common
    element's strongest line lies inside it, else ``""`` (nothing is guessed
    from an ambiguous window; rare elements are never suggested, a Cl 2s window
    would read "Fr 4f"). Meant as an editable suggestion, never applied
    silently."""
    try:
        import xpslines
        lines = xpslines.load_lines()
    except Exception:                       # noqa: BLE001
        return ""
    lo, hi = min(be_lo, be_hi), max(be_lo, be_hi)
    if hi - lo > MAX_SUGGEST_SPAN or not lines:
        return ""
    centre, half = (lo + hi) / 2, (hi - lo) / 2
    found = []
    for _d, e in xpslines.candidates(centre, half, lines, hv, split=False):
        if (e.get("rank", 1) == 1 and e["el"] in xpslines.COMMON
                and lo <= xpslines.line_be(e, hv) <= hi):
            found.append(e)
    labels = {xpslines.label_of(e) for e in found}
    return labels.pop() if len(labels) == 1 else ""


# -- options ---------------------------------------------------------------------
def guess_options(table, filename=""):
    """The options a file reads with unless the user says otherwise."""
    ecol, _ok = find_energy_col(table)
    hdr = table.header
    einfo = header_info(hdr[ecol]) if hdr else {}
    icols = [c for c in range(table.n_cols) if c != ecol]
    names, units = [], None
    for c in icols:
        info = header_info(hdr[c]) if hdr else {}
        nm = info.get("name") or einfo.get("name") or ""
        if not nm and hdr and len(icols) > 1:
            nm = hdr[c].strip()
        names.append(nm)
        units = units or info.get("units")
    return sanitise_options({
        "energy_col": ecol, "intensity_cols": icols,
        "axis": einfo.get("axis") or "BE",
        "units": units or einfo.get("units") or "a.u.",
        "names": names})


def _opt_float(v, lo=0.0):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) and x > lo else None


def sanitise_options(opts, n_cols=None):
    """A full, valid options dict from anything (unknown keys dropped, bad
    values replaced by the default). ``n_cols`` (when known) bounds the column
    numbers: an intensity column outside the table, or equal to the energy
    column, is dropped."""
    o = opts if isinstance(opts, dict) else {}
    out = {"energy_col": 0, "intensity_cols": [1], "axis": "BE",
           "units": "a.u.", "names": [], "sample": "", "photon_energy": None,
           "anode": "", "pass_energy": None}
    try:
        e = int(o.get("energy_col", 0))
    except (TypeError, ValueError):
        e = 0
    e = max(0, e)
    if n_cols is not None and e >= n_cols:
        e = 0
    out["energy_col"] = e
    cols, names = [], []
    raw_cols = o.get("intensity_cols")
    raw_names = o.get("names")
    for k, c in enumerate(raw_cols if isinstance(raw_cols, list) else [1]):
        try:
            c = int(c)
        except (TypeError, ValueError):
            continue
        if c < 0 or c == e or c in cols or (n_cols is not None
                                            and c >= n_cols):
            continue
        cols.append(c)
        nm = raw_names[k] if isinstance(raw_names, list) and k < len(
            raw_names) else ""
        names.append(str(nm).strip()[:80] if isinstance(nm, str) else "")
    if not cols:
        cols = [c for c in (range(n_cols) if n_cols else [1]) if c != e][:1]
        names = [""] * len(cols)
    out["intensity_cols"], out["names"] = cols, names
    if o.get("axis") in AXES:
        out["axis"] = o["axis"]
    if o.get("units") in UNITS:
        out["units"] = o["units"]
    if isinstance(o.get("sample"), str):
        out["sample"] = o["sample"].strip()[:120]
    out["photon_energy"] = _opt_float(o.get("photon_energy"))
    out["pass_energy"] = _opt_float(o.get("pass_energy"))
    if isinstance(o.get("anode"), str) and o["anode"] in ANODES:
        out["anode"] = o["anode"]
    return out


def photon_label(hv):
    """A short text for a photon energy: "Al Kα (1486.6 eV)" for a known
    anode, else "<hv> eV", else "unknown"."""
    if not hv:
        return "unknown"
    for name, e in ANODES.items():
        if abs(e - hv) < 0.05:
            return f"{name} ({e:g} eV)"
    return f"{hv:g} eV"


def parse_photon(text):
    """A photon energy in eV from what a user typed ("Al Kα (1486.6 eV)",
    "1486.6", "Mg"), or ``None`` for blank / unknown / not a number."""
    t = (text or "").strip()
    if not t or t.lower() in ("unknown", "not set", "(not set)") \
            or t.startswith("-"):
        return None
    for name, e in ANODES.items():
        if t.lower().startswith(name.lower()[:2]) and not re.search(r"\d", t):
            return e
        if name.lower() in t.lower():
            return e
    m = re.search(r"\d+(?:\.\d+)?", t.replace(",", "."))
    return _opt_float(m.group(0)) if m else None
