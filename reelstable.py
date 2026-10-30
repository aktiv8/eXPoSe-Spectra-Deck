"""The REELS π / (σ+π) ratios of a workbook as one table, without Tk: the rows of
the "REELS analysis" section of the PDF, the slides and the Word document, and
the CSV of the hand-over package and *Export REELS ratios*.

Nothing is computed here that ``reels`` does not already: a row is a stored ratio
(``Annotations.pi_ratio``), its sp² fraction comes from ``reels.sp2_for`` against
the workbook's HOPG / diamond references, and what is wrong with it comes from
``reels.ratio_flags``. A value that is not there is left blank, never filled in
(a ratio without a σ+π area, a fraction with no reference), and a number is
never clipped. The CSV follows the layout of HarwellXPS's
``reels_pi_sigma_ratio.csv`` (sample, elastic peak, π and σ+π positions, the two
areas, the ratio, the windows) with the file, the sp² fraction, the end-point
channels, the two references and the notes added.
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field

import reels
import repeats

DASH = "–"

# What the section says under its table: the one definition, and the limits of
# the number (the last sentence is HarwellXPS's own, in its REELS card).
CAVEAT = ("The ratio is not an sp² percentage: it depends on the windows, the "
          "beam energy and the geometry, so on its own it compares spectra "
          "measured together. The sp² fraction is relative to the HOPG (= 1) "
          "and, if one was measured, diamond (= 0) of this workbook, and is "
          "only given where they were measured with the same windows.")


@dataclass
class Row:
    """One spectrum's ratio, as numbers (``None`` where there is none)."""
    file: str
    sample: str
    spectrum: str
    elastic: float | None
    pos_pi: float | None
    pos_sp: float | None
    i_pi: float
    i_sp: float
    ratio: float | None
    fraction: float | None
    pi_window: list
    sp_window: list
    avg: int
    key: str = ""                                  # the repeat-group key of its sample
    flags: list = field(default_factory=list)      # what is wrong with it
    problem: str = ""                              # why there is no fraction
    notes: list = field(default_factory=list)      # about the fraction


@dataclass
class Table:
    rows: list
    refs: dict                       # the workbook's {"hopg": ..., "diamond": ...}
    unit: str = ""                   # the intensity unit of the areas ("" = mixed)
    settings: list = field(default_factory=list)   # distinct window sentences
    groups: list = field(default_factory=list)     # [RepGroup]: the repeat statistics

    def __bool__(self):
        return bool(self.rows)

    def __len__(self):
        return len(self.rows)


def collect(items, refs=None, groups=None):
    """The table of ``items`` = ``[{file, sample, spectrum, unit, res, key}]``,
    where ``res`` is a stored ratio result (``reels.pi_ratio`` plus
    ``elastic``) and ``key`` the repeat-group key of the sample
    (``<file key>/<sample>``, as ``resultspages.Sample.key``); ``refs`` is
    ``Annotations.reels_refs`` and ``groups`` ``Annotations.repeat_groups``
    (the repeat statistics are worked out for the groups that have a REELS
    spectrum). An item with no result is skipped. Returns a :class:`Table`
    (falsy when there are no rows)."""
    refs = refs or {}
    rows, units, settings = [], set(), []
    for it in items or ():
        res = it.get("res")
        if not res:
            continue
        sp = reels.sp2_for(res, refs)
        rows.append(Row(
            file=str(it.get("file") or ""), sample=str(it.get("sample") or ""),
            spectrum=str(it.get("spectrum") or ""),
            elastic=res.get("elastic"), pos_pi=res.get("pos_pi"),
            pos_sp=res.get("pos_sp"), i_pi=res["i_pi"], i_sp=res["i_sp"],
            ratio=res.get("ratio"), fraction=sp["fraction"],
            pi_window=list(res["pi_window"]), sp_window=list(res["sp_window"]),
            avg=res["avg"], key=str(it.get("key") or ""),
            flags=reels.ratio_flags(res),
            problem=sp["problem"], notes=list(sp["notes"])))
        units.add(str(it.get("unit") or ""))
        text = reels.ratio_text(res)
        if text not in settings:
            settings.append(text)
    unit = units.pop() if len(units) == 1 else ""
    table = Table(rows, refs, unit, settings)
    table.groups = repeat_groups(table, groups)
    return table


# -- repeat groups ------------------------------------------------------------------------
# The mean, sample standard deviation (n - 1) and n of the ratio, the sp2 fraction
# and the two peak positions over the spectra of samples that repeat one another.
# A group is the user's statement (``Annotations.repeat_groups``, the same groups
# as the Quantification tab's *Groups…*: a name, a kind and sample keys); nothing
# is re-measured here, each value is the one in the table.
QUANTITIES = (("Iπ/Iσ+π", "ratio", ".4f"),
              ("sp² fraction", "fraction", ".2f"),
              ("π position (eV)", "pos_pi", ".2f"),
              ("σ+π position (eV)", "pos_sp", ".2f"))


@dataclass
class QRow:
    """One quantity of one group: its statistics and each spectrum's value."""
    name: str
    spec: str                                    # number format
    mean: float | None
    sd: float | None
    n: int
    values: list = field(default_factory=list)   # [(spectrum label, value)]


@dataclass
class RepGroup:
    name: str
    kind: str
    members: list = field(default_factory=list)  # spectrum labels counted
    rows: list = field(default_factory=list)     # [QRow]
    notes: list = field(default_factory=list)

    @property
    def kind_text(self):
        return repeats.KINDS.get(self.kind, self.kind)


def _setting(r):
    return (tuple(r.pi_window), tuple(r.sp_window), r.avg)


def repeat_groups(table, groups):
    """``[RepGroup]`` for the saved ``groups`` over the rows of ``table``. A
    group none of whose samples has a REELS spectrum is about something else
    (the atomic % of the Quantification tab) and is left out; a member with no
    ratio is named and left out. A group whose spectra were measured with other
    windows or end-point channels has no statistics (the ratios are not
    comparable), and a flagged spectrum is counted as measured and named."""
    out = []
    by_key = {}
    for r in table.rows:
        by_key.setdefault(r.key, []).append(r)
    for name, g in (groups or {}).items():
        keys = list(g.get("samples", []))
        if not any(k in by_key for k in keys):
            continue
        grp = RepGroup(name, g.get("kind", repeats.DEFAULT_KIND))
        chosen = []
        for k in keys:
            rows = by_key.get(k)
            if not rows:
                grp.notes.append(f"{k.split('/', 1)[-1]}: no REELS ratio "
                                 f"measured; left out.")
                continue
            for r in rows:
                label = r.sample if len(rows) == 1 else (
                    f"{r.sample} / {r.spectrum}")
                chosen.append((label, r))
        grp.members = [lab for lab, _r in chosen]
        if len({_setting(r) for _l, r in chosen}) > 1:
            grp.notes.append(
                "The spectra were measured with different windows or "
                "end-point channels, so their ratios are not comparable: no "
                "statistics.")
            out.append(grp)
            continue
        for qname, attr, spec in QUANTITIES:
            vals = [(lab, getattr(r, attr)) for lab, r in chosen
                    if getattr(r, attr) is not None]
            mean, sd, n = repeats.stats(v for _l, v in vals)
            if n:
                grp.rows.append(QRow(qname, spec, mean, sd, n, vals))
        if len(chosen) < 2:
            grp.notes.append("Only one spectrum: there is no spread to "
                             "report.")
        marked = [lab for lab, r in chosen if r.flags or r.problem or r.notes]
        if marked:
            grp.notes.append("Counted as measured, but flagged in the table: "
                             + ", ".join(marked) + ".")
        out.append(grp)
    return out


def group_lines(table):
    """The groups in one line each, for the dialog: "name (n spectra): ratio
    0.0414 ± 0.0021, sp² 0.98 ± 0.05"."""
    lines = []
    for g in table.groups:
        parts = []
        for r in g.rows[:2]:
            sd = f" ± {format(r.sd, r.spec)}" if r.sd is not None else ""
            label = "ratio" if r.name.startswith("Iπ") else r.name.split()[0]
            parts.append(f"{label} {format(r.mean, r.spec)}{sd}")
        lines.append(f"{g.name} ({len(g.members)} spectra)"
                     + (": " + ", ".join(parts) if parts else
                        ": no statistics"))
    return lines


def repeat_header():
    return ["Quantity", "Mean", "SD (n−1)", "n", "Values"]


def repeat_report_rows(group):
    """One list of strings per quantity; the Values column is each spectrum's
    value, the labels shortened by ``repeats.short_labels``."""
    out = []
    for r in group.rows:
        labels = repeats.short_labels([lab for lab, _v in r.values])
        out.append([r.name, format(r.mean, r.spec),
                    DASH if r.sd is None else format(r.sd, r.spec), str(r.n),
                    "; ".join(f"{lab}={format(v, r.spec)}"
                              for lab, (_x, v) in zip(labels, r.values))])
    return out


REPEAT_WEIGHTS = [26, 17, 17, 7, 100]
REPEAT_INTRO = ("Mean, sample standard deviation (n − 1) and n of each "
                "quantity over the spectra of samples that repeat one another "
                "(the groups of the Quantification tab, Groups…). Each value "
                "is the one in the table above; nothing is re-measured. A "
                "standard deviation needs two spectra.")


def repeat_csv_header():
    return ["Group", "Kind", "Quantity", "Mean", "SD (n-1)", "n",
            "Individual values"]


def repeat_csv_rows(table):
    """The repeat statistics as CSV rows (header first), in the layout of the
    atomic % repeat-group CSV (group, kind, quantity, mean, sd, n, each
    value); an empty list when no group has statistics."""
    rows = [repeat_csv_header()]
    for g in table.groups:
        for r in g.rows:
            rows.append([g.name, g.kind_text, r.name, _g(r.mean, ".6g"),
                         _g(r.sd, ".6g"), str(r.n),
                         "; ".join(f"{lab}={format(v, '.6g')}"
                                   for lab, v in r.values)])
    return rows if len(rows) > 1 else []


def repeat_csv_bytes(table):
    """The repeat-group CSV as UTF-8 with a byte-order mark (b"" when none)."""
    rows = repeat_csv_rows(table)
    if not rows:
        return b""
    buf = io.StringIO(newline="")
    csv.writer(buf).writerows(rows)
    return buf.getvalue().encode("utf-8-sig")


def suggest_samples(table):
    """One stand-in per sample with a ratio, with the attributes
    ``repeats.suggest`` reads, so the samples of a REELS-only workbook can be
    grouped by name too."""
    from types import SimpleNamespace
    seen, out = set(), []
    for r in table.rows:
        if r.key and r.key not in seen:
            seen.add(r.key)
            out.append(SimpleNamespace(key=r.key, label=r.sample,
                                       is_profile=False, kind="regions",
                                       levels=[1]))
    return out


# -- words ------------------------------------------------------------------------------
def reference_text(table):
    """The references in a sentence (what the sp² fraction is relative to)."""
    bits = []
    for which, name in (("hopg", "HOPG"), ("diamond", "diamond")):
        ref = (table.refs or {}).get(which)
        if ref:
            how = ("typed" if ref.get("kind") == "entered" else "measured")
            src = f", {ref['source']}" if ref.get("source") else ""
            bits.append(f"{name} {ref['ratio']:.4f} "
                        f"({how}{src}) = {1 if which == 'hopg' else 0}")
    if not bits:
        return ("No HOPG reference is stored, so there is no sp² fraction "
                "(ISS / REELS ▸ REELS π / (σ+π) ratio).")
    return "sp² fraction relative to " + "; ".join(bits) + "."


_PCT = re.compile(r"(\d+(?:\.\d+)?) %")
_EV = re.compile(r"\s*\(\d+ eV\)")


def issues(table):
    """What is wrong, one line per kind of problem and the spectra it applies
    to: ``[(text, [spectrum labels])]`` in the order met. Flags that differ only
    in their numbers are one line: a percentage that varies is given as its
    range ("dips 22–378 %"), and the energies of the beam-energy note, which
    differ per spectrum, are left to the CSV."""
    found = {}
    for r in table.rows:
        label = r.sample + (f" / {r.spectrum}" if r.spectrum
                            and r.spectrum != r.sample else "")
        for text in r.flags + ([f"no sp² fraction: {r.problem}"]
                               if r.problem else []) + r.notes:
            key = _EV.sub("", _PCT.sub("N %", text))
            found.setdefault(key, []).append((text, label))
    out = []
    for entries in found.values():
        text, labels = entries[0][0], [l for _t, l in entries]
        if len(entries) > 1:
            nums = sorted({float(m.group(1)) for t, _l in entries
                           for m in [_PCT.search(t)] if m})
            if len(nums) > 1:
                rng = f"{nums[0]:g}–{nums[-1]:g} %"
                text = _PCT.sub(rng, text, count=1)
            text = _EV.sub("", text)
        out.append((text, labels))
    return out


def issue_lines(table):
    """``issues`` as sentences: "text — all N spectra" / "text — a, b"."""
    n = len(table.rows)
    lines = []
    for text, labels in issues(table):
        which = (f"all {n} spectra" if len(labels) == n and n > 1
                 else ", ".join(labels))
        lines.append(f"{text} — {which}")
    return lines


# -- the report table -------------------------------------------------------------------
def _num(v, spec):
    """``v`` for the report table: a dash for none, and no minus sign on a
    value that rounds to zero ("-0.0" is float noise, not a measurement); the
    CSV keeps the value as measured."""
    if v is None:
        return DASH
    text = format(v, spec)
    return text.lstrip("-") if float(text) == 0 else text


def header(table):
    """The report's column headings (short: the unit is in the text under the
    table, ``units_text``)."""
    return ["Sample", "Spectrum", "Elastic (eV)", "π (eV)", "σ+π (eV)",
            "Iπ", "Iσ+π", "Iπ/Iσ+π", "sp²"]


def units_text(table):
    """What the columns are, in one sentence (the areas' unit when the spectra
    agree on one)."""
    area = (f"The areas Iπ and Iσ+π are in {table.unit}·eV"
            if table.unit else "The areas Iπ and Iσ+π are in the spectra's "
                               "intensity unit times eV")
    return (f"Elastic is the kinetic energy of the elastic peak; π and σ+π "
            f"are the positions of the two peaks in eV of energy loss. "
            f"{area}.")


# column widths: PDF in mm (180 across), slides and Word as weights
WEIGHTS = [40, 24, 17, 13, 14, 15, 17, 19, 13]


def report_rows(table):
    """The report table, one list of strings per spectrum (a spectrum with a
    problem has a † after its sample's name; the problems follow the table)."""
    flagged = {id(r) for r in table.rows if r.flags or r.problem or r.notes}
    return [[r.sample + (" †" if id(r) in flagged else ""), r.spectrum,
             _num(r.elastic, ".2f"), _num(r.pos_pi, ".2f"),
             _num(r.pos_sp, ".2f"), _num(r.i_pi, ".1f"), _num(r.i_sp, ".1f"),
             _num(r.ratio, ".4f"), _num(r.fraction, ".2f")]
            for r in table.rows]


def paragraphs(table):
    """The text that goes with the table: the definition, the references, the
    caveat, then one line per kind of problem (marked † in the table)."""
    out = ["; ".join(table.settings) + "." if table.settings else "",
           units_text(table), reference_text(table), CAVEAT]
    lines = issue_lines(table)
    if lines:
        out.append("† " + " · ".join(lines))
    return [p for p in out if p]


# -- the CSV ----------------------------------------------------------------------------
def csv_header(table):
    unit = f" ({table.unit}·eV)" if table.unit else " (·eV)"
    return ["File", "Sample", "Spectrum", "Elastic peak KE (eV)",
            "pi position (eV)", "sigma+pi position (eV)",
            "I_pi" + unit, "I_sigma+pi" + unit, "I_pi/I_sigma+pi",
            "sp2 fraction", "pi window (eV)", "sigma+pi window (eV)",
            "End-point channels", "HOPG ratio", "Diamond ratio", "Notes"]


def _g(v, spec):
    return "" if v is None else format(v, spec)


def csv_rows(table):
    """The CSV as lists of strings: a header and one row per spectrum, in
    XPSView's ``reels_pi_sigma_ratio.csv`` order with the additions named in the
    module text. Numbers are written to the precision they are known to."""
    refs = table.refs or {}
    hopg = _g((refs.get("hopg") or {}).get("ratio"), ".6g")
    dia = _g((refs.get("diamond") or {}).get("ratio"), ".6g")
    rows = [csv_header(table)]
    for r in table.rows:
        rows.append([
            r.file, r.sample, r.spectrum, _g(r.elastic, ".3f"),
            _g(r.pos_pi, ".3f"), _g(r.pos_sp, ".3f"), _g(r.i_pi, ".6g"),
            _g(r.i_sp, ".6g"), _g(r.ratio, ".6g"), _g(r.fraction, ".4g"),
            f"{r.pi_window[0]:g}-{r.pi_window[1]:g}",
            f"{r.sp_window[0]:g}-{r.sp_window[1]:g}", str(r.avg), hopg, dia,
            "; ".join(r.flags + ([f"no sp2 fraction: {r.problem}"]
                                 if r.problem else []) + r.notes)])
    return rows


def csv_bytes(table):
    """The CSV as UTF-8 with a byte-order mark (Excel then reads the Greek
    letters), like the other tables the app writes."""
    buf = io.StringIO(newline="")
    csv.writer(buf).writerows(csv_rows(table))
    return buf.getvalue().encode("utf-8-sig")
