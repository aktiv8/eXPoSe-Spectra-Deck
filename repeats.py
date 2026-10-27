"""Repeat groups: the mean, standard deviation and number of the atomic percent
of a region (and of each chemical state of it) over samples that are repeats of
one another -- several locations on one sample, a spot measured again, pieces of
one material. Tk-free, no numpy.

A group is the user's statement (``Annotations.repeat_groups``): a name, a kind
and the keys of the samples in it (``resultspages.Sample.key``). The statistics
are worked out from each member's own normalised composition, so they follow
the same rules as every other table (a region counted once, the user's ticks,
CasaXPS's own numbers when those are in use): nothing is requantified here.

The layout follows the table XPSView / HarwellXPS writes with its reports
(``quantification_repeat_groups.csv``: group, kind, region, state, mean, sd
(n - 1), n, individual values). A region or state a member does not have is
absent from that member and ``n`` says how many members had it; a mean over one
value has no standard deviation. Depth profiles are not repeats and are left out.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import quant

KINDS = {"locations": "different locations on the sample",
         "repeats": "repeat measurements of one spot",
         "samples": "separate pieces of one material"}
DEFAULT_KIND = "locations"

HEADER = ("Group", "Kind", "Region", "State", "at % mean", "at % sd (n-1)",
          "n", "Individual values")


def sanitise_group(g):
    """``{"kind", "samples"}`` from a saved group, or None when it names no
    sample. A duplicate sample key is kept once, in order."""
    if not isinstance(g, dict):
        return None
    keys = g.get("samples")
    if not isinstance(keys, (list, tuple)):
        return None
    seen = []
    for k in keys:
        if isinstance(k, str) and k and k not in seen:
            seen.append(k)
    if not seen:
        return None
    kind = g.get("kind")
    return {"kind": kind if kind in KINDS else DEFAULT_KIND, "samples": seen}


def sanitise_groups(d):
    """A clean ``{name: group}`` from saved data (bad entries dropped)."""
    out = {}
    if isinstance(d, dict):
        for name, g in d.items():
            clean = sanitise_group(g)
            if clean and str(name).strip():
                out[str(name).strip()] = clean
    return out


def stats(values):
    """``(mean, sd, n)``; the sample standard deviation (n - 1), None below
    two values. ``(None, None, 0)`` for none."""
    vals = [float(v) for v in values if v is not None and math.isfinite(v)]
    n = len(vals)
    if not n:
        return None, None, 0
    mean = sum(vals) / n
    if n < 2:
        return mean, None, n
    return mean, math.sqrt(sum((v - mean) ** 2 for v in vals) / (n - 1)), n


@dataclass
class Row:
    region: str
    state: str = ""
    mean: float | None = None
    sd: float | None = None
    n: int = 0
    of: int = 0                                 # members in the group
    values: list = field(default_factory=list)  # [(member label, at %)]


@dataclass
class Group:
    name: str
    kind: str = DEFAULT_KIND
    members: list = field(default_factory=list)   # labels counted
    rows: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def kind_text(self):
        return KINDS.get(self.kind, self.kind)


def member_values(sample):
    """``[((region, state), at %)]`` for a one-level sample, in table order:
    each region that counts, then its states (a single state is listed too:
    ``compute`` drops them when no member has more than one). A region listed
    twice keeps its first value."""
    if not sample.levels:
        return []
    level = sample.levels[0]
    out, seen = [], set()
    for e, x in zip(level.entries, level.res):
        pct = x.get("at_pct")
        row = e["row"]
        name = row["region"]
        if pct is None or (name, "") in seen:
            continue
        seen.add((name, ""))
        out.append(((name, ""), pct))
        out += [((name, st["name"]), st["at_pct"])
                for st in quant.states(row, pct)]
    return out


def compute(samples, groups):
    """``[Group]`` for the saved ``groups`` ({name: {"kind", "samples"}}) over
    the ``resultspages.Sample`` list. A key that is not among the samples, a
    depth profile and a member with nothing quantified are left out and a note
    says which; a group with no member left has no rows."""
    by_key = {s.key: s for s in samples}
    out = []
    for name, g in groups.items():
        grp = Group(name, g.get("kind", DEFAULT_KIND))
        cols = []
        for key in g.get("samples", []):
            s = by_key.get(key)
            if s is None:
                grp.notes.append(f"{key}: not among the quantified samples.")
            elif s.is_profile:
                grp.notes.append(
                    f"{s.label}: a depth profile is not a repeat; left out.")
            else:
                vals = member_values(s)
                if vals:
                    cols.append((s.label, vals))
                    grp.members.append(s.label)
                else:
                    grp.notes.append(
                        f"{s.label}: nothing quantified; left out.")
        order, per = [], {}
        for label, vals in cols:
            for k, v in vals:
                if k not in per:
                    per[k] = []
                    order.append(k)
                per[k].append((label, v))
        names = {}
        for reg, st in order:
            if st:
                names.setdefault(reg, set()).add(st)
        # a region with one state in every member says it already: no rows
        order = [k for k in order if not k[1] or len(names[k[0]]) > 1]
        for k in order:
            mean, sd, n = stats(v for _, v in per[k])
            grp.rows.append(Row(k[0], k[1], mean, sd, n, len(cols), per[k]))
        if 0 < len(cols) < 2:
            grp.notes.append("Only one member: there is no spread to report.")
        out.append(grp)
    return out


MAX_MEMBER_COLUMNS = 6


def table(group, max_members=MAX_MEMBER_COLUMNS):
    """``(header, rows)`` of one repeat group for display: ``rows`` are
    ``(kind, cells, partial)`` with ``kind`` "region" or "state" and
    ``partial`` true when fewer members than the group has contributed (the
    row's ``n`` is below its size). Each member gets its own column up to
    ``max_members``; with more, the tab lists the values in one column
    (``max_members`` None) and a report leaves them out (``max_members`` 0
    gives no member columns at all)."""
    heads = short_labels(group.members)
    wide = len(group.members) <= max_members
    header = ["Region", "State", "Mean at %", "SD (n-1)", "n"]
    header += heads if wide else []
    rows = []
    for r in group.rows:
        by = dict(r.values)
        cells = [("    " if r.state else "") + r.region, r.state,
                 f"{r.mean:.2f}", "–" if r.sd is None else f"{r.sd:.2f}",
                 str(r.n)]
        if wide:
            cells += [f"{by[m]:.2f}" if m in by else "–"
                      for m in group.members]
        rows.append(("state" if r.state else "region", cells, r.n < r.of))
    return header, rows


def _g(v, nd=3):
    return "" if v is None else f"{v:.{nd}f}"


def csv_rows(groups):
    """The table as rows of cells, header first (XPSView's layout)."""
    out = [list(HEADER)]
    for g in groups:
        for r in g.rows:
            out.append([g.name, g.kind_text, r.region, r.state, _g(r.mean),
                        _g(r.sd), str(r.n),
                        "; ".join(f"{lab}={v:.3f}" for lab, v in r.values)])
    return out


_TAIL = re.compile(
    r"(\s*<[^>]*>|\s+(?:location|spot|point|pt|repeat|rep)\s*#?\d+\w?)\s*$", re.I)


def base_label(label):
    """A sample label without what tells repeats of it apart: the text from
    its first " - " on and a trailing "<Pt #001a>", "location 2" or
    "repeat 3"."""
    s = re.split(r"\s+[-–—]\s+", label or "", maxsplit=1)[0].strip()
    prev = None
    while prev != s:
        prev, s = s, _TAIL.sub("", s).strip()
    return s


def suggest(samples):
    """Groups the labels suggest: ``{name: {"kind", "samples"}}`` for each
    base label shared by two or more one-level samples (survey scans are
    grouped apart from the high-resolution regions, as they are quantified
    apart). A starting point for the user, never applied by itself."""
    found = {}
    for s in samples:
        if s.is_profile or s.kind not in ("regions", "survey") or not s.levels:
            continue
        label = re.sub(r"\s*\(survey\)\s*$", "", s.label)
        name = base_label(label) + (" (survey)" if s.kind == "survey" else "")
        found.setdefault(name, []).append(s.key)
    return {n: {"kind": DEFAULT_KIND, "samples": keys}
            for n, keys in found.items() if n and len(keys) > 1}


def short_labels(labels, limit=26):
    """Column heads for the members of a group: each label without the words
    every label has (the sample name, "Spot sequence", a shared date), the
    separators that are left over trimmed. A label with nothing left keeps its
    last word; one longer than ``limit`` characters is cut."""
    labels = list(labels)
    words = [lab.split() for lab in labels]
    common = set(words[0]) if words else set()
    for w in words[1:]:
        common &= set(w)
    out = []
    for lab, w in zip(labels, words):
        kept = [x for x in w if x not in common]
        t = " ".join(kept).strip(" -–—,") if len(labels) > 1 else lab
        if not t:
            t = w[-1] if w else lab
        out.append(t if len(t) <= limit else t[:limit - 1] + "…")
    return out
