"""Element identification from a survey: an approximate table of XPS lines
(``assets/xps_lines.json``, editable), candidate lookup near a binding energy
and automatic peak labelling. Pure functions; no Tk."""

from __future__ import annotations

import json
import os
import re

DEFAULT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "assets", "xps_lines.json")
DEFAULT_HV = 1486.6                      # Al K-alpha
# Elements met on almost every sample: they win close calls (~1 eV) against
# rarer elements when labelling.
COMMON = {"C", "O", "N", "Si", "Na", "Cl", "S", "Ca", "Al", "F", "K", "P",
          "Fe", "Zn", "Cu", "Ti", "Mg", "B"}
COMMON_BONUS = 1.0
RARE_SECONDARY_PENALTY = 1.5
# A photoelectron line only exists when the photon can reach it: a line whose
# binding energy is within this many eV of hν (or beyond it) is not offered.
# Mg Kα (1253.6) so loses a 1 300 eV level that Al Kα (1486.6) shows, and the
# deep levels of the table appear only for Ag Lα and harder sources.
REACH_MARGIN = 5.0
_SPIN_ORBIT = re.compile(r"^(\d[spdf])(\d/2)$")        # 2p3/2, 3d5/2, 4f7/2


def load_lines(path=None):
    """The line table as a list of dicts (``el``, ``line``, ``be`` or
    ``ke``, ``rank``); [] if the file is missing or unreadable."""
    try:
        with open(path or DEFAULT_PATH, encoding="utf-8") as fh:
            data = json.load(fh)
        out = []
        for e in data.get("lines", []):
            if not isinstance(e, dict):
                continue
            if e.get("el") and e.get("line") and (
                    isinstance(e.get("be"), (int, float))
                    or isinstance(e.get("ke"), (int, float))):
                out.append(e)
        return out
    except (OSError, ValueError, AttributeError):
        return []


def line_be(entry, hv=None):
    """Binding energy of a table entry (Auger lines follow the photon
    energy: BE = hv - KE)."""
    if isinstance(entry.get("be"), (int, float)):
        return float(entry["be"])
    return float(hv or DEFAULT_HV) - float(entry["ke"])


def label_of(entry):
    return f"{entry['el']} {entry['line']}"


def split_line(line):
    """``("2p", "3/2")`` for a spin-orbit component such as ``2p3/2``;
    ``(line, "")`` for anything else (``1s``, ``3p``, an Auger line)."""
    m = _SPIN_ORBIT.match(line or "")
    return (m.group(1), m.group(2)) if m else (line, "")


def base_label(label):
    """A label without its spin-orbit component: ``"Ti 2p3/2"`` and
    ``"Ti 2p"`` both give ``"Ti 2p"`` (what a marker named either way means)."""
    head, _sep, tail = str(label).rpartition(" ")
    return f"{head} {split_line(tail)[0]}" if head else str(label)


def reachable(entry, hv=None):
    """False for a photoelectron line the photon energy cannot excite."""
    if "be" not in entry:
        return True                      # Auger lines follow hν through ke
    return float(entry["be"]) < float(hv or DEFAULT_HV) - REACH_MARGIN


def candidates(be, window, lines, hv=None, split=False):
    """Lines within ``window`` eV of ``be``: ``[(delta, entry)]`` with the
    most plausible first: nearest, with secondary lines (rank > 1) needing to
    be about 0.8 eV closer per rank step to win, and common elements
    (``COMMON``, main lines only) given a 1 eV head start. Lines the photon
    energy cannot reach are left out.

    The two components of a spin-orbit pair are **one candidate by default**
    (``Ti 2p``, not ``Ti 2p3/2`` and ``Ti 2p1/2``): the nearest component
    gives the delta and binding energy, the entry carries ``line`` = the pair's
    name, ``component`` = the one that matched, and the pair is ranked by its
    stronger component (a click on the 2p1/2 peak still finds ``Ti 2p``
    first). ``split=True`` lists every component on its own."""
    out = []
    for e in lines:
        if not reachable(e, hv):
            continue
        d = line_be(e, hv) - be
        if abs(d) <= window:
            out.append((d, e))
    if not split:
        out = _merge_doublets(out, lines)
    out.sort(key=_plausibility)
    return out


def _merge_doublets(found, lines):
    """``found`` with each element's spin-orbit components (and an unsplit
    entry of the same orbital) merged into the one nearest to the peak, in
    the table order of the entry kept."""
    best_rank = {}
    for e in lines:
        if "be" in e:
            k = (e["el"], split_line(e["line"])[0])
            best_rank[k] = min(best_rank.get(k, 99), e.get("rank", 1))
    nearest = {}
    for i, (d, e) in enumerate(found):
        if "be" in e:
            k = (e["el"], split_line(e["line"])[0])
            if k not in nearest or abs(d) < abs(found[nearest[k]][0]):
                nearest[k] = i
    keep = set(nearest.values())
    out = []
    for i, (d, e) in enumerate(found):
        if "be" in e:
            if i not in keep:
                continue
            k = (e["el"], split_line(e["line"])[0])
            if k[1] != e["line"]:
                e = dict(e, line=k[1], component=e["line"], rank=best_rank[k])
            elif best_rank[k] != e.get("rank", 1):
                e = dict(e, rank=best_rank[k])
        out.append((d, e))
    return out


def _plausibility(item):
    """Sort key of a ``(delta, entry)`` candidate (smaller wins): distance,
    plus 0.8 eV per rank step, minus the head start of a common element's
    main line, plus a penalty on a rarer element's secondary line (a rank-2
    line of an element nobody expects must be clearly closer to beat a main
    line of a common one: an O 1s peak 2 eV off no longer reads "At 4d3/2")."""
    d, e = item
    rank = e.get("rank", 1)
    common = e["el"] in COMMON
    key = abs(d) + 0.8 * (rank - 1)
    if common and rank == 1:
        key -= COMMON_BONUS
    elif not common and rank > 1:
        key += RARE_SECONDARY_PENALTY
    return key


def find_peaks(energy, counts, prominence=0.04, min_sep=3.0, smooth=5,
               limit=30):
    """Peak positions of a spectrum, strongest first: ``[(be, height)]``.

    A peak is a local maximum of the lightly smoothed data standing above the
    lowest point within +/- ``min_sep`` * 4 eV by at least ``prominence`` of
    the full range; no two peaks closer than ``min_sep`` eV are kept."""
    pts = sorted(zip(energy, counts))
    if len(pts) < 5:
        return []
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    n = len(ys)
    half = max(0, smooth // 2)
    sm = []
    for i in range(n):
        a, b = max(0, i - half), min(n, i + half + 1)
        sm.append(sum(ys[a:b]) / (b - a))
    span = max(sm) - min(sm)
    if span <= 0:
        return []
    reach = min_sep * 4
    found = []
    for i in range(1, n - 1):
        if sm[i] < sm[i - 1] or sm[i] <= sm[i + 1]:
            continue
        lo = min(sm[j] for j in range(n) if abs(xs[j] - xs[i]) <= reach)
        if sm[i] - lo >= prominence * span:
            found.append((xs[i], sm[i], sm[i] - lo))
    found.sort(key=lambda t: -t[2])
    kept = []
    for x, y, _p in found:
        if all(abs(x - k[0]) >= min_sep for k in kept):
            kept.append((x, y))
        if len(kept) >= limit:
            break
    return kept


def auto_label(energy, counts, lines, hv=None, window=2.0, split=False, **kw):
    """Peaks of a survey each given the best-matching line:
    ``[(be, "C 1s")]`` (peaks with no candidate are left out). Unless
    ``split``, a doublet is named once (``Ti 2p``) at its strongest peak, not
    at both of its components."""
    out, seen = [], set()
    for be, _h in find_peaks(energy, counts, **kw):
        cand = candidates(be, window, lines, hv, split=split)
        if cand:
            label = label_of(cand[0][1])
            if not split and label in seen:
                continue                 # strongest peak first: keep that one
            seen.add(label)
            out.append((be, label))
    return sorted(out)


def nearby_lines(be, window, lines, hv=None, exclude=None, secondary=False,
                 auger=False, max_extra=2, split=False):
    """Other candidate lines near ``be`` besides the primary match
    (``exclude``, its label from :func:`label_of`; a spin-orbit pair counts
    as one, so ``Ti 2p`` and ``Ti 2p3/2`` both leave out the whole doublet),
    for showing alongside an already-identified peak: up to ``max_extra``
    other photoelectron lines (nearest first) when ``secondary`` is set, and
    every Auger line in the window when ``auger`` is set. ``[(be, "El line",
    "secondary"|"auger")]``, each tier nearest first."""
    skip = base_label(exclude) if exclude else None
    extra_secondary, extra_auger = [], []
    for _d, e in candidates(be, window, lines, hv, split=split):
        label = label_of(e)
        if skip is not None and base_label(label) == skip:
            continue
        if "ke" in e:
            if auger:
                extra_auger.append((line_be(e, hv), label, "auger"))
        elif secondary:
            extra_secondary.append((line_be(e, hv), label, "secondary"))
    return extra_secondary[:max_extra] + extra_auger
