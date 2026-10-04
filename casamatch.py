"""Tie CasaXPS's own exported quantification (``casaquant.py``) to the fitted
regions the user ticked (Tk-free, nothing guessed).

``Quant_regions.txt`` lists one row per fitted component of every region
CasaXPS had ticked (names such as "Cl 2p x" or "Pt 4f loss", a name may repeat,
only a position tells rows apart); ``Quant_survey.txt`` lists one row per
element of the survey scan. This module gives each ticked fit region its rows
and tags the region's ``quant.fit_rows`` row dict so ``quant.normalise`` uses
CasaXPS's number instead of recomputing one:

* ``row["casa_pct"]`` -- the CasaXPS %At of the region (the sum of its rows) and
  ``row["casa_rows"]`` -- the file rows behind it, or
* ``row["casa_why"]`` -- why the region has none ("no CasaXPS quantification
  found for this region"). Such a region is left out of the total rather than
  filled with a recomputed value: CasaXPS's %At carries transmission and
  mean-free-path terms ours does not, so the two scales are never mixed.

Percentages are renormalised over the regions that are counted
(``at% = casa_pct / sum``), which keeps CasaXPS's ratios.
"""

from __future__ import annotations

POSITION_TOL = 0.5          # eV: component positions vs. file positions
LIMIT_MARGIN = 1.0          # eV: a file row this far outside a fit region's
                            # limits does not belong to it

WHY_NONE = "no CasaXPS quantification found for this region"
WHY_NO_SURVEY = "no CasaXPS survey quantification found for this region"


def norm(text):
    """A name with its spaces collapsed and case folded ("Cl 2p " -> "cl 2p")."""
    return " ".join(str(text or "").split()).lower()


def region_of(file_name, region_names):
    """The region a file row's name belongs to: the longest of
    ``region_names`` the name starts with, as a whole word ("Cl 2p x" and
    "Cl 2p #" -> "Cl 2p", "Pt 4f loss" -> "Pt 4f"), else None."""
    n = norm(file_name)
    best = None
    for r in region_names:
        rn = norm(r)
        if rn and (n == rn or n.startswith(rn + " ")):
            if best is None or len(rn) > len(norm(best)):
                best = r
    return best


def blocks(rows, region_names, name_key="name"):
    """File rows as ``[(region, [row, ...])]``: consecutive rows of one region
    form one block, so a region listed again after another region (a second
    acquisition of the sample) starts a new one. Rows that belong to none of
    ``region_names`` are dropped."""
    out = []
    for row in rows:
        region = region_of(row.get(name_key), region_names)
        if region is None:
            continue
        if out and norm(out[-1][0]) == norm(region):
            out[-1][1].append(row)
        else:
            out.append((region, [row]))
    return out


def _cost(entry_row, file_rows):
    """Mean distance (eV) from each file row's position to the nearest fitted
    component of the region, or None when the region has no components."""
    comps = entry_row.get("components") or []
    if not comps or not file_rows:
        return None
    pos = [r.get("position") for r in file_rows if r.get("position") is not None]
    if not pos:
        return None
    return sum(min(abs(p - c["be"]) for c in comps) for p in pos) / len(pos)


def _inside(entry_row, file_rows):
    """False when a file row lies clearly outside the fit region's limits."""
    lo, hi = entry_row.get("be_lo"), entry_row.get("be_hi")
    if lo is None or hi is None:
        return True
    return all(lo - LIMIT_MARGIN <= r["position"] <= hi + LIMIT_MARGIN
               for r in file_rows if r.get("position") is not None)


def _tag(row, rows):
    pct = sum(r["at_pct"] for r in rows if r.get("at_pct") is not None)
    row.pop("casa_why", None)
    row["casa_pct"] = pct
    row["casa_rows"] = [dict(r) for r in rows]


def _untag(row, why):
    row.pop("casa_pct", None)
    row.pop("casa_rows", None)
    row["casa_why"] = why


def tag_regions(entries, file_rows, label=""):
    """Tag the high-resolution ``entries`` (``[{"spectrum", "row"}]`` of one
    sample) from that sample's ``Quant_regions.txt`` rows (``[{"name",
    "position", "at_pct"}]``). Returns the notes to show. With no file rows the
    entries are left untouched (the recomputed numbers stand)."""
    if not file_rows:
        return []
    notes = []
    names = sorted({e["row"]["region"] for e in entries})
    by_region = {}
    for region, rows in blocks(file_rows, names):
        by_region.setdefault(norm(region), []).append(rows)
    by_entry = {}
    for e in entries:
        by_entry.setdefault(norm(e["row"]["region"]), []).append(e)
    for key, ents in by_entry.items():
        bl = by_region.get(key, [])
        pairs = _assign(ents, bl)
        for e, rows in pairs:
            row = e["row"]
            if rows is None:
                _untag(row, WHY_NONE)
                continue
            if not _inside(row, rows):
                _untag(row, "the CasaXPS rows for this name lie outside the "
                            "fitted region")
                continue
            _tag(row, rows)
            comps = len(row.get("components") or [])
            if comps and comps != len(rows):
                notes.append(
                    f"CasaXPS's file lists {len(rows)} component(s) for "
                    f"{row['region']}" + (f" in {label}" if label else "")
                    + f", the fit here has {comps}: the file may come from a "
                    "different fit.")
    return notes


def _assign(ents, bl):
    """``[(entry, rows or None)]``: blocks of file rows paired with entries of
    one region name. In file order when the counts agree; otherwise by how well
    the component positions agree (a block goes to one entry only)."""
    if not bl:
        return [(e, None) for e in ents]
    if len(bl) == len(ents):
        return list(zip(ents, bl))
    scored = []
    for i, e in enumerate(ents):
        for j, rows in enumerate(bl):
            c = _cost(e["row"], rows)
            scored.append((c if c is not None else POSITION_TOL * 0.99, i, j,
                           c is None))
    scored.sort()
    used_e, used_b, got = set(), set(), {}
    for c, i, j, blind in scored:
        if i in used_e or j in used_b:
            continue
        if c > POSITION_TOL and not (blind and len(ents) == 1
                                     and len(bl) == 1):
            continue
        used_e.add(i)
        used_b.add(j)
        got[i] = bl[j]
    return [(e, got.get(i)) for i, e in enumerate(ents)]


def tag_survey(entries, survey_rows):
    """Tag the survey ``entries`` of one sample from its ``Quant_survey.txt``
    rows (``[{"element", "pct"}]``): an element's k-th listing goes to the k-th
    survey fit region of that name. With no file rows the entries are left
    untouched."""
    if not survey_rows:
        return []
    listed = {}
    for r in survey_rows:
        listed.setdefault(norm(r.get("element")), []).append(r)
    seen = {}
    for e in entries:
        row = e["row"]
        key = norm(row["region"])
        k = seen[key] = seen.get(key, -1) + 1
        rows = listed.get(key, [])
        if k < len(rows) and rows[k].get("pct") is not None:
            row.pop("casa_why", None)
            row["casa_pct"] = rows[k]["pct"]
            row["casa_rows"] = [{"name": rows[k]["element"], "position": None,
                                 "at_pct": rows[k]["pct"]}]
        else:
            _untag(row, WHY_NO_SURVEY)
    return []
