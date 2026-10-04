"""Quantification as report pages and slides (Tk-free; matplotlib only to draw
a chart).

What the pages say is what the desktop and the HTML browser say: atomic percent
from CasaXPS's own areas and sensitivity factors (``quant``), each sample at
each depth level normalised on its own. Nothing is fitted or guessed here.
:func:`collect` reads the loaded files into a :class:`Results`; the report and
the deck then lay out the same numbers (``composition_cells`` for a sample at
one level, ``profile_cells`` for a depth profile, ``profile_png`` for its
chart), so the PDF, the slides and the numbers on screen agree.

Rules for what counts, stated on the pages (``Results.method`` / ``notes``):

* the areas are the integral of data minus background over the fit region, in
  counts/s x eV, as CasaXPS reports them; the transmission function is not
  applied (it is not known whether CasaXPS applies it to every instrument);
* a region name that appears in more than one spectrum at the same position and
  depth is counted once (the dedicated scan, else the first), so the total is
  not inflated by a survey and its own high-resolution scan of the same line;
  an element counted from two of its lines (2p and 1s) is counted from both and
  the page says so -- unless one of the two is that element's known "standard"
  line (``_PREFERRED_LINE``, e.g. Pt 4f over Pt 4d), in which case only that
  one is counted and the page says so instead;
* a region with no sensitivity factor or no area is listed and says why; an
  optional reference-table fallback (``collect(..., rsf_table=)``, off by
  default -- ``quant.py``'s own "nothing is guessed" stance) can supply one
  when the file records none, and the page marks it as a substitute, never
  as though it were the file's own recorded value; when the fallback is off
  and it would have applied, a note points at the Report generator option
  instead of leaving "no RSF" looking like a dead end (``_rsf_hint_note``).
"""

from __future__ import annotations

import io
import os
import re
from dataclasses import dataclass, field

import casamatch
import casaquant
import quant
import reportspec
import rsf as rsf_lib

METHOD = ("Atomic percent is each fitted region's area divided by its "
          "sensitivity factor (RSF), as a share of the sum over the regions "
          "of the same sample and depth level. Areas are the integral of the "
          "data minus the background over the fit region, in counts/s·eV, "
          "with the RSFs and fits recorded by CasaXPS; no transmission "
          "correction is applied.")

METHOD_CASA = (METHOD + " Where CasaXPS's own quantification files "
               "(Quant_regions.txt, Quant_survey.txt) cover a sample, its "
               "percentages are used instead, shared out over the regions "
               "counted, and the sample says so; a survey scan is always "
               "its own total.")
CASA_FOOT = ("Atomic % = CasaXPS's own quantification (Quant_regions.txt / "
             "Quant_survey.txt), shared out over the regions counted; "
             "nothing recomputed.")

COMPOSITION_HEADER = ("Region", "Background", "RSF", "Area (counts/s·eV)",
                      "Area / RSF", "at %", "Fit RMS", "Reduced χ²")
CASA_HEADER = tuple("CasaXPS %At" if h == "Area / RSF" else h
                    for h in COMPOSITION_HEADER)
PROFILE_AXES = (("depth", "Depth (nm)"), ("etch", "Etch time (s)"),
                ("fluence", "Ion fluence (ions/cm²)"), ("level", "Level"))


def entry_key(sample_key, level, entries, ei):
    """A region's identity by what it is, not where it sits in a list:
    ``(sample key, level number, spectrum name, region name, occurrence)``
    where ``occurrence`` counts earlier entries of the level that share the
    spectrum and region names (two spectra with one display name). Stable
    across a reload (``Sample.key`` is built from the workbook file id) and
    when other spectra are ticked or unticked, so a user's choice of what
    counts can be saved and handed to the report builders. ``entries`` is
    ``Level.entries``."""
    e = entries[ei]
    spectrum, region = e["spectrum"], e["row"]["region"]
    occ = sum(1 for o in entries[:ei] if o["spectrum"] == spectrum
              and o["row"]["region"] == region)
    return (sample_key, level, spectrum, region, occ)


@dataclass
class Level:
    """One sample at one depth level (``level`` None: not a depth profile)."""
    level: int | None
    etch: float | None = None
    depth: float | None = None
    fluence: float | None = None
    entries: list = field(default_factory=list)   # [{"spectrum", "row"}]
    include: list = field(default_factory=list)   # counts in the total?
    why: list = field(default_factory=list)       # reason it does not
    res: list = field(default_factory=list)       # quant.normalise output

    def group(self):
        """The dict ``quant.profile`` takes."""
        return {"level": self.level, "entries": self.entries}


@dataclass
class Sample:
    key: str
    label: str
    levels: list = field(default_factory=list)
    notes: list = field(default_factory=list)     # what a reader should know
    casaxps: object = None    # casaquant.SampleQuant: CasaXPS's own export,
                              # shown as exported for a sample that has no
                              # fit to tie it to (see collect(); levels is
                              # then empty)
    kind: str = "regions"     # "regions": the ticked high-resolution regions;
                              # "survey": the survey scan, always its own
                              # total; "casaxps": the files as exported
    numbers: str = "fits"     # where the at % come from: "fits" (recomputed
                              # from the fits) or "casaxps" (CasaXPS's own
                              # quantification files, ``casamatch``)
    casa_name: str = ""       # the name the CasaXPS files know the sample by
    dparam: list = field(default_factory=list)   # [{"name","fwhm"}]: the
                              # sample's D-parameter rows from CasaXPS's file
                              # (not a composition: shown beside it)
    rsf_table: list | None = None    # quant.normalise's RSF-fallback args,
    rsf_library: str = "scofield"    # carried so profile_series() can reuse
                                     # them when it re-normalises fresh
    by_hand: tuple = ((), ())        # (left out, counted anyway): the
                                     # regions whose automatic choice the
                                     # user reversed (collect(overrides=))

    @property
    def is_profile(self):
        return len(self.levels) >= 2


@dataclass
class Results:
    samples: list = field(default_factory=list)
    method: str = METHOD

    def __bool__(self):
        return bool(self.samples)

    def children(self):
        """``[(sample key, label)]`` for the Report generator."""
        return [(s.key, s.label) for s in self.samples]

    def chosen(self, skip=()):
        return [s for s in self.samples if s.key not in skip]

    @staticmethod
    def notes_for(samples):
        """The notes of these samples, each said once."""
        return list(dict.fromkeys(n for s in samples for n in s.notes))


def _num(md, key):
    try:
        return float(str((md or {}).get(key, "")).strip())
    except ValueError:
        return None


CASAXPS_NOTE = ("Quantification for this sample is CasaXPS's own exported "
               "result (Quant_survey.txt / Quant_regions.txt / "
               "Quant_Dparam.txt), not recomputed from an embedded fit.")


def collect(docs, display=None, key_of=None, casa_quant=None, ticked=None,
           rsf_table=None, rsf_library="scofield", prefer_csv=False,
           overrides=None, casa_numbers=True):
    """Read the CasaXPS fits of the loaded files (``display`` maps a region to
    the copy that is drawn, with the user's names and binding-energy shift).

    ``rsf_table``/``rsf_library`` are ``quant.normalise``'s own RSF-fallback
    args (a list from ``rsf.load_rsf()``, off by default) -- carried onto
    each ``Sample`` so ``profile_series`` can reuse them later.

    ``overrides`` is the user's own choice of which regions count,
    ``{entry_key: bool}`` (``quantview.ViewState.include``): a region ticked
    here counts although the automatic rules ("counted once", "not the
    preferred line") left it out, an unticked one is left out. It is applied
    after those rules, the automatic notes about a region the user touched
    are not written (they would contradict the table), and one note per
    sample says what was chosen by hand (``Sample.by_hand``). An entry equal
    to the automatic choice, or naming a region that is not here, changes
    nothing. Not the transmission choice: the report applies none.

    ``prefer_csv`` quantifies every fit region that has a complete CasaXPS
    CSV match (``casacsv.py``) from CasaXPS's own exported background and
    curves rather than the reconstruction; a sample note says so.

    A survey scan and the high-resolution regions are **never one total**: a
    sample's survey fits form their own ``Sample`` (``kind == "survey"``,
    labelled "<sample> (survey)"), the ticked high-resolution regions another
    (``kind == "regions"``), each normalised on its own.

    ``casa_quant`` (a ``casaquant.CasaQuant``, see that module) supplies the
    numbers: for a sample it names, each fit region's atomic percent is
    CasaXPS's own (``casamatch``: ``Quant_regions.txt`` for the ticked
    high-resolution regions, ``Quant_survey.txt`` for the survey), renormalised
    over the regions counted; a region the file has no row for is left out and
    says so, never filled with a recomputed value. ``casa_numbers=False`` (the
    Quantification tab's "recompute from fits") ignores the files and
    recomputes every region as before. A sample with no entry in the files is
    recomputed and a note says so. A sample the files name that has no fit at
    all keeps ``Sample.casaxps``, CasaXPS's tables as exported (three
    independent flat lists -- see ``casaquant.py`` for why they are not
    nested). The match strips a "Sample Name: " prefix from the region's own
    sample identifier before comparing: a VAMAS file CasaXPS itself exported
    can carry that same cosmetic prefix in its SAMPLE IDENTIFIER field, while
    the quant text files never repeat it in the samples they list.

    ``ticked`` (a region predicate, e.g. ``lambda r: id(r) in app.checked``)
    restricts the report to what is ticked in the tree: a fit-derived region
    that fails it is skipped, and a ``casa_quant`` sample name is only kept
    when at least one ticked region's sample matches it. ``None`` (the
    default) reads every loaded region, ticked or not, as before."""
    key_of = key_of or reportspec.doc_key
    out = Results()
    by_sample, order = {}, []
    approx, csv_used = set(), set()
    casa_names = set(casa_quant.samples) if casa_quant else set()
    ticked_casa_names = casa_names if ticked is None else set()
    if casa_names and ticked is not None:
        for p in docs:
            for r in p.regions:
                if not ticked(r):
                    continue
                d = display(r) if display else r
                label = (d.sample or r.sample
                        or os.path.basename((p.path or "").rstrip("\\/"))
                        or "sample")
                name = casaquant.strip_sample_prefix(label)
                if name in casa_names:
                    ticked_casa_names.add(name)
    fitted_names = set()          # casa sample names that got a fit-derived total
    for p in docs:
        for r in p.regions:
            if ticked is not None and not ticked(r):
                continue
            if getattr(r, "fit", None) is None or not r.decodable:
                continue
            d = display(r) if display else r
            label = (d.sample or r.sample
                    or os.path.basename((p.path or "").rstrip("\\/"))
                    or "sample")
            rows = quant.fit_rows(d, prefer_csv=prefer_csv)
            if not rows:
                continue
            survey = bool(r.is_survey)
            k = (key_of(p), r.sample, survey)
            if k not in by_sample:
                by_sample[k] = Sample(
                    f"{k[0]}/{k[1]}" + ("#survey" if survey else ""),
                    label + (" (survey)" if survey else ""),
                    rsf_table=rsf_table, rsf_library=rsf_library,
                    kind="survey" if survey else "regions")
                by_sample[k].casa_name = casaquant.strip_sample_prefix(label)
                order.append(k)
            fitted_names.add(by_sample[k].casa_name)
            sample = by_sample[k]
            lv = r.etch_level
            level = next((x for x in sample.levels if x.level == lv), None)
            if level is None:
                md = p.region_metadata(r)
                level = Level(lv, r.etch_time, _num(md, "Depth (nm)"),
                              _num(md, "Fluence (ions/cm²)"))
                sample.levels.append(level)
            if prefer_csv and any(fr.csv_curves is not None
                                  for fr in d.fit.regions):
                csv_used.add(label)
            for row in rows:
                level.entries.append({"spectrum": d.name, "row": row})
                if row["basis"] != "data" or row["approximate"]:
                    approx.add(sample.label)
    for k in order:
        sample = by_sample[k]
        sample.levels.sort(key=lambda x: (x.level is None, x.level or 0))
        _use_casa_numbers(sample, casa_quant if casa_numbers else None)
        settle_sample(sample, overrides, rsf_table, rsf_library,
                      prefer=sample.numbers != "casaxps")
        _element_note(sample)
        _rsf_note(sample)
        _rsf_hint_note(sample)
        if sample.label in csv_used:
            sample.notes.append(
                f"Where CasaXPS's own exported fit curves were imported for "
                f"{sample.label}, the background and areas come from them "
                "rather than from the reconstruction.")
        if sample.label in approx:
            sample.notes.append(
                "The background under some of the fits of "
                f"{sample.label} is not reproduced exactly, so their areas "
                "are approximate.")
        out.samples.append(sample)
    if any(s.numbers == "casaxps" for s in out.samples):
        out.method = METHOD_CASA
    for name in sorted(casa_names):
        if name not in ticked_casa_names or name in fitted_names:
            continue
        sample = Sample(f"casaxps:{name}", name, kind="casaxps",
                        casaxps=casa_quant.samples[name],
                        notes=[CASAXPS_NOTE])
        out.samples.append(sample)
    return out


def _num_cell(v, spec=".2f"):
    return "" if v is None else format(v, spec)


def casaxps_tables(cq):
    """``[(title, header, rows, weights)]`` of CasaXPS's tables as exported
    (survey, regions, D parameter; the ones that have rows), cells already
    formatted, for the reports of a sample that has no fit to tie them to."""
    out = []
    if cq.survey:
        out.append(("Survey (% concentration)", ("Element", "%Conc"),
                    [(r["element"], _num_cell(r["pct"])) for r in cq.survey],
                    [2.0, 1.0]))
    if cq.regions:
        out.append(("Regions (% atomic concentration)",
                    ("Name", "Position (eV)", "%At Conc"),
                    [(r["name"], _num_cell(r["position"], "g"),
                      _num_cell(r["at_pct"])) for r in cq.regions],
                    [2.4, 1.4, 1.2]))
    if cq.dparam:
        out.append(dparam_table(cq.dparam))
    return out


def dparam_table(rows):
    """``(title, header, rows, weights)`` of an Auger D-parameter list."""
    return ("D parameter", ("Name", "FWHM (eV)"),
            [(r["name"], _num_cell(r["fwhm"], "g")) for r in rows],
            [2.0, 1.0])


def _use_casa_numbers(sample, casa_quant):
    """Tag the fit rows of ``sample`` with CasaXPS's own percentages (see
    ``casamatch``) when ``casa_quant`` has a table for it, and say what was
    done. Leaves the sample on recomputed numbers (and says why) when the files
    do not name it or hold no rows for this kind of scan."""
    if casa_quant is None:
        return
    base = sample.label.removesuffix(" (survey)")
    sq = casa_quant.samples.get(sample.casa_name)
    if sq is None:
        sample.notes.append(
            f"CasaXPS's quantification files have no entry for {base}: the "
            "numbers are recomputed from the fits.")
        return
    entries = [e for lv in sample.levels for e in lv.entries]
    if sample.kind == "survey":
        rows, fname = sq.survey, "Quant_survey.txt"
        notes = casamatch.tag_survey(entries, rows)
    else:
        rows, fname = sq.regions, "Quant_regions.txt"
        notes = casamatch.tag_regions(entries, rows, base)
    if not rows:
        sample.notes.append(
            f"CasaXPS's {fname} has no rows for {base}: the numbers are "
            "recomputed from the fits.")
        return
    if sample.kind == "regions":
        sample.dparam = list(sq.dparam)
    sample.numbers = "casaxps"
    sample.notes.append(
        f"Atomic percent for {sample.label} is CasaXPS's own ({fname}), "
        "shared out over the regions counted here; nothing is recomputed.")
    sample.notes.extend(notes)
    lost = sorted({e["row"]["region"] for e in entries
                   if e["row"].get("casa_why")})
    if lost:
        sample.notes.append(
            f"No CasaXPS quantification was found in {fname} for "
            + ", ".join(lost) + f" in {base}: left out of the total.")


def settle_sample(sample, overrides=None, rsf_table=None,
                  rsf_library="scofield", prefer=True):
    """Decide what counts at every level of ``sample`` and normalise it: the
    automatic rules (``_settle``, ``_prefer_lines`` unless ``prefer`` is False:
    with CasaXPS's own numbers the user's tick is the choice), then the user's
    own ticks (``overrides``, see ``collect``). Fills ``Level.include`` /
    ``why`` / ``res``, ``Sample.by_hand`` and the notes that say what was
    done."""
    left_out, counted = [], []
    for level in sample.levels:
        hand = _hand_lookup(overrides, sample.key, level)
        _settle(level, sample.label, sample.notes, hand)
        if prefer:
            _prefer_lines(level, sample.label, sample.notes, rsf_table,
                          rsf_library, hand)
        for ei, now in _apply_overrides(level, hand, rsf_table, rsf_library):
            name = _entry_name(level, ei, len(sample.levels) > 1)
            (counted if now else left_out).append(name)
    sample.notes = list(dict.fromkeys(sample.notes))
    sample.by_hand = (left_out, counted)
    if left_out or counted:
        sample.notes.append(_hand_note(sample.label, left_out, counted))


def _same(a, b):
    """Two names alike apart from case and spaces ("C 1s" and "C1s")."""
    strip = lambda t: "".join(str(t).lower().split())        # noqa: E731
    return strip(a) == strip(b)


_LINE = re.compile(r"^([A-Z][a-z]?)\s*\d+\s*[spdf]")


def element_of(region):
    """The element of a core-level name ("Ti 2p" -> "Ti"), else ''."""
    m = _LINE.match(str(region).strip())
    return m.group(1) if m else ""


def _hand_lookup(overrides, skey, level):
    """``hand(ei)`` -> the user's tick for entry ``ei`` of ``level`` (True /
    False), or None when they made none; None for no overrides at all."""
    if not overrides:
        return None
    return lambda ei: overrides.get(entry_key(skey, level.level,
                                              level.entries, ei))


def _touched(hand, entries, match):
    """True when the user ticked or unticked any entry that ``match``
    selects (an automatic note about those would contradict them)."""
    return hand is not None and any(
        hand(i) is not None for i, e in enumerate(entries) if match(e))


def _apply_overrides(level, hand, rsf_table=None, rsf_library="scofield"):
    """Apply the user's ticks to a settled level and normalise again.
    Returns ``[(entry index, now counted)]`` for the entries whose automatic
    choice was reversed. A reversed exclusion has no reason; a row the user
    unticked says "not included"."""
    if hand is None:
        return []
    changed = []
    for i in range(len(level.entries)):
        want = hand(i)
        if want is None or bool(want) == level.include[i]:
            continue
        level.include[i] = bool(want)
        level.why[i] = "" if want else "not included"
        changed.append((i, bool(want)))
    if changed:
        level.res = quant.normalise(
            [e["row"] for e in level.entries], level.include,
            rsf_table=rsf_table, rsf_library=rsf_library)
        for res, why in zip(level.res, level.why):
            if why:
                res["why"] = why
    return changed


def _entry_name(level, ei, several_levels):
    """"C 1s" (+ " in Survey" when it is a survey's fit, + " at level 3")."""
    e = level.entries[ei]
    name = e["row"]["region"]
    if not _same(e["spectrum"], name):
        name += f" in {e['spectrum']}"
    if several_levels:
        n = level.level if level.level is not None else "-"
        name += f" at level {n}"
    return name


def _hand_note(label, left_out, counted):
    bits = []
    if left_out:
        bits.append("left out: " + ", ".join(left_out))
    if counted:
        bits.append("counted although the automatic rules would not: "
                    + ", ".join(counted))
    return (f"Which regions count in {label} was set by hand ("
            + "; ".join(bits) + ").")


def _settle(level, label, notes, hand=None):
    """Decide what counts at one level, then normalise it. A region name that
    is fitted in more than one spectrum counts once: the dedicated scan (the
    spectrum named after the region) if there is one, else the first."""
    entries = level.entries
    chosen = {}                                  # region -> index counted
    for i, e in enumerate(entries):
        name = e["row"]["region"]
        if name not in chosen or (_same(e["spectrum"], name) and not _same(
                entries[chosen[name]]["spectrum"], name)):
            chosen[name] = i
    level.include, level.why = [], []
    for i, e in enumerate(entries):
        name = e["row"]["region"]
        keep = chosen[name]
        level.include.append(i == keep)
        if i == keep:
            level.why.append("")
        else:
            src = entries[keep]["spectrum"]
            level.why.append("counted once")
            if not _touched(hand, entries,
                            lambda o: o["row"]["region"] == name):
                notes.append(f"{name} is fitted in more than one spectrum "
                             f"of {label}; only the fit in {src} is "
                             "counted.")
    level.res = quant.normalise([e["row"] for e in entries], level.include)
    for res, why in zip(level.res, level.why):
        if why:
            res["why"] = why


# The standard "workhorse" line for an element that has a real, commonly-
# fitted choice between more than one core level -- checked against the
# curated RSF library's own numbers (see rsf.py): the 4f/4d "heavy metal"
# block, where 4f is the accepted standard line even though 4d's RSF can be
# higher (Au: 4f 17.12 vs 4d 19.80, Scofield/Al -- the user's own example).
# Pd/Ag are deliberately NOT here: their libraries show 3d overwhelmingly
# dominant in real use (Ag 3d 18.04 vs 4d 1.55), so there is no genuine
# ambiguity to resolve for them. An element absent from this table is left
# to _element_note's existing "counted from both, with a note" handling.
_PREFERRED_LINE = {"Pt": "4f", "Au": "4f", "Ir": "4f", "Os": "4f", "Re": "4f",
                   "W": "4f", "Ta": "4f", "Hf": "4f", "Hg": "4f", "Tl": "4f",
                   "Pb": "4f", "Bi": "4f"}


def _prefer_lines(level, label, notes, rsf_table=None, rsf_library="scofield",
                  hand=None):
    """After ``_settle``: when more than one of an element's distinct region
    names is still included and one of them is that element's known
    preferred line (``_PREFERRED_LINE``), the others are excluded (``why =
    "not the preferred line"``) and a note is added -- an element absent
    from the table, or whose preferred line was not itself fitted here, is
    left untouched (``_element_note`` still says both are counted). Re-runs
    ``quant.normalise`` the same way ``_settle`` does, this time with
    ``rsf_table``/``rsf_library`` (``quant.normalise``'s own RSF-fallback
    args, off by default) so the final numbers reflect both changes."""
    entries = level.entries
    by_element = {}                   # element -> {region name: entry index}
    for i, e in enumerate(entries):
        if not level.include[i]:
            continue
        region = e["row"]["region"]
        el = element_of(region)
        if el:
            by_element.setdefault(el, {})[region] = i
    for el, regions in by_element.items():
        preferred_orbital = _PREFERRED_LINE.get(el)
        if not preferred_orbital or len(regions) < 2:
            continue
        preferred_name = f"{el} {preferred_orbital}"
        if preferred_name not in regions:
            continue
        losers = sorted(name for name in regions if name != preferred_name)
        if not losers:
            continue
        for name in losers:
            level.include[regions[name]] = False
            level.why[regions[name]] = "not the preferred line"
        if _touched(hand, entries, lambda o: element_of(
                o["row"]["region"]) == el):
            continue
        notes.append(
            f"{el} is fitted from more than one line ({preferred_name}, "
            + ", ".join(losers) + f") in {label}; only {preferred_name}, "
            "the standard line for quantification, is counted.")
    level.res = quant.normalise([e["row"] for e in entries], level.include,
                                rsf_table=rsf_table, rsf_library=rsf_library)
    for res, why in zip(level.res, level.why):
        if why:
            res["why"] = why


def _element_note(sample):
    """One note when an element is counted from two of its lines (its 2p and
    its 1s, say): its share then includes both."""
    lines = {}
    for level in sample.levels:
        for e, inc in zip(level.entries, level.include):
            el = element_of(e["row"]["region"])
            if inc and el and e["row"]["region"] not in lines.setdefault(
                    el, []):
                lines[el].append(e["row"]["region"])
    for el, names in lines.items():
        if len(names) > 1:
            sample.notes.append(
                f"{el} is counted from more than one line ("
                + ", ".join(names) + "), so its share includes both.")


def _rsf_note(sample):
    """One note per region whose atomic percent used a reference-table
    substitute RSF (``quant.normalise``'s ``rsf_source`` -- "component" is
    not noted here, since that is the file's own cross-referenced data, just
    read from a different place, not a substitute). When the substitute is
    one of the IMFP-corrected Scofield tiers (``x["imfp_nm"]``/
    ``x["ke_power_factor"]`` set), the note also says what kinetic-energy
    factor was used, so a reader sees exactly what was assumed, not just a
    combined number."""
    seen = set()
    for level in sample.levels:
        for e, x in zip(level.entries, level.res):
            src = x.get("rsf_source")
            if not src or src == "component":
                continue
            region = e["row"]["region"]
            if region in seen:
                continue
            seen.add(region)
            lib_label = rsf_lib.LIBRARY_SHORT.get(src, src)
            extra = ""
            if x.get("imfp_nm"):
                extra = (f", scaled by a TPP-2M mean free path of "
                        f"{x['imfp_nm']:.3g} nm at that line's own "
                        "kinetic energy")
            elif x.get("ke_power_factor"):
                extra = (", scaled by (kinetic energy)^0.6 -- Thermo "
                        "Avantage's own convention")
            sample.notes.append(
                f"{region}'s sensitivity factor is not recorded in "
                f"{sample.label}; the {lib_label} library's value "
                f"({x['rsf_value']:g}, {x['rsf_anode']} Kα) is used "
                f"instead{extra}.")


def _rsf_hint_note(sample):
    """One note when the RSF reference-table fallback was off (``sample.
    rsf_table`` falsy) and at least one region was left out purely for lack
    of a recorded RSF -- points at the option that might resolve it, so "no
    RSF" in the table doesn't read as a dead end. Never fires when the
    fallback was on (that region simply has no entry in the library either
    -- turning it on again would not help, so no further hint is useful)."""
    if sample.rsf_table:
        return
    if any(x.get("why") == "no RSF"
          for level in sample.levels for x in level.res):
        sample.notes.append(
            f"{sample.label} has a region with no recorded sensitivity "
            "factor, left out of the total. The Report generator's RSF "
            "fallback option (Scofield or Kratos Axis F1s) may be able to "
            "supply one.")


# -- cells ----------------------------------------------------------------------------
def _fmt(v, spec):
    return "" if v is None else format(v, spec)


def _fmt_rms(v):
    """Residual RMS (a fraction of the region's data range) as a percentage,
    or "" when there is no envelope to compare against (a component-less
    survey region, or one whose background is not reproduced)."""
    return "" if v is None else f"{100 * v:.1f}%"


def _fmt_chi2(v):
    """Poisson-weighted reduced chi-square, or "" when it could not be
    computed (no envelope, or dwell/scans -- and so true counts -- unknown)."""
    return "" if v is None else f"{v:.2f}"


def _fmt_rsf(row, x):
    """The RSF cell: the row's own recorded value normally, or, when
    ``quant.normalise``'s reference-table fallback supplied a substitute
    (``rsf_source`` is the library name, not "component"), that value
    labelled with its library and anode so it never reads like a real
    recorded RSF -- e.g. "15.45 (Scofield, Al Kα)"."""
    src = x.get("rsf_source")
    if src and src != "component":
        label = rsf_lib.LIBRARY_SHORT.get(src, src)
        return f"{x['rsf_value']:.4g} ({label}, {x['rsf_anode']} Kα)"
    return _fmt(row.get("rsf"), ".4g")


def composition_header(level):
    """The column titles of ``composition_cells(level)``: the "Area / RSF"
    column holds CasaXPS's own %At when the level's numbers come from its
    quantification files (``casamatch``)."""
    from_casa = any("casa_pct" in e["row"] or "casa_why" in e["row"]
                    for e in level.entries)
    return CASA_HEADER if from_casa else COMPOSITION_HEADER


def composition_cells(level):
    """``[(kind, [cells])]`` for one sample at one level: a "region" row
    (region, background, RSF, area, area / RSF, at %, fit RMS, reduced
    chi-square; the reason instead of the percent when it is left out) with
    a "state" row under it for each chemical state. A region quantified from
    a survey scan is its own sample ("<sample> (survey)"), never a row of the
    high-resolution total. Fit RMS is the residual between the data and the
    fitted envelope as a percentage of the region's data range; reduced
    chi-square is the same residual weighted by Poisson counting statistics
    (data - envelope)^2 / max(data, 1), summed over the region's points and
    divided by points - 1 -- the standard XPS goodness-of-fit figure, and
    comparable across regions in a way the RMS fraction is not. Both are
    blank where there is no envelope (a component-less survey region, or an
    unreproduced background); chi-square is also blank when dwell/scans are
    unknown, since there are then no true counts to weight by. The RSF cell
    shows a substitute's own value and provenance (e.g. "15.45 (Scofield,
    Al Kα)") in place of the row's own (blank/zero) recorded RSF when
    ``quant.normalise``'s reference-table fallback supplied one -- see
    ``_fmt_rsf``; a component-level fallback (the file's own data, just
    summed a different way) leaves the cell as the plain, unlabelled area
    ratio, matching the file's own honest style."""
    rows = []
    for e, x in zip(level.entries, level.res):
        row = e["row"]
        pct = f"{x['at_pct']:.1f}" if x["at_pct"] is not None else x["why"]
        name = row["region"]
        casa = "casa_pct" in row or "casa_why" in row
        rows.append(("region", [name, row.get("background") or "",
                                _fmt_rsf(row, x),
                                _fmt(row.get("area"), ".4g"),
                                (_fmt(row.get("casa_pct"), ".2f") if casa
                                 else _fmt(x["corrected"], ".4g")), pct,
                                _fmt_rms(row.get("rms")),
                                _fmt_chi2(row.get("chi2_red"))]))
        states = (quant.states(row, x["at_pct"])
                  if x["at_pct"] is not None else [])
        if len(states) > 1:                 # one state says nothing more
            for st in states:
                rows.append(("state", ["    " + st["name"], "", "", "",
                                       f"{100 * st['frac']:.0f}% of region",
                                       f"{st['at_pct']:.1f}", "", ""]))
    return rows


def composition_png(level, size=(7.0, 3.0), dpi=200):
    """The composition of one sample at one level as a bar chart (PNG bytes),
    or None without matplotlib: one bar per region with a usable at %, in the
    same order as ``composition_cells``. Chemical states are not broken out
    (the table already does that)."""
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
    except ImportError:
        return None
    import themes
    names, values = [], []
    for e, x in zip(level.entries, level.res):
        if x["at_pct"] is None:
            continue
        row = e["row"]
        names.append(row["region"])
        values.append(x["at_pct"])
    if not names:
        return None
    fig = Figure(figsize=size, dpi=dpi)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0.09, 0.22, 0.88, 0.72))
    cycle = themes.PALETTES["Light"]["cycle"]
    ax.bar(names, values,
                  color=[cycle[i % len(cycle)] for i in range(len(names))])
    ax.set_ylabel("Atomic %", fontsize=9)
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    for i, v in enumerate(values):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=7)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    return buf.getvalue()


def profile_axis(sample):
    """``(label, [x per level])``: depth if every level has a distinct one,
    else etch time, else fluence, else the level number."""
    for attr, label in PROFILE_AXES:
        vals = [getattr(lv, attr) for lv in sample.levels]
        if all(v is not None for v in vals) and len(set(vals)) > 1:
            return label, vals
    return "Level", list(range(1, len(sample.levels) + 1))


def profile_series(sample):
    """``quant.profile`` (atomic percent of each region) for a sample, with
    the same RSF fallback (if any) ``collect()`` used for it."""
    return quant.profile([lv.group() for lv in sample.levels], "element",
                         [lv.include for lv in sample.levels],
                         rsf_table=sample.rsf_table,
                         rsf_library=sample.rsf_library)


def profile_cells(sample):
    """``(header, rows)`` of the depth profile: level, the axis when it is not
    the level, then a column of atomic percent per region."""
    label, xs = profile_axis(sample)
    series = profile_series(sample)["series"]
    header = ["Level"] + ([label] if label != "Level" else []) \
        + [s["name"] for s in series]
    rows = []
    for i, lv in enumerate(sample.levels):
        row = [str(lv.level if lv.level is not None else i + 1)]
        if label != "Level":
            row.append(f"{xs[i]:g}")
        row += [_fmt(s["values"][i], ".1f") for s in series]
        rows.append(row)
    return header, rows


def profile_png(sample, size=(7.0, 3.4), dpi=200, accent=None):
    """The depth profile as a PNG (bytes), or None without matplotlib."""
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
    except ImportError:
        return None
    import themes
    label, xs = profile_axis(sample)
    series = profile_series(sample)["series"]
    if not series:
        return None
    fig = Figure(figsize=size, dpi=dpi)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0.09, 0.17, 0.72, 0.78))
    cycle = themes.PALETTES["Light"]["cycle"]
    for i, s in enumerate(series):
        pts = [(x, v) for x, v in zip(xs, s["values"]) if v is not None]
        if pts:
            ax.plot(*zip(*pts), marker="o", markersize=3.5, linewidth=1.4,
                    color=cycle[i % len(cycle)], label=s["name"])
    ax.set_xlabel(label, fontsize=9)
    ax.set_ylabel("Atomic %", fontsize=9)
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              fontsize=8)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    return buf.getvalue()
