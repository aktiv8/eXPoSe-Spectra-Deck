"""What kind of measurement a spectrum belongs to, as the instrument file says
(Tk-free; the readers fill ``Region.extra["block"]`` and
``SpectrumFile.region_metadata`` shows it as "Block type").

A block is a label, never an inference about purpose: it comes from what the
file records and the reason is given with it.

* **Kratos ESCApe** writes each run with a class name in front of its spectra
  (``Kernel.Core.Instrument.Acquisition.<Name>Data``); ``RUN_CLASSES`` names the
  five seen in 13 files. Every run of the SR and EA sessions of D:\\Temp\\hxps
  has exactly one such marker and the number of spectra after it equals
  XPSView's own run (6 + 8 and 56 runs; the 112-spectrum depth profiles, the
  7-spectrum sequences and the 35-spectrum angle-resolved runs). XPSView calls
  the angle-resolved runs a possible "X-ray exposure test"; the file says
  angle-resolved.
* **Thermo Avantage** stores a repeated scan, a depth profile and a tilt series
  alike as levels of one file; which it is shows in the level axis: an etch time
  is a depth profile, a stage tilt an angle-resolved series, and levels with
  neither are the same scan repeated in succession (Avantage's "Iteration").
"""

from __future__ import annotations

import re

LABELS = {
    "sequence": "Acquisition sequence",
    "depth_profile": "Depth profile",
    "angle_resolved": "Angle-resolved series",
    "etch": "Etch",
    "repeat_sequence": "Repeated sequence",
    "iteration_series": "Iteration series (the same scan repeated)",
    "multi_point": "Multi-point scan",
}

# Kratos ESCApe run classes -> (kind, what the instrument recorded)
RUN_CLASSES = {
    "UserSequenceData": ("sequence", "recorded by the instrument as a user "
                                     "sequence"),
    "DepthProfileData": ("depth_profile", "recorded by the instrument as a "
                                          "depth-profile run"),
    "AngleResolvedData": ("angle_resolved", "recorded by the instrument as an "
                                            "angle-resolved run"),
    "EtchData": ("etch", "recorded by the instrument as an etch run"),
    "RepeatSequenceData": ("repeat_sequence", "recorded by the instrument as "
                                              "a repeated sequence"),
}
RUN_MARKER = re.compile(
    rb"Kernel\.Core\.Instrument\.Acquisition\.(\w+Data)\b")


def make(kind, reason, run=None, size=None):
    """The ``Region.extra["block"]`` dict: ``kind`` (a ``LABELS`` key), the
    ``reason`` it was given, the run's number (from 1) and how many spectra it
    holds, when known."""
    b = {"kind": kind, "reason": reason}
    if run is not None:
        b["run"] = run
    if size is not None:
        b["size"] = size
    return b


def text(b):
    """"Depth profile (recorded by the instrument as a depth-profile run; run
    3, 39 spectra)" for the metadata row, or "" for no block."""
    if not isinstance(b, dict) or b.get("kind") not in LABELS:
        return ""
    bits = [b.get("reason", "")]
    if b.get("run") is not None:
        bits.append(f"run {b['run']}" + (f", {b['size']} spectra"
                                         if b.get("size") else ""))
    inner = "; ".join(x for x in bits if x)
    return LABELS[b["kind"]] + (f" ({inner})" if inner else "")


def kratos_runs(raw, regions):
    """Give each region of a Kratos ESCApe ``.experiment`` the block of the run it
    sits in: the last run marker in the file before its spectrum. A region
    before any marker, or a file with none, gets nothing. Returns the runs that
    have spectra, in file order: ``[{"run", "kind", "start", "end",
    "regions"}]`` with ``start`` / ``end`` the byte range from its marker to the
    next one (the run's own records, such as its etch durations, lie there)."""
    marks = [(m.start(), m.group(1).decode()) for m in RUN_MARKER.finditer(raw)
             if m.group(1).decode() in RUN_CLASSES]
    if not marks:
        return []
    starts = [o for o, _n in marks]
    import bisect
    members = {}
    for r in regions:
        i = bisect.bisect_right(starts, r.offset) - 1
        if i >= 0:
            members.setdefault(i, []).append(r)
    runs = []
    for i in sorted(members):
        rs = members[i]
        kind, why = RUN_CLASSES[marks[i][1]]
        for r in rs:
            r.extra["block"] = make(kind, why, run=i + 1, size=len(rs))
        runs.append({"run": i + 1, "kind": kind, "regions": rs,
                     "start": starts[i],
                     "end": starts[i + 1] if i + 1 < len(starts) else len(raw)})
    return runs


def levels_block(n_levels, etch=False, tilt=False):
    """The block of an Avantage file whose scan is repeated as levels."""
    if etch:
        return make("depth_profile", "the levels have etch times",
                    size=n_levels)
    if tilt:
        return make("angle_resolved", "the levels are stage tilts",
                    size=n_levels)
    return make("iteration_series",
                f"{n_levels} levels with no etch time and no stage tilt",
                size=n_levels)
