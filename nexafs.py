"""NEXAFS processing (Tk-free, no numpy): ring-current scaling and edge
normalisation.

**Ring-current scaling.** A beamline NEXAFS file records the storage-ring
current at every point (``Region.extra["ring_current_points"]``, set by
``readers/nexus_nexafs.py``). The files do not say what unit it is in, so the
current is not divided by it (that would claim a unit such as A/mA): each
point is scaled to the **mean** ring current of the scan,
``I x mean(ring) / ring``, which keeps the unit of the signal (A) and corrects
the decay of the beam. A point whose ring current is zero, negative or not a
number becomes NaN (it is not guessed). The user chooses this once per
workbook (``Annotations.nexafs_ring``).

**Edge normalisation.** Two methods, chosen per spectrum
(``Annotations.nexafs_edge``; parameters below), applied after the ring-current
scaling:

* ``"step"`` (the standard XANES one): a polynomial (order 0-2, default a
  line) is fitted to a *pre-edge* window and another to a *post-edge* window;
  the pre-edge fit is subtracted and the result divided by the **edge step**,
  ``post(E0) - pre(E0)``, so the step is 1. The sign of the step is kept, so a
  spectrum of negative currents normalises to a positive edge.
* ``"max"`` (the author's own *NEXAFS Processing Template*: a line fitted over
  a pre-edge window is subtracted, then the spectrum is scaled so that its
  minimum is 0 and its maximum 1). The template's last cell multiplies by the
  maximum where this divides by the range; the two differ by the factor
  ``(max - min) x max`` of the background-subtracted curve (0.2 % for the
  spectrum it was checked on) and agree exactly once that factor is applied
  (``tests/test_nexafs_edge.py`` reproduces its cached cells).

Nothing is guessed: ``suggest_edge`` proposes windows (the template's own
automatic rule for the pre-edge: the first 4 eV of the scan) and an edge
energy (the steepest rise), the user confirms them, and a spectrum with no
saved parameters is never normalised. A choice that cannot be applied
(``check_edge``) leaves the spectrum as it was.
"""

from __future__ import annotations

import math

# --- ring-current scaling --------------------------------------------------


def ring_points(region):
    """The per-point ring current of a NEXAFS region (a list as long as its
    data), or None when the file did not record one."""
    ring = (region.extra or {}).get("ring_current_points")
    n = region.n_points
    return ring if ring and n and len(ring) == n else None


def _good(x):
    return isinstance(x, (int, float)) and math.isfinite(x) and x > 0


def mean_ring(ring):
    """Mean of the usable (finite, positive) ring-current values; None when
    there is none."""
    ok = [x for x in ring if _good(x)]
    return sum(ok) / len(ok) if ok else None


def scale_to_mean(counts, ring):
    """``(scaled counts, mean ring current)``, or ``(None, None)`` when the
    ring current has no usable value."""
    mean = mean_ring(ring)
    if mean is None or len(counts) != len(ring):
        return None, None
    nan = float("nan")
    return [c * mean / r if _good(r) else nan
            for c, r in zip(counts, ring)], mean


# --- edge normalisation ------------------------------------------------------
MODES = ("step", "max")
MODE_NAMES = {"step": "Edge step (pre- and post-edge lines, step at E0 = 1)",
              "max": "Pre-edge line, maximum = 1 (the template)"}
ORDER_NAMES = {0: "constant", 1: "linear", 2: "quadratic"}
AUTO_PRE_WIDTH = 4.0        # eV: the template's automatic pre-edge window
AUTO_POST_FRACTION = 0.10   # of the scan: the suggested post-edge window


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _window(w):
    if not isinstance(w, (list, tuple)) or len(w) != 2:
        return None
    lo, hi = _num(w[0]), _num(w[1])
    if lo is None or hi is None:
        return None
    return [min(lo, hi), max(lo, hi)]


def sanitise_edge(d):
    """A valid parameter dict or {}: ``mode`` ("step" / "max"), ``pre``
    ``[lo, hi]`` and ``pre_order`` (0-2, default 1); for "step" also ``post``,
    ``post_order`` and ``e0``. Anything else is dropped."""
    if not isinstance(d, dict) or d.get("mode") not in MODES:
        return {}
    pre = _window(d.get("pre"))
    if pre is None:
        return {}
    out = {"mode": d["mode"], "pre": pre,
           "pre_order": _order(d.get("pre_order"))}
    if d["mode"] == "step":
        post, e0 = _window(d.get("post")), _num(d.get("e0"))
        if post is None or e0 is None:
            return {}
        out.update(post=post, post_order=_order(d.get("post_order")), e0=e0)
    return out


def _order(v):
    v = _num(v)
    return int(v) if v is not None and int(v) in (0, 1, 2) else 1


def fit_poly(xs, ys, order):
    """Least-squares polynomial of ``order`` (0-2) through the points, as
    coefficients lowest power first, or None when there are too few points or
    the system is singular. The abscissa is centred first, so an energy of
    1500 eV does not cost precision."""
    pts = [(x, y) for x, y in zip(xs, ys)
           if math.isfinite(x) and math.isfinite(y)]
    if len(pts) < order + 1:
        return None
    x0 = sum(p[0] for p in pts) / len(pts)
    n = order + 1
    a = [[sum((p[0] - x0) ** (i + j) for p in pts) for j in range(n)]
         for i in range(n)]
    b = [sum(p[1] * (p[0] - x0) ** i for p in pts) for i in range(n)]
    for i in range(n):                          # Gaussian elimination
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        if abs(a[piv][i]) < 1e-300:
            return None
        a[i], a[piv], b[i], b[piv] = a[piv], a[i], b[piv], b[i]
        for r in range(i + 1, n):
            f = a[r][i] / a[i][i]
            for c in range(i, n):
                a[r][c] -= f * a[i][c]
            b[r] -= f * b[i]
    coef = [0.0] * n
    for i in range(n - 1, -1, -1):
        coef[i] = (b[i] - sum(a[i][c] * coef[c]
                              for c in range(i + 1, n))) / a[i][i]
    return {"x0": x0, "c": coef}


def poly_at(fit, x):
    """The fitted polynomial at ``x``."""
    t = x - fit["x0"]
    return sum(c * t ** i for i, c in enumerate(fit["c"]))


def _inside(energy, y, w):
    return [(e, v) for e, v in zip(energy, y)
            if w[0] <= e <= w[1] and math.isfinite(v)]


def check_edge(p, energy, y=None):
    """What stops ``p`` being applied to a scan over ``energy``: a list of
    reasons, empty when it can be. Rules: a valid dict; each window holds more
    points than its order needs; for the step method the pre-edge window ends
    at or below E0, the post-edge window starts at or above it and E0 lies in
    the data; the step is not zero."""
    p = sanitise_edge(p)
    if not p:
        return ["the settings are incomplete"]
    if not energy:
        return ["there are no data"]
    lo, hi = min(energy), max(energy)
    out = []
    yy = y if y is not None else [0.0] * len(energy)
    wins = [("pre-edge", p["pre"], p["pre_order"])]
    if p["mode"] == "step":
        wins.append(("post-edge", p["post"], p["post_order"]))
    for name, w, order in wins:
        if len(_inside(energy, yy, w)) < order + 2:
            out.append(f"the {name} window {w[0]:g}-{w[1]:g} eV holds fewer "
                       f"than {order + 2} points")
        if w[1] < lo or w[0] > hi:
            out.append(f"the {name} window is outside the data "
                       f"({lo:g}-{hi:g} eV)")
    if p["mode"] == "step":
        if not lo <= p["e0"] <= hi:
            out.append(f"E0 {p['e0']:g} eV is outside the data")
        if p["pre"][1] > p["e0"]:
            out.append("the pre-edge window must end at or below E0")
        if p["post"][0] < p["e0"]:
            out.append("the post-edge window must start at or above E0")
    return out


def edge_normalise(energy, y, p):
    """``(normalised y, info)`` for a scan, or ``(None, reasons)`` when ``p``
    cannot be applied. ``info`` has ``mode``, ``pre`` / ``post`` fits, ``e0``
    and for "step" the ``step`` (signed). A point whose value is not finite
    stays NaN."""
    p = sanitise_edge(p)
    problems = check_edge(p, energy, y)
    if problems:
        return None, problems
    pre = fit_poly(*zip(*_inside(energy, y, p["pre"])), p["pre_order"])
    if pre is None:
        return None, ["the pre-edge fit is singular"]
    nan = float("nan")
    z = [v - poly_at(pre, e) if math.isfinite(v) else nan
         for e, v in zip(energy, y)]
    info = {"mode": p["mode"], "pre": pre, "pre_window": p["pre"],
            "pre_order": p["pre_order"]}
    if p["mode"] == "step":
        post = fit_poly(*zip(*_inside(energy, y, p["post"])), p["post_order"])
        if post is None:
            return None, ["the post-edge fit is singular"]
        step = poly_at(post, p["e0"]) - poly_at(pre, p["e0"])
        scale = max((abs(v) for v in y if math.isfinite(v)), default=0.0)
        if not math.isfinite(step) or abs(step) <= 1e-12 * scale or step == 0:
            return None, ["the edge step is zero"]
        info.update(post=post, post_window=p["post"],
                    post_order=p["post_order"], e0=p["e0"], step=step)
        return [v / step for v in z], info
    good = [v for v in z if math.isfinite(v)]
    lo, hi = min(good), max(good)
    size = max((abs(v) for v in y if math.isfinite(v)), default=0.0)
    if hi - lo <= 1e-9 * size:          # only rounding noise is left
        return None, ["the pre-edge-subtracted spectrum is flat"]
    info.update(minimum=lo, maximum=hi)
    return [(v - lo) / (hi - lo) for v in z], info


def suggest_edge(energy, y, mode="step"):
    """Proposed parameters for a scan (a suggestion for the user to confirm,
    never applied by itself): the pre-edge window is the first 4 eV of the
    scan (the template's automatic rule, at most a fifth of the scan); for
    "step" the post-edge window is the last tenth (at least 5 points) and E0
    is the energy of the steepest rise (in the direction the curve ends
    relative to its start) of a 5-point moving average between the two
    windows. Returns {} for fewer than 20 points."""
    pts = [(e, v) for e, v in zip(energy, y)
           if math.isfinite(e) and math.isfinite(v)]
    if len(pts) < 20:
        return {}
    pts.sort()
    es, ys = [p[0] for p in pts], [p[1] for p in pts]
    lo, hi = es[0], es[-1]
    span = hi - lo
    pre = [lo, lo + min(AUTO_PRE_WIDTH, span / 5.0)]
    if mode == "max":
        return {"mode": "max", "pre": pre, "pre_order": 1}
    post_lo = min(hi - span * AUTO_POST_FRACTION, es[-5])
    sm = [sum(ys[max(0, i - 2):i + 3]) / len(ys[max(0, i - 2):i + 3])
          for i in range(len(ys))]
    direction = 1.0 if sm[-1] >= sm[0] else -1.0
    best, e0 = None, (pre[1] + post_lo) / 2.0
    for i in range(1, len(es) - 1):
        if not pre[1] <= es[i] <= post_lo:
            continue
        slope = direction * (sm[i + 1] - sm[i - 1]) / (es[i + 1] - es[i - 1])
        if best is None or slope > best:
            best, e0 = slope, es[i]
    return {"mode": "step", "pre": pre, "pre_order": 1,
            "post": [post_lo, hi], "post_order": 1, "e0": e0}


def edge_text(p):
    """The windows as the metadata and the methods text state them."""
    p = sanitise_edge(p)
    if not p:
        return ""
    pre = (f"pre-edge {p['pre'][0]:g}-{p['pre'][1]:g} eV "
           f"({ORDER_NAMES[p['pre_order']]})")
    if p["mode"] == "step":
        return (f"edge step: {pre}, post-edge {p['post'][0]:g}-"
                f"{p['post'][1]:g} eV ({ORDER_NAMES[p['post_order']]}), "
                f"E0 {p['e0']:g} eV; step set to 1")
    return (f"{pre} subtracted, then scaled to a minimum of 0 and a "
            "maximum of 1")


# --- both together, as a spectrum is drawn and described ---------------------
class Result:
    """What ``process`` made of a NEXAFS spectrum: ``counts`` (a list), the
    ``mean`` ring current it was scaled to (None when it was not), ``edge``
    (the ``edge_normalise`` info, None when no edge normalisation applies),
    and ``problems`` (why a chosen edge normalisation could not be applied)."""

    def __init__(self, counts, mean=None, edge=None, problems=()):
        self.counts, self.mean, self.edge = counts, mean, edge
        self.problems = list(problems)


def process(region, ring_on, edge_params):
    """The ring-current scaling (``ring_on``) and then the edge normalisation
    (``edge_params``, saved by the user) of a NEXAFS ``region``; None when
    neither applies."""
    counts, mean = region.counts, None
    if ring_on:
        ring = ring_points(region)
        if ring is not None:
            scaled, m = scale_to_mean(region.counts, ring)
            if scaled is not None:
                counts, mean = scaled, m
    edge, problems = None, []
    if edge_params:
        out, info = edge_normalise(region.energy, counts, edge_params)
        if out is not None:
            counts, edge = out, info
        else:
            problems = info
    if mean is None and edge is None:
        return None
    return Result(counts, mean, edge, problems)
