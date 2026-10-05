"""A saved figure's zoom (Tk-free).

A figure stores the look (``Workspace.capture_state``) and, since a zoom is not
part of that look, the energy / intensity limits of each panel that was zoomed
when it was saved, as ``state["zoom"]``::

    {"scale": "Binding", "panels": {panel key: {"x": [a, b], "y": [c, d]}}}

A figure saved with nothing zoomed stores an empty ``panels`` (the full range,
and a Recall clears an older zoom); a figure without a ``zoom`` key predates
this and leaves the live zoom alone. Limits are in the drawn frame, so they are
only used under the energy scale they were taken in.
"""

from __future__ import annotations

import math

REL_TOL = 1e-6


def _pair(v):
    if (isinstance(v, (list, tuple)) and len(v) == 2
            and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                    and math.isfinite(x) for x in v)):
        return [float(v[0]), float(v[1])]
    return None


def encode(panels, scale):
    """``panels``: ``{key: (xlim, ylim)}`` of the zoomed panels -> the value
    stored in a figure's state."""
    out = {}
    for key, (x, y) in panels.items():
        xs, ys = _pair(x), _pair(y)
        if xs is not None and ys is not None:
            out[str(key)] = {"x": xs, "y": ys}
    return {"scale": str(scale), "panels": out}


def sanitise(zoom, scale):
    """The saved zoom as ``{key: {"x": [a, b], "y": [c, d]}}``; ``None`` when
    the state has none (an older figure, or malformed data), ``{}`` for a
    figure saved unzoomed *or* taken under another energy scale (its limits
    would mean something else)."""
    if not isinstance(zoom, dict) or not isinstance(zoom.get("panels"), dict):
        return None
    if zoom.get("scale") != scale:
        return {}
    out = {}
    for key, v in zoom["panels"].items():
        if not isinstance(v, dict):
            continue
        xs, ys = _pair(v.get("x")), _pair(v.get("y"))
        if xs is not None and ys is not None and xs[0] != xs[1] \
                and ys[0] != ys[1]:
            out[str(key)] = {"x": xs, "y": ys}
    return out


def _same(a, b):
    span = max(abs(a[0] - a[1]), abs(b[0] - b[1]), 1e-300)
    return all(abs(p - q) <= REL_TOL * span for p, q in zip(a, b))


def is_zoomed(auto, current):
    """True when a panel's current ``(xlim, ylim)`` differ from what it was
    drawn with (``auto``)."""
    return not (_same(auto[0], current[0]) and _same(auto[1], current[1]))


def fits(auto_x, saved_x):
    """A saved x range is usable when it overlaps the panel's own (auto)
    range: if the data have changed since the figure was saved so that they
    no longer do, the panel keeps its full range instead of showing nothing."""
    lo, hi = sorted(auto_x)
    a, b = sorted(saved_x)
    return a < hi and b > lo
