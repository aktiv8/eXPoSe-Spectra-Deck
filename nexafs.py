"""NEXAFS ring-current scaling (Tk-free, no numpy).

A beamline NEXAFS file records the storage-ring current at every point
(``Region.extra["ring_current_points"]``, set by ``readers/nexus_nexafs.py``).
The files do not say what unit it is in, so the current is not divided by it
(that would claim a unit such as A/mA): each point is scaled to the **mean**
ring current of the scan, ``I x mean(ring) / ring``, which keeps the unit of
the signal (A) and corrects the decay of the beam. A point whose ring current
is zero, negative or not a number becomes NaN (it is not guessed). The user
chooses this once per workbook (``Annotations.nexafs_ring``); it is applied on
the way out by ``Workspace._display`` and never to what a reader returned.
"""

from __future__ import annotations

import math


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
