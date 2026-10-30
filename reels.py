"""REELS (reflection electron energy loss spectroscopy): the band gap from
the onset of inelastic losses.

A REELS spectrum is intensity against kinetic energy, with the elastic peak
at the primary energy; ``loss = E_elastic - KE`` is how much energy each
electron gave up. The band gap is read from the *onset*: a straight line (the
tangent) is drawn through two points the user picks on the rising edge of the
loss spectrum, and the gap is where that line meets the baseline level of the
flat region between the elastic peak and the onset.

Pure functions, no Tk.
"""

from __future__ import annotations

import math

ELASTIC_TAIL = 0.5          # eV: the elastic peak's own tail is not baseline


def elastic_peak(ke, counts):
    """Kinetic energy of the strongest point (the elastic peak), or None."""
    if not ke or not counts or len(ke) != len(counts):
        return None
    i = max(range(len(counts)), key=counts.__getitem__)
    return ke[i]


def loss_axis(ke, elastic):
    """Energy loss for each kinetic energy."""
    return [elastic - k for k in ke]


def _median(vals):
    vals = sorted(vals)
    n = len(vals)
    if not n:
        return None
    mid = n // 2
    return vals[mid] if n % 2 else (vals[mid - 1] + vals[mid]) / 2.0


def baseline(loss, y, onset):
    """The flat level before the onset: the median intensity between the end
    of the elastic tail and the lower of the picked points; with too few
    points there, the lowest intensity up to the higher pick."""
    lo_x = min(onset)
    sel = [v for x, v in zip(loss, y) if ELASTIC_TAIL < x < lo_x]
    if len(sel) > 3:
        return _median(sel)
    hi_x = max(onset)
    sel = [v for x, v in zip(loss, y) if 0 < x < hi_x]
    return min(sel) if sel else None


def band_gap(loss, y, p1, p2, base=None):
    """``{"gap", "slope", "base"}`` from the tangent through ``p1`` and
    ``p2`` (each ``(loss, intensity)``), or None when the two points do not
    define a rising edge that meets the baseline at a positive loss."""
    (x1, y1), (x2, y2) = p1, p2
    if not all(math.isfinite(v) for v in (x1, y1, x2, y2)) or x1 == x2:
        return None
    slope = (y2 - y1) / (x2 - x1)
    if slope <= 0:
        return None
    if base is None:
        base = baseline(loss, y, (x1, x2))
    if base is None:
        return None
    gap = x1 + (base - y1) / slope
    if not math.isfinite(gap) or gap <= 0:
        return None
    return {"gap": gap, "slope": slope, "base": base}


BACKGROUNDS = ("Linear", "Shirley", "None")
METHODS = ("gauss", "d2")
FEATURES = {"π+σ": (15.0, 45.0), "π": (3.0, 12.0)}   # default loss windows (eV)
GAUSS_FRAC = 0.9            # Filippi & Calliari 2006: the data above 90 % of the peak
FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))
PLASMA_EV = 28.816          # hbar*omega_p (eV) = PLASMA_EV * sqrt(Z rho / A) for free electrons
BG_POINTS = 80              # points kept of the background curve for drawing


def _sorted(loss, y):
    """``(loss, y)`` as float arrays, ascending loss, non-finite points dropped."""
    import numpy as np
    x = np.asarray(loss, dtype=float)
    v = np.asarray(y, dtype=float)
    if x.shape != v.shape or x.ndim != 1:
        return None, None
    ok = np.isfinite(x) & np.isfinite(v)
    x, v = x[ok], v[ok]
    order = np.argsort(x, kind="stable")
    return x[order], v[order]


def window_background(xs, ys, kind):
    """The background under ``ys`` (ascending loss) between the window's ends.

    ``Linear`` joins the two end levels; ``Shirley`` is
    ``lineshapes.shirley`` run the way XPS runs it, i.e. on ascending kinetic
    energy, which is descending loss (the high-loss end is the higher one).
    A Shirley background whose high-loss end lies *below* its low-loss end is
    refused (None): the step would be the wrong way round. ``None`` is zero."""
    import numpy as np
    import lineshapes
    if kind == "None":
        return np.zeros(len(ys))
    if kind == "Linear":
        return lineshapes.linear_bg(ys, avg=1)
    if kind == "Shirley":
        if ys[-1] < ys[0]:
            return None
        return lineshapes.shirley(ys[::-1], avg=1)[::-1]
    return None


def _gaussian(x, y):
    """``(centre, sigma, amplitude)`` of the least-squares Gaussian through
    ``(x, y)``: a log-parabola start refined by Gauss-Newton (numpy only)."""
    import numpy as np
    pos = y > 0
    if pos.sum() < 3:
        return None
    c2, c1, c0 = np.polyfit(x[pos], np.log(y[pos]), 2)
    if c2 >= 0:
        return None
    centre = -c1 / (2.0 * c2)
    sigma = math.sqrt(-1.0 / (2.0 * c2))
    amp = math.exp(c0 - c1 * c1 / (4.0 * c2))
    p = np.array([amp, centre, sigma])
    for _ in range(30):
        a, c, s = p
        g = a * np.exp(-0.5 * ((x - c) / s) ** 2)
        r = y - g
        J = np.column_stack([g / a, g * (x - c) / s ** 2,
                             g * (x - c) ** 2 / s ** 3])
        try:
            step, *_ = np.linalg.lstsq(J, r, rcond=None)
        except np.linalg.LinAlgError:
            break
        p = p + step
        if p[2] <= 0:
            return None
        if float(np.abs(step).max()) < 1e-10 * max(1.0, float(np.abs(p).max())):
            break
    a, c, s = p
    if not all(math.isfinite(v) for v in (a, c, s)) or s <= 0 or a <= 0:
        return None
    return float(c), float(abs(s)), float(a)


def _gauss_peak(xs, ys, kind, frac):
    bg = window_background(xs, ys, kind)
    if bg is None:
        return None
    d = ys - bg
    top = int(d.argmax())
    if d[top] <= 0:
        return None
    keep = d >= frac * d[top]
    lo = hi = top                       # the run of points round the maximum
    while lo > 0 and keep[lo - 1]:
        lo -= 1
    while hi < len(d) - 1 and keep[hi + 1]:
        hi += 1
    if hi - lo + 1 < 3:
        return None
    fit = _gaussian(xs[lo:hi + 1], d[lo:hi + 1])
    if fit is None:
        return None
    centre, sigma, _amp = fit
    return centre, sigma * FWHM_PER_SIGMA, bg


def _second_derivative_peak(x, y, window, smooth):
    """Position of the maximum of -d2y/dx2 of the Gaussian-smoothed spectrum
    inside ``window`` (the smoothing width ``smooth`` is a sigma in eV)."""
    import numpy as np
    lo, hi = window
    if len(x) < 5 or smooth <= 0:
        return None
    h = float(np.median(np.diff(x)))
    if not h > 0:
        return None
    grid = np.arange(x[0], x[-1] + h / 2, h)
    v = np.interp(grid, x, y)
    half = max(1, int(math.ceil(4.0 * smooth / h)))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) * h / smooth) ** 2)
    k /= k.sum()
    pad = np.pad(v, half, mode="edge")
    sm = np.convolve(pad, k, mode="valid")
    d2 = -np.gradient(np.gradient(sm, h), h)
    sel = (grid >= lo) & (grid <= hi)
    if sel.sum() < 3:
        return None
    idx = np.flatnonzero(sel)
    i = int(idx[d2[idx].argmax()])
    if d2[i] <= 0:
        return None
    if 0 < i < len(grid) - 1:                 # parabolic refinement
        a, b, c = d2[i - 1], d2[i], d2[i + 1]
        den = a - 2 * b + c
        if den < 0:
            return float(grid[i] + 0.5 * h * (a - c) / den)
    return float(grid[i])


def plasmon(loss, y, window, background="Linear", method="gauss",
            frac=GAUSS_FRAC, smooth=1.0):
    """The energy of a loss peak (a plasmon), or None when it is not defined.

    ``window`` is ``(low loss, high loss)`` in eV. ``method`` ``"gauss"``
    follows Filippi & Calliari (*Surf. Interface Anal.* 38, 595, 2006): a
    background (``background``) is taken from the window's two ends and the
    centre of the least-squares Gaussian through the data above ``frac`` of
    the peak height is Ep (with the Gaussian's FWHM). ``"d2"`` takes the
    maximum of the negated second derivative of the spectrum smoothed with a
    Gaussian of sigma ``smooth`` eV: it needs no background and is reported
    with background "None" (they find 26.9 eV for graphite by it, against
    27.1 / 26.2 eV with a linear / Shirley background).

    Returns ``{ep, fwhm, method, background, window, frac, smooth, bg}``;
    ``fwhm`` is None for ``d2`` and ``bg`` is ``[[loss, level], ...]`` for
    drawing (empty without a background)."""
    import numpy as np
    try:
        lo, hi = sorted((float(window[0]), float(window[1])))
    except (TypeError, ValueError, IndexError):
        return None
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        return None
    if method not in METHODS or background not in BACKGROUNDS:
        return None
    x, v = _sorted(loss, y)
    if x is None or len(x) < 5:
        return None
    sel = (x >= lo) & (x <= hi)
    if sel.sum() < 5:
        return None
    if method == "gauss":
        if not 0.0 < frac < 1.0:
            return None
        got = _gauss_peak(x[sel], v[sel], background, frac)
        if got is None:
            return None
        ep, fwhm, bg = got
        step = max(1, len(bg) // BG_POINTS)
        curve = [[float(a), float(b)] for a, b in
                 zip(x[sel][::step], bg[::step])]
        used_bg = background
    else:
        ep = _second_derivative_peak(x, v, (lo, hi), float(smooth))
        if ep is None:
            return None
        fwhm, curve, used_bg = None, [], "None"
    if not math.isfinite(ep):
        return None
    return {"ep": float(ep), "fwhm": fwhm, "method": method,
            "background": used_bg, "window": [lo, hi],
            "frac": float(frac) if method == "gauss" else None,
            "smooth": float(smooth) if method == "d2" else None,
            "bg": curve}


def method_text(res):
    """One sentence saying how a plasmon energy was found."""
    lo, hi = res["window"]
    if res["method"] == "gauss":
        return (f"Gaussian fitted above {res['frac'] * 100:g} % of the peak "
                f"height, {res['background'].lower()} background, "
                f"window {lo:g}–{hi:g} eV")
    return (f"maximum of the negated second derivative (smoothing "
            f"{res['smooth']:g} eV), window {lo:g}–{hi:g} eV")


def density_from_plasmon(ep, z_valence=4, a=12.011, m_eff=0.87):
    """Mass density (g/cm3) from a bulk plasmon energy (eV):
    ``hbar*omega_p = sqrt(n e^2 / (eps0 m*))`` with ``n`` the valence electron
    density, ``m* = m_eff`` electron masses. The defaults are the amorphous
    carbon of Yang et al. (*Phys. Chem. Chem. Phys.* 23, 25335, 2021: 25.5 eV
    gives 2.046 g/cm3). A raw REELS peak sits higher than the bulk-loss peak
    (surface excitation, multiple scattering), so it overstates the density."""
    if not (ep and ep > 0 and z_valence > 0 and a > 0 and m_eff > 0):
        return None
    return ep * ep * a * m_eff / (z_valence * PLASMA_EV ** 2)


def tangent_points(result, p1, p2, span=1.25):
    """Two ``(loss, intensity)`` ends for drawing: from the baseline crossing
    to a little beyond the higher pick."""
    gap = result["gap"]
    x_end = max(p1[0], p2[0])
    x_end = gap + (x_end - gap) * span
    y_end = result["base"] + result["slope"] * (x_end - gap)
    return (gap, result["base"]), (x_end, y_end)
