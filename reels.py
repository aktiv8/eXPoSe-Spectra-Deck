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


# -- the pi+sigma loss region as Gaussians, and a layer count from it ------------------------
# Stobinski et al., J. Electron Spectrosc. Relat. Phenom. 195, 145 (2014): the
# loss region of graphite-like carbon is fitted with Gaussians at the loss
# energies they list (their text and Fig. 5); the pi peak has a width of its
# own, the other components share one. The positions are an assignment taken
# from the literature, so they are only starting values, free within ``tol``.
LOSS_COMPONENTS = (            # key, name, loss (eV), width group
    ("pi", "π", 6.0, "pi"),
    ("inter", "interband", 12.5, "sigma"),
    ("sp2S", "C sp2 surface", 19.5, "sigma"),
    ("sp3S", "C sp3 surface", 23.0, "sigma"),
    ("sp2B", "C sp2 bulk", 27.0, "sigma"),
    ("sp3B", "C sp3 bulk", 33.8, "sigma"),
    ("ss35", "σ→σ* 35", 35.0, "sigma"),
    ("ss42", "σ→σ* 42", 42.0, "sigma"),
)
LOSS_WINDOW = (3.0, 50.0)      # eV: below that is the elastic peak's tail
LOSS_TOL = 1.0                 # eV a position may move from its starting value
LOSS_WIDTHS = {"pi": (2.5, 0.8, 6.0), "sigma": (7.0, 2.0, 16.0)}  # start, low, high
GRAPHITE_R = 4.2               # sp2 bulk / sp2 surface of their graphite
LAYER_NM = 0.335               # graphite interlayer distance (nm)
_UNIT_AREA = 2.0 * math.sqrt(math.log(2.0) / math.pi)


def gaussian_unit(x, centre, fwhm):
    """A Gaussian of unit area (per eV) with the given full width at half
    maximum, on the array ``x``."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    return (_UNIT_AREA / fwhm) * np.exp(-4.0 * math.log(2.0)
                                        * ((x - centre) / fwhm) ** 2)


def _golden(f, lo, hi, x0, evals=14):
    """The minimiser of ``f`` on ``[lo, hi]`` by golden section, or ``x0`` when
    that is no worse (so a sweep never makes the fit worse)."""
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(evals):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    best = c if fc < fd else d
    return best if f(best) < f(x0) else x0


def decompose(loss, y, window=LOSS_WINDOW, background="Linear", tol=LOSS_TOL,
              components=None):
    """Fit the loss region with the Gaussians of :data:`LOSS_COMPONENTS`.

    ``components`` is ``[(key, name, loss, group)]`` (default the table above);
    ``tol`` is how far each position may move from its starting value (0 holds
    it). The widths are refined too: one for the pi peak, one shared by the
    rest. Amplitudes (the areas) and nothing else are solved by non-negative
    least squares, inside a bounded one-parameter search of each position and
    width, so the fit has no starting-amplitude dependence. ``background``
    is taken from the window's two ends as for :func:`plasmon`.

    Returns ``{window, background, tol, components: [{key, name, loss, fwhm,
    area}], bg, rms, r}`` or None when the window holds too few points, the
    background cannot be taken, or nothing fits. ``rms`` is the residual's
    root mean square as a fraction of the data's range in the window and ``r``
    is sp2 bulk / sp2 surface (None when there is no sp2 surface area)."""
    import numpy as np
    import snapfit
    comps = [tuple(c) for c in
             (LOSS_COMPONENTS if components is None else components)]
    try:
        lo, hi = sorted((float(window[0]), float(window[1])))
        tol = max(0.0, float(tol))
    except (TypeError, ValueError, IndexError):
        return None
    if not (math.isfinite(lo) and math.isfinite(hi) and math.isfinite(tol)) \
            or hi <= lo or background not in BACKGROUNDS or not comps:
        return None
    x, v = _sorted(loss, y)
    if x is None:
        return None
    sel = (x >= lo) & (x <= hi)
    xs, ys = x[sel], v[sel]
    if len(xs) < max(20, 3 * len(comps)):
        return None
    bg = window_background(xs, ys, background)
    if bg is None:
        return None
    d = ys - bg
    centres = [float(c[2]) for c in comps]
    groups = sorted({c[3] for c in comps})
    widths = {g: LOSS_WIDTHS.get(g, LOSS_WIDTHS["sigma"])[0] for g in groups}

    def solve(cs, ws):
        A = np.column_stack([gaussian_unit(xs, c, ws[comp[3]])
                             for c, comp in zip(cs, comps)])
        amp = snapfit.nnls(A, d)
        r = d - A @ amp
        return float(r @ r), amp

    def cost_of(cs, ws):
        return solve(cs, ws)[0]

    for _ in range(3):                                   # sweeps
        for i, c0 in enumerate(centres):
            if tol <= 0:
                break
            centres[i] = _golden(
                lambda c: cost_of(centres[:i] + [c] + centres[i + 1:], widths),
                comps[i][2] - tol, comps[i][2] + tol, c0)
        for g in groups:
            _s, wl, wh = LOSS_WIDTHS.get(g, LOSS_WIDTHS["sigma"])
            widths[g] = _golden(
                lambda w: cost_of(centres, {**widths, g: w}), wl, wh,
                widths[g])
    sse, amp = solve(centres, widths)
    if not np.isfinite(amp).all() or not (amp > 0).any():
        return None
    span = float(d.max() - d.min()) or 1.0
    step = max(1, len(xs) // BG_POINTS)
    out = [{"key": c[0], "name": c[1], "loss": float(cs), "area": float(a),
            "fwhm": float(widths[c[3]])}
           for c, cs, a in zip(comps, centres, amp)]
    area = {c["key"]: c["area"] for c in out}
    s2 = area.get("sp2S", 0.0)
    return {"window": [lo, hi], "background": background, "tol": tol,
            "components": out,
            "bg": [[float(a), float(b)] for a, b in zip(xs[::step],
                                                        bg[::step])],
            "rms": math.sqrt(sse / len(xs)) / span,
            "r": (area.get("sp2B", 0.0) / s2
                  if s2 > 1e-9 * (sum(area.values()) or 1.0) else None)}


def loss_curves(fit, x):
    """The fitted components and their sum on the loss array ``x``:
    ``(sum including the background, [(component, curve)])``."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    bg = np.asarray(fit.get("bg") or [], dtype=float)
    base = (np.interp(x, bg[:, 0], bg[:, 1]) if len(bg) > 1
            else np.zeros(len(x)))
    parts = [(c, c["area"] * gaussian_unit(x, c["loss"], c["fwhm"]))
             for c in fit["components"]]
    return base + sum((p for _c, p in parts), np.zeros(len(x))), parts


# -- the pi / (sigma+pi) intensity ratio and an sp2 fraction against HOPG ----------------
# The ratio of the areas of the pi and the (sigma+pi) losses, each above a
# straight background between its own window's limits, as HarwellXPS's REELS
# analysis defines it (windows 3.5-8.0 and 10.0-35.0 eV of loss, the line
# anchored on the mean of 3 channels centred on each limit, trapezoidal areas,
# the loss zero at each spectrum's own parabola-refined elastic peak; checked
# digit for digit on the six REELS spectra of D:\Temp\hxps\BM-20263874-25-55
# against the rows its .xpsv stores). It rises with sp2 content but
# depends on the windows, the beam energy and the geometry, so on its own it
# only compares spectra measured together. Against a measured HOPG (sp2, set to
# 1) and, optionally, diamond (no pi, set to 0) it becomes a fraction on a
# scale that removes those dependences *if the references were measured the
# same way*; nothing here is a published reference value.
PI_WINDOW = (3.5, 8.0)
SP_WINDOW = (10.0, 35.0)
EDGE_AVG = 3
SCALE_LIMITS = (0.0, 1.15)     # a fraction outside this is flagged, never clipped
ELASTIC_TOL = 0.01             # relative difference in elastic-peak energy that is noted
SENSITIVITY_SHIFT = 0.5        # eV a window limit is moved by window_sensitivity


def elastic_peak_refined(ke, counts):
    """Kinetic energy of the elastic peak: the strongest point, refined by a
    parabola through it and its two neighbours; the plain maximum at an end
    of the spectrum or on a flat top. None for no data."""
    if not ke or not counts or len(ke) != len(counts):
        return None
    i = max(range(len(counts)), key=counts.__getitem__)
    e = ke[i]
    if 0 < i < len(ke) - 1:
        a, b, c = counts[i - 1], counts[i], counts[i + 1]
        den = a - 2 * b + c
        if den < 0:
            step = (ke[i + 1] - ke[i - 1]) / 2.0
            return float(ke[i] + 0.5 * step * (a - c) / den)
    return float(e)


def _window_area(x, v, window, avg):
    """``(area, loss, background)`` of ``v`` above the straight line between the
    window's first and last points, by the trapezoid rule, or None when the
    window holds too few points. Each end of the line is the mean of ``avg``
    channels **centred on** that end point, so the channel just outside the
    window counts (for ``avg`` 3: the point and one either side): that is what
    reproduces HarwellXPS's areas to every digit (the six REELS spectra of
    ``D:\\Temp\\hxps\\BM-20263874-25-55``); the first ``avg`` channels inside
    the window do not."""
    import numpy as np
    try:
        lo, hi = sorted((float(window[0]), float(window[1])))
    except (TypeError, ValueError, IndexError):
        return None
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        return None
    idx = np.flatnonzero((x >= lo) & (x <= hi))
    if len(idx) < max(5, avg):
        return None
    i0, i1 = int(idx[0]), int(idx[-1])

    def level(i):
        a = max(0, i - avg // 2)
        return float(v[a:a + avg].mean())
    xs, ys = x[i0:i1 + 1], v[i0:i1 + 1]
    lo_y, hi_y = level(i0), level(i1)
    bg = lo_y + (hi_y - lo_y) * (xs - xs[0]) / (xs[-1] - xs[0])
    d = ys - bg
    return float(np.sum((d[1:] + d[:-1]) * np.diff(xs)) / 2.0), xs, bg


def pi_ratio(loss, y, pi_window=PI_WINDOW, sp_window=SP_WINDOW,
             avg=EDGE_AVG):
    """The π / (σ+π) ratio of a REELS spectrum, or None when a window holds too
    few points.

    ``loss`` is energy loss (eV) from the elastic peak, ``y`` the intensity.
    Returns ``{pi_window, sp_window, avg, i_pi, i_sp, ratio, bg_pi, bg_sp}``:
    the two areas (intensity units x eV), their ratio (None when the σ+π area
    is not positive) and each window's background as ``[loss, level]`` pairs
    for drawing. A negative π area is returned as it is, never clipped (see
    :func:`ratio_flags`)."""
    try:
        avg = int(avg)
    except (TypeError, ValueError):
        return None
    if avg < 1:
        return None
    x, v = _sorted(loss, y)
    if x is None:
        return None
    a = _window_area(x, v, pi_window, avg)
    b = _window_area(x, v, sp_window, avg)
    if a is None or b is None:
        return None
    step = lambda n: max(1, n // BG_POINTS)
    return {"pi_window": [float(min(pi_window)), float(max(pi_window))],
            "sp_window": [float(min(sp_window)), float(max(sp_window))],
            "avg": avg, "i_pi": a[0], "i_sp": b[0],
            "ratio": a[0] / b[0] if b[0] > 0 else None,
            "bg_pi": [[float(p), float(q)] for p, q in
                      zip(a[1][::step(len(a[1]))], a[2][::step(len(a[2]))])],
            "bg_sp": [[float(p), float(q)] for p, q in
                      zip(b[1][::step(len(b[1]))], b[2][::step(len(b[2]))])]}


def ratio_flags(res):
    """What is wrong with a measured ratio, in words (an empty list when
    nothing is): a negative π area or a σ+π area that is not positive."""
    out = []
    if res.get("i_pi") is not None and \
            res["i_pi"] < -1e-9 * max(abs(res.get("i_sp") or 0.0), 1.0):
        out.append("the π area is negative (no peak above the background)")
    if res.get("i_sp") is not None and res["i_sp"] <= 0:
        out.append("the σ+π area is not positive, so there is no ratio")
    return out


def sp2_fraction(ratio, ref, zero=0.0):
    """The sp2 fraction on the scale where the reference (HOPG) ratio ``ref`` is
    1 and ``zero`` (diamond, 0 unless measured) is 0:
    ``(ratio - zero) / (ref - zero)``. The raw value; None when any input is
    missing or the reference does not lie above the zero."""
    try:
        ratio, ref, zero = float(ratio), float(ref), float(zero or 0.0)
    except (TypeError, ValueError):
        return None
    if not all(math.isfinite(v) for v in (ratio, ref, zero)) or ref <= zero:
        return None
    return (ratio - zero) / (ref - zero)


def _same_setting(res, ref):
    """None when a reference was measured with the windows and end-point
    channels of ``res`` (or does not say), else what differs."""
    if ref.get("pi_window") is None or ref.get("sp_window") is None:
        return None
    same = lambda a, b: all(abs(p - q) < 1e-9 for p, q in zip(a, b))
    if not (same(res["pi_window"], ref["pi_window"])
            and same(res["sp_window"], ref["sp_window"])
            and res["avg"] == ref.get("avg")):
        return "windows or end-point channels"
    return None


def sp2_for(res, refs):
    """The sp2 fraction of a stored ratio result ``res`` against the workbook's
    references ``refs`` (``{"hopg": {...}, "diamond": {...}}``, each with a
    ``ratio``, ``label`` and the windows and elastic energy it was measured
    with). Returns ``{fraction, problem, notes}``: ``fraction`` is None and
    ``problem`` says why when there is no reference or it was measured another
    way; ``notes`` flag a fraction outside the scale or an elastic-peak energy
    more than 1 % from the reference's."""
    out = {"fraction": None, "problem": "", "notes": []}
    hopg = (refs or {}).get("hopg")
    if not hopg or res.get("ratio") is None:
        return out
    dia = (refs or {}).get("diamond")
    for ref in (hopg, dia):
        diff = _same_setting(res, ref) if ref else None
        if diff:
            out["problem"] = (f"the {ref.get('label', 'reference')} was "
                              f"measured with other {diff}")
            return out
    zero = dia["ratio"] if dia else 0.0
    f = sp2_fraction(res["ratio"], hopg["ratio"], zero)
    if f is None:
        out["problem"] = "the reference ratio does not lie above the zero"
        return out
    out["fraction"] = f
    lo, hi = SCALE_LIMITS
    if f < lo or f > hi:
        out["notes"].append(f"the fraction is outside the {lo:g}–{hi:g} scale")
    for ref in (hopg, dia):
        e0, e1 = (ref or {}).get("e0"), res.get("elastic")
        if e0 and e1 and abs(e1 - e0) > ELASTIC_TOL * e0:
            out["notes"].append(
                f"the elastic peak ({e1:.0f} eV) differs from the "
                f"{ref.get('label', 'reference')}'s ({e0:.0f} eV) by more "
                f"than {ELASTIC_TOL * 100:g} %: the beam energy is not the "
                f"same")
    return out


def window_sensitivity(loss, y, pi_window=PI_WINDOW, sp_window=SP_WINDOW,
                       avg=EDGE_AVG, shift=SENSITIVITY_SHIFT):
    """How the ratio changes when each window limit is moved by ``shift`` eV
    either way: ``[{window, limit, shift, ratio, change_pct}]`` (empty when
    the ratio itself is not defined)."""
    base = pi_ratio(loss, y, pi_window, sp_window, avg)
    if not base or base["ratio"] is None or base["ratio"] == 0:
        return []
    out = []
    for name, win in (("π", list(pi_window)), ("σ+π", list(sp_window))):
        for k, limit in enumerate(("start", "end")):
            for sgn in (-1, 1):
                w = list(win)
                w[k] = w[k] + sgn * shift
                moved = (pi_ratio(loss, y, w, sp_window, avg) if name == "π"
                         else pi_ratio(loss, y, pi_window, w, avg))
                r = moved["ratio"] if moved else None
                out.append({"window": name, "limit": limit,
                            "shift": sgn * shift, "ratio": r,
                            "change_pct": (100.0 * (r - base["ratio"])
                                           / abs(base["ratio"])
                                           if r is not None else None)})
    return out


def ratio_text(res):
    """The windows and background of a ratio, in one sentence."""
    p, s = res["pi_window"], res["sp_window"]
    return (f"π {p[0]:g}–{p[1]:g} eV and σ+π {s[0]:g}–{s[1]:g} eV of loss, a "
            f"straight background between each window's limits (the mean of "
            f"{res['avg']} channels centred on each), trapezoidal areas")


def ratio_metadata(res, refs=None):
    """The metadata rows of a stored ratio (only what it holds)."""
    out = {"REELS π / (σ+π) windows": ratio_text(res)}
    notes = ratio_flags(res)
    if res.get("ratio") is not None:
        out["REELS π / (σ+π) ratio"] = f"{res['ratio']:.4f}"
    sp = sp2_for(res, refs)
    if sp["fraction"] is not None:
        hopg, dia = refs["hopg"], refs.get("diamond")
        out["REELS sp2 fraction"] = (
            f"{sp['fraction']:.2f} (HOPG {hopg['ratio']:.4f} = 1"
            + (f", {dia.get('label', 'diamond')} {dia['ratio']:.4f} = 0"
               if dia else "") + ")")
    elif sp["problem"]:
        out["REELS sp2 fraction"] = "not calculated: " + sp["problem"]
    notes += sp["notes"]
    if notes:
        out["REELS π / (σ+π) notes"] = "; ".join(notes)
    return out


def loss_fit_metadata(fit):
    """The metadata rows of a stored loss-region fit (only what it holds)."""
    out = {}
    lo, hi = fit["window"]
    by_key = {c["key"]: c for c in fit["components"]}
    text = (f"Gaussians at the loss energies of Stobinski et al. (2014), "
            f"{fit['background'].lower()} background, window {lo:g}–{hi:g} eV, "
            f"positions free within {fit['tol']:g} eV")
    out["REELS loss fit"] = text
    if "pi" in by_key and by_key["pi"]["area"] > 0:
        out["REELS π loss (eV)"] = f"{by_key['pi']['loss']:.2f}"
    if "sp2S" in by_key and by_key["sp2S"]["area"] > 0:
        out["REELS C sp2 surface loss (eV)"] = f"{by_key['sp2S']['loss']:.2f}"
    if fit.get("r") is not None:
        out["REELS sp2 bulk / sp2 surface (R)"] = f"{fit['r']:.2f}"
    g, ly = fit.get("geometry"), fit.get("layers")
    if g and ly and (ly.get("n") is not None or ly.get("at_least") is not None):
        n = (f"{ly['n']:.1f}" if ly.get("n") is not None
             else f"{ly['at_least']:.0f} or more")
        out["REELS layers (estimate)"] = (
            f"{n} (mean free path {g['lam']:g} nm, beam {g['a_in']:g}° in, "
            f"{g['a_out']:g}° out, {g['d']:g} nm a layer, graphite R "
            f"{g['r_graphite']:g}, information depth {ly['id']:.1f} nm)")
    return out


def information_depth(lam, a_in, a_out, p=99.0):
    """The depth (nm) from which ``p`` per cent of a REELS signal comes, in the
    straight-line approximation of Stobinski et al.: ``lam cos(a_in)
    cos(a_out) / (cos(a_in) + cos(a_out)) ln(1 / (1 - p/100))`` for the
    inelastic mean free path ``lam`` (nm) and the angles of the incident and
    emitted beams to the surface normal (degrees). None for an angle at or
    beyond 90 degrees, ``lam`` <= 0 or ``p`` outside (0, 100)."""
    try:
        ci, co = math.cos(math.radians(a_in)), math.cos(math.radians(a_out))
        lam, p = float(lam), float(p)
    except (TypeError, ValueError):
        return None
    if not (ci > 1e-9 and co > 1e-9 and lam > 0 and 0 < p < 100):
        return None
    return lam * ci * co / (ci + co) * math.log(1.0 / (1.0 - p / 100.0))


def layers(r, lam, a_in, a_out, d=LAYER_NM, r_graphite=GRAPHITE_R, p=99.0):
    """An estimate of the number of graphene layers from ``r`` = sp2 bulk /
    sp2 surface, or None when it is not defined.

    One layer has ``r`` = 0 (all surface) and graphite ``r_graphite``; between
    them the bulk signal grows as ``1 - exp(-z / (lam cos a_out))`` with ``z``
    the depth below the first layer, ``(n - 1) d``, and graphite is the
    thickness the information depth reaches (``n_graphite`` = ID / d):
    ``r / r_graphite = (1 - exp(-(n - 1) d / L)) / (1 - exp(-(n_g - 1) d / L))``.
    **This closed form is a reading of the curve of Fig. 7 of Stobinski et al.
    (J. Electron Spectrosc. Relat. Phenom. 195, 145, 2014), which they do not
    write as an equation.** With ``lam`` chosen so that the information
    depth is their 6.5 nm (2 keV, 0 and 35 degrees) it gives 7.9 and 4.1
    layers for their R = 2.72 and 1.55 (they state 8 and 4).

    Returns ``{n, at_least, id, n_graphite, L}``: ``n`` is None and
    ``at_least`` the graphite figure when ``r`` reaches ``r_graphite`` (the
    layer count is then beyond what the information depth can tell)."""
    depth = information_depth(lam, a_in, a_out, p)
    try:
        r, d, rg = float(r), float(d), float(r_graphite)
    except (TypeError, ValueError):
        return None
    if depth is None or not (math.isfinite(r) and r >= 0 and d > 0
                             and rg > 0):
        return None
    n_g = depth / d
    if n_g <= 1.0:
        return None
    big_l = float(lam) * math.cos(math.radians(a_out))
    out = {"n": None, "at_least": n_g, "id": depth, "n_graphite": n_g,
           "L": big_l}
    if r >= rg:
        return out
    frac = r / rg * (1.0 - math.exp(-(n_g - 1.0) * d / big_l))
    out["n"] = 1.0 - (big_l / d) * math.log(1.0 - frac)
    out["at_least"] = None
    return out


def tangent_points(result, p1, p2, span=1.25):
    """Two ``(loss, intensity)`` ends for drawing: from the baseline crossing
    to a little beyond the higher pick."""
    gap = result["gap"]
    x_end = max(p1[0], p2[0])
    x_end = gap + (x_end - gap) * span
    y_end = result["base"] + result["slope"] * (x_end - gap)
    return (gap, result["base"]), (x_end, y_end)
