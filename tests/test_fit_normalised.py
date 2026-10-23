"""A CasaXPS fit keeps its components when the spectrum is normalised
(the visibility floor of a component once used the unnormalised peak and
blanked every point)."""

import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    HAVE_MPL = True
except Exception:                                   # pragma: no cover
    HAVE_MPL = False


@unittest.skipUnless(HAVE_MPL, "matplotlib not available")
class TestFitNormalised(unittest.TestCase):
    def test_components_survive_every_normalisation(self):
        import casafit
        import plots
        from test_casafit import fitted_region
        r = fitted_region(n=61)
        cv = casafit.curves(r.fit, r.energy, r.counts, r.photon_energy,
                            *r.dwell_and_scans())
        raw_peak = max(v for c in cv for _k, vals in c.components
                       for v in vals if v == v)
        for norm in ("None", "Max = 1", "Area = 1"):
            ax = Figure().add_subplot()
            fit = {"curves": cv, "colours": ["#111", "#c00", "#0a0", "#00c"],
                   "show": {"envelope": True, "components": True,
                            "background": True}}
            plots.draw_stack(ax, [r], norm=norm, fit=fit)
            tops = [max(v for v in line.get_ydata() if v == v)
                    for line in ax.lines
                    if any(v == v for v in line.get_ydata())]
            # data, background, component (twice: top line and outline), envelope
            self.assertEqual(len(tops), len(ax.lines), norm)
            nf = plots.norm_factor(r, norm)
            self.assertTrue(any(t >= 0.5 * raw_peak / nf for t in tops), norm)


if __name__ == "__main__":
    unittest.main()
