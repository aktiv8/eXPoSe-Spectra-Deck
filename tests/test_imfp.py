"""imfp.py: the TPP-2M average-matrix IMFP and the KE^0.6 Avantage-style
proxy used by quant.py's "scofield_tpp2m"/"scofield_ke06" RSF-fallback
tiers.

Run:  python -m unittest discover tests
"""

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import imfp  # noqa: E402


class TestImfpNm(unittest.TestCase):
    # hand-computed from the same formula independently, at the module's
    # fixed average-matrix constants -- see imfp.py's own docstring
    REFERENCE = {50: 0.5066507922262123, 100: 0.577623423099999,
                200: 0.7947830339179165, 500: 1.4263303268788368,
                1000: 2.368001423854056, 2000: 4.060332098481516}

    def test_reference_values(self):
        for ke, expected in self.REFERENCE.items():
            self.assertAlmostEqual(imfp.imfp_nm(ke), expected, places=9,
                                   msg=f"KE={ke}")

    def test_increases_with_kinetic_energy(self):
        values = [imfp.imfp_nm(ke) for ke in
                 (50, 100, 200, 500, 1000, 2000)]
        self.assertEqual(values, sorted(values))

    def test_just_inside_the_valid_range(self):
        self.assertIsNotNone(imfp.imfp_nm(50.0))
        self.assertIsNotNone(imfp.imfp_nm(2000.0))

    def test_outside_the_valid_range_gives_none(self):
        self.assertIsNone(imfp.imfp_nm(49.9))
        self.assertIsNone(imfp.imfp_nm(2000.1))
        self.assertIsNone(imfp.imfp_nm(0))
        self.assertIsNone(imfp.imfp_nm(None))


class TestKePowerFactor(unittest.TestCase):
    def test_matches_the_plain_power_law(self):
        for ke in (50, 200, 954.7, 1419, 2000):
            self.assertAlmostEqual(imfp.ke_power_factor(ke), ke ** 0.6,
                                   places=9, msg=f"KE={ke}")

    def test_no_upper_or_lower_bound_unlike_tpp2m(self):
        # not a published-range-restricted formula the way TPP-2M is
        self.assertIsNotNone(imfp.ke_power_factor(10))
        self.assertIsNotNone(imfp.ke_power_factor(5000))

    def test_non_positive_or_missing_gives_none(self):
        self.assertIsNone(imfp.ke_power_factor(0))
        self.assertIsNone(imfp.ke_power_factor(-5))
        self.assertIsNone(imfp.ke_power_factor(None))


if __name__ == "__main__":
    unittest.main()


class TestCasaEscapeDepth(unittest.TestCase):
    """CasaXPS's own escape-depth types, as read from CasaXPS.exe (Ghidra,
    function 0x1403cc280). Each expectation is written out from the
    decompiled arithmetic, not from imfp.py."""

    def test_tpp2m_matrix_and_no_energy_limit(self):
        import math
        u = 4.0 * 3.513 / 12.01
        ep2 = 829.4 * u
        beta = -0.1 + 0.944 / math.sqrt(ep2 + 5.4 ** 2) + 0.069 * 3.513 ** 0.1
        gamma = 0.191 / math.sqrt(3.513)
        for ke in (30.0, 950.0, 1420.0, 5000.0):
            den = ep2 * (beta * math.log(gamma * ke) - (1.97 - 0.91 * u) / ke
                         + (53.4 - 20.8 * u) / ke ** 2)
            self.assertAlmostEqual(imfp.casa_imfp_nm(ke), ke / den / 10,
                                   places=12)
        self.assertIsNone(imfp.casa_imfp_nm(0))
        self.assertIsNone(imfp.casa_imfp_nm(None))

    def test_relative_factor_agrees_with_the_average_matrix(self):
        # the absolute values differ by ~19 %, a ratio between two lines by
        # 0.02 % -- all that atomic percent depends on
        a = imfp.casa_imfp_nm(1420) / imfp.casa_imfp_nm(950)
        b = imfp.imfp_nm(1420) / imfp.imfp_nm(950)
        self.assertAlmostEqual(a / b, 1.0, delta=5e-4)
        self.assertGreater(imfp.casa_imfp_nm(950) / imfp.imfp_nm(950), 1.15)

    def test_power_law_uses_the_stored_exponent(self):
        self.assertAlmostEqual(imfp.casa_escape_factor(1, 400, 100, -0.6),
                               4 ** -0.6)
        self.assertEqual(imfp.casa_escape_factor(1, 400, 100, 0.0), 1.0)

    def test_eal_and_imfp_ratios(self):
        eal = lambda e: 0.65 + 0.007 * e ** 0.93
        self.assertAlmostEqual(imfp.casa_escape_factor(2, 900, 1), eal(1) / eal(900))
        line = lambda e: 0.73 + 0.0095 * e ** 0.872
        self.assertAlmostEqual(imfp.casa_escape_factor(3, 900, 1), line(1) / line(900))

    def test_tpp2m_factor_is_the_inverse_length_ratio(self):
        self.assertAlmostEqual(
            imfp.casa_escape_factor(5, 1420, 950),
            imfp.casa_imfp_nm(950) / imfp.casa_imfp_nm(1420), places=12)

    def test_overlayer_types(self):
        import math
        eal = lambda e: 0.65 + 0.007 * e ** 0.93
        l, lr, t = eal(900), eal(300), 2.0
        c = 4.0 ** 0.38
        self.assertAlmostEqual(imfp.casa_escape_factor(4, 900, 300, t),
                               lr / l * math.exp(t * c / l - t * c / lr))
        self.assertAlmostEqual(
            imfp.casa_escape_factor(6, 900, 300, t),
            lr * (1 - math.exp(-t / lr)) / (l * (1 - math.exp(-t / l))))

    def test_none_and_unknown_types_do_nothing(self):
        self.assertEqual(imfp.casa_escape_factor(0, 900, 1, -0.6), 1.0)
        self.assertEqual(imfp.casa_escape_factor(9, 900, 1, -0.6), 1.0)
