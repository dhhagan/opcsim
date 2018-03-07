import unittest
import opcsim
import numpy as np
import random

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_add_modes(self):
        tmp = opcsim.AerosolDistribution()

        n, gm, gsd = random.random(), random.random(), random.random()

        tmp.add_mode(n, gm, gsd, "label")

        self.assertEqual(len(tmp.modes), 1)

    def test_get_modes(self):
        tmp = opcsim.AerosolDistribution()

        label = 'mode1'
        n, gm, gsd = random.random(), random.random(), random.random()

        tmp.add_mode(n, gm, gsd, label)

        # Get the mode
        m = tmp._get_mode(label)

        self.assertEqual(m['N'], n)

        # Test gettin a non-existing mode
        ne = tmp._get_mode('none')

        self.assertEqual(ne, None)

    def test_pdf(self):
        # Use a sample Urban distribution
        d = opcsim.load_distribution("Urban")

        dp = 0.1

        dps = np.linspace(0.01, 1., 100)

        # Test to make sure the evaluation works for an individual diameter
        # Number-Weighted
        pdf = d.pdf(dp, base=None)

        # Evaluate an array
        pdf_arr = d.pdf(dps)

        # Evaluate for a single mode
        pdf_1   = d.pdf(dp, mode="Mode I")

        self.assertGreaterEqual(pdf, 0.0)
        self.assertGreaterEqual(pdf, pdf_1)
        self.assertEqual(len(pdf_arr), len(dps))

        # Check the log and log10 weighted versions
        pdf_log = d.pdf(dp, base='log')
        pdf_log10 = d.pdf(dp, base='log10')

        # Make sure the pdf functions of various weights are correct
        self.assertEqual(round(pdf * dp, 3), round(pdf_log, 3))
        self.assertEqual(round(pdf * dp * np.log(10), 3), round(pdf_log10, 3))

        # Various Weights for Number-Weighted
        pdf_s = d.pdf(dp, weight='surface', base=None)
        pdf_v = d.pdf(dp, weight='volume', base=None)
        pdf_m = d.pdf(dp, weight='mass', base=None)

        # Various Weights for log-weighted
        pdf_s = d.pdf(dp, weight='surface', base='log')
        pdf_v = d.pdf(dp, weight='volume', base='log')
        pdf_m = d.pdf(dp, weight='mass', base='log')

        # Various Weights for log10-weighted
        pdf_s = d.pdf(dp, weight='surface', base='log10')
        pdf_v = d.pdf(dp, weight='volume', base='log10')
        pdf_m = d.pdf(dp, weight='mass', base='log10')

        with self.assertRaises(Exception):
            d.pdf(dp, weight='error')

        with self.assertRaises(Exception):
            d.pdf(dp, base='error')

    def test_cdf_number(self):
        # Use a sample Urban distribution
        d = opcsim.load_distribution("Urban")

        cdf_1 = d.cdf(dmax=1.0)
        cdf_25 = d.cdf(dmax=2.5)
        cdf_diff = d.cdf(dmin=1.0, dmax=2.5)

        self.assertGreaterEqual(cdf_25, cdf_1)
        self.assertEqual(round(cdf_diff, 3), round(cdf_25 - cdf_1, 3))

        with self.assertRaises(Exception):
            d.cdf(0.1, weight='error')

    def test_cdf_surface(self):
        d = opcsim.load_distribution("Urban")

        # Test the surface area weighted versions
        cdf_sa = d.cdf(dmax=1.0, weight='surface')
        cdf_sa2 = d.cdf(dmax=2.5, weight='surface')
        cdf_sa_diff = d.cdf(dmin=1.0, dmax=2.5, weight='surface')

        self.assertGreaterEqual(cdf_sa2, cdf_sa)
        self.assertEqual(round(cdf_sa_diff, 3), round(cdf_sa2 - cdf_sa, 3))

    def test_cdf_volume(self):
        d = opcsim.load_distribution("Urban")

        # Test the surface area weighted versions
        cdf_v = d.cdf(dmax=1.0, weight='volume')
        cdf_v2 = d.cdf(dmax=2.5, weight='volume')
        cdf_v_diff = d.cdf(dmin=1.0, dmax=2.5, weight='volume')

        self.assertGreaterEqual(cdf_v2, cdf_v)
        self.assertEqual(round(cdf_v_diff, 3), round(cdf_v2 - cdf_v, 3))

        with self.assertRaises(ValueError):
            d.cdf(dmin=2.5, dmax=1.)

        # Test a single mode
        cdf = d.cdf(dmax=2.5, mode='Mode I')

        self.assertIsNotNone(cdf)

    def test_cdf_mass(self):
        d = opcsim.load_distribution("Urban")

        # Test the mass weighted versions
        cdf_m = d.cdf(dmax=1.0, weight='mass')
        cdf_m2 = d.cdf(dmax=2.5, weight='mass')
        cdf_m_diff = d.cdf(dmin=1.0, dmax=2.5, weight='mass')

        self.assertGreaterEqual(cdf_m2, cdf_m)
        self.assertEqual(round(cdf_m_diff, 3), round(cdf_m2 - cdf_m, 3))

        # Make sure the mass is roughly correct
        self.assertGreaterEqual(cdf_m, 3.)
        self.assertLessEqual(cdf_m, 20.)

    def test_bad_distribution(self):
        with self.assertRaises(ValueError):
            d = opcsim.load_distribution("None")

    def test_repr(self):
        d = opcsim.load_distribution("Urban")

        self.assertTrue(repr(d) == "AerosolDistribution: urban")
