import unittest
import opcsim
import pandas as pd
import numpy as np
import os
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

        dps = np.linspace(0.01, 1., 100)

        # Test to make sure the evaluation works for an individual diameter
        # Number-Weighted
        # Evaluate at 0.1 microns
        pdf     = d.pdf(0.1)

        # Evaluate an array
        pdf_arr = d.pdf(dps)

        # Evaluate for a single mode
        pdf_1   = d.pdf(0.1, mode="Mode I")

        self.assertGreaterEqual(pdf, 0.0)
        self.assertGreaterEqual(pdf, pdf_1)
        self.assertEqual(len(pdf_arr), len(dps))

        #lnpdf   = self._d_.pdf(0.1, base = 'log')
        #logpdf  = self._d_.pdf(0.1, base = 'log10')

        #self.assertEqual(pdf * 0.1, lnpdf)
        #self.assertEqual(pdf * 0.1 * np.log(10), logpdf)

        # Various Weights for Number-Weighted
        #pdf_s   = self._d_.pdf(0.1, weight = 'surface')
        #pdf_v   = self._d_.pdf(0.1, weight = 'volume')

        # Various Weights for log-weighted
        #pdf_s   = self._d_.pdf(0.1, weight = 'surface', base = 'log')
        #pdf_v   = self._d_.pdf(0.1, weight = 'volume', base = 'log')

        # Various Weights for log10-weighted
        #pdf_s   = self._d_.pdf(0.1, weight = 'surface', base = 'log10')
        #pdf_v   = self._d_.pdf(0.1, weight = 'volume', base = 'log10')

        #pdf     = self._d_.pdf(0.1, weight = 'surface', base = 'log10', mode = _mode_label)

        #with self.assertRaises(Exception):
        #    self._d_.pdf(0.1, weight = 'error')

        #with self.assertRaises(Exception):
        #    self._d_.pdf(0.1, base = 'error')

    """
    def test_cdf(self):
        # Test the cdf functionality for both eval and integrations
        cdf1    = self._d_.cdf(0.1)
        cdf     = self._d_.cdf(0.3)
        cdf2    = self._d_.cdf(dmax = 0.3, dmin = 0.1)
        cdf_s   = self._d_.cdf(0.1, weight = 'surface')
        cdf_v   = self._d_.cdf(0.1, weight = 'volume')

        self.assertEqual(cdf2, cdf - cdf1)
        self.assertGreater(cdf, cdf2)

        with self.assertRaises(Exception):
            self._d_.cdf(0.1, weight = 'error')

        res = self._d_.__repr__()

        # Test an individual mode
        cdf     = self._d_.cdf(0.1, mode = _mode_label)

        self.assertEqual(cdf, cdf1)
    """
    """

    def test_dist_mean(self):
        # Test the mean diameter
        mean_d  = self._d_.mean(weight = 'number', diameter = True)
        mean_db = self._d_.mean(weight = 'number', diameter = False)

        mean_sd = self._d_.mean(weight = 'surface', diameter = True)
        mean_s  = self._d_.mean(weight = 'surface', diameter = False)

        mean_vd = self._d_.mean(weight = 'volume', diameter = True)
        mean_v  = self._d_.mean(weight = 'volume', diameter = False)

        self.assertEqual(mean_d, mean_db)

        # Test bad weight
        with self.assertRaises(Exception):
            self._d_.mean(weight = 'error')

        # Test getting just one node
        mean    = self._d2_.mean(weight = 'number', mode = _mode_label)

        self.assertEqual(mean, mean_d)

    def test_dist_median(self):
        med     = self._d_.median(weight = 'number')
        med_s   = self._d_.median(weight = 'surface')
        med_v   = self._d_.median(weight = 'volume')

        self.assertGreater(med_s, med)
        self.assertGreater(med_v, med_s)

        # Test bad weighting
        with self.assertRaises(Exception):
            self._d_.median(weight = 'error')

        # Test first mode
        med2    = self._d2_.median(weight = 'number', mode = _mode_label)

        self.assertEqual(med2, med)

    """
