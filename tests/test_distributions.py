import unittest
import opcsim
import pandas as pd
import os

from opcsim.distributions import *

_mode_          = (1000, 0.05, 1.25)
_mode2_         = (1000, 0.1, 1.5)

_mode_label     = 'ex_mode'
_mode_label2    = 'ex_mode2'

def fake_dist():
    dist = AerosolDistribution()

    dist.add_mode(_mode_, _mode_label)

    return dist

def fake_dist_2_modes():
    dist = AerosolDistribution()

    dist.add_mode(_mode_, _mode_label)
    dist.add_mode(_mode2_, _mode_label2)

    return dist

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        self._d_    = fake_dist()
        self._d2_   = fake_dist_2_modes()

    def tearDown(self):
        pass

    def test_add_modes(self):
        _tmp = AerosolDistribution()

        _tmp.add_mode(_mode_, _mode_label)

        self.assertEqual(len(_tmp.modes), 1)

    def test_get_modes(self):
        _tmp = AerosolDistribution()

        _tmp.add_mode(_mode_, _mode_label)

        data = _tmp._get_mode(_mode_label)

        self.assertEqual(data['N'], _mode_[0])

        bad_data = _tmp._get_mode('nopes')

        self.assertEqual(bad_data, None)

    def test_pdf(self):
        # Test to make sure the evaluation works for an individual diameter
        # Number-Weighted
        pdf     = self._d_.pdf(0.1)
        lnpdf   = self._d_.pdf(0.1, base = 'log')
        logpdf  = self._d_.pdf(0.1, base = 'log10')

        self.assertEqual(pdf * 0.1, lnpdf)
        self.assertEqual(pdf * 0.1 * np.log(10), logpdf)

        # Various Weights for Number-Weighted
        pdf_s   = self._d_.pdf(0.1, weight = 'surface')
        pdf_v   = self._d_.pdf(0.1, weight = 'volume')

        # Various Weights for log-weighted
        pdf_s   = self._d_.pdf(0.1, weight = 'surface', base = 'log')
        pdf_v   = self._d_.pdf(0.1, weight = 'volume', base = 'log')

        # Various Weights for log10-weighted
        pdf_s   = self._d_.pdf(0.1, weight = 'surface', base = 'log10')
        pdf_v   = self._d_.pdf(0.1, weight = 'volume', base = 'log10')

        with self.assertRaises(Exception):
            self._d_.pdf(0.1, weight = 'error')

        with self.assertRaises(Exception):
            self._d_.pdf(0.1, base = 'error')

    def test_cdf(self):
        # Test the cdf functionality for both eval and integrations
        cdf1    = self._d_.cdf(0.1)
        cdf     = self._d_.cdf(0.3)
        cdf2    = self._d_.cdf(Dmax = 0.3, Dmin = 0.1)
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
