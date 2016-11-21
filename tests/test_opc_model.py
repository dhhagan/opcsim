import unittest
import opcsim
import pandas as pd
import os

from opcsim.distributions import *
from opcsim.models import *

_mode_          = (1000, 0.15, 1.5)

_mode_label     = 'ex_mode'

def fake_dist():
    dist = AerosolDistribution()

    dist.add_mode(_mode_, _mode_label)

    return dist

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        self._d_    = fake_dist()

    def tearDown(self):
        pass

    def test_opc_model(self):
        res = OPC(num_bins = 1, dmin = 0.5, dmax = 2.5)

        self.assertEqual(res.num_bins, 1)

        res = OPC(num_bins = 2, dmin = 0.5, dmax = 2.5)

        self.assertEqual(res.num_bins, 2)

        # Make sure the distance in log space is equal for both bins
        self.assertEqual(np.log10(res.bins[0, 2]) - np.log10(res.bins[0, 0]),
                np.log10(res.bins[1, 2]) - np.log10(res.bins[1, 0]))

    def test_2_column_bins(self):
        bins = np.array([[0.5, 2.5]])

        res = OPC(bins = bins)

        self.assertEqual(res.num_bins, 1)
        self.assertEqual(res.dmin, 0.5)
        self.assertEqual(res.dmax, 2.5)
        self.assertEqual(round(res.bins[0, 1], 3), 1.118)

    def test_3_column_bins(self):
        bins = np.array([[0.5, 1.18, 2.5], [2.5, 5, 10.0]])

        res = OPC(bins = bins)

        self.assertEqual(res.num_bins, 2)
        self.assertEqual(res.dmin, 0.5)
        self.assertEqual(res.dmax, 10)
        self.assertEqual(round(res.bins[1, 1], 3), 5.000)

    def test_histogram(self):
        res     = OPC(dmin = 0.5, dmax = 2.5, num_bins = 1)

        hist_n  = res.histogram(distribution = self._d_, weight = 'number', base = 'log10')
        hist_s  = res.histogram(distribution = self._d_, weight = 'surface', base = 'log10')
        hist_v  = res.histogram(distribution = self._d_, weight = 'volume', base = 'log10')

        hist_n2 = res.histogram(distribution = self._d_, weight = 'number', base = 'log')
        hist_s2 = res.histogram(distribution = self._d_, weight = 'surface', base = 'log')
        hist_v2 = res.histogram(distribution = self._d_, weight = 'volume', base = 'log')

        self.assertEqual(len(hist_n), 1)
        self.assertNotEqual(hist_n, hist_s)
        self.assertNotEqual(hist_n, hist_v)
        self.assertNotEqual(hist_n, hist_n2)
        self.assertNotEqual(hist_s, hist_s2)
        self.assertNotEqual(hist_v, hist_v2)

        res     = OPC(dmin = 0.5, dmax = 2.5, num_bins = 5)

        hist_n  = res.histogram(distribution = self._d_, weight = 'number', base = 'log10')
        hist_s  = res.histogram(distribution = self._d_, weight = 'surface', base = 'log10')
        hist_v  = res.histogram(distribution = self._d_, weight = 'volume', base = 'log10')

        hist_n2 = res.histogram(distribution = self._d_, weight = 'number', base = 'log')
        hist_s2 = res.histogram(distribution = self._d_, weight = 'surface', base = 'log')
        hist_v2 = res.histogram(distribution = self._d_, weight = 'volume', base = 'log')

        self.assertEqual(len(hist_n), 5)
        self.assertNotEqual(hist_n[0], hist_s[0])

    def test_evaluate(self):
        res = OPC(dmin = 0.5, dmax = 2.5, num_bins = 2)

        # Try to evaluate
        ans     = res.evaluate(self._d_)
        ansb    = res.evaluate(self._d_, param = 'na')
        ansc    = res.evaluate(self._d_, param = 'vm')
        ansd    = res.evaluate(self._d_, param = 'va')
        anse    = res.evaluate(self._d_, param = 'nm_total')
        ansf    = res.evaluate(self._d_, dmax = 2.5, param = 'na_total')
        ansg    = res.evaluate(self._d_, param = 'vm_total')
        ansh    = res.evaluate(self._d_, dmax = 2.5, param = 'va_total')

        nm_na   = res.evaluate(self._d_, param = 'nm_na')
        nm_va   = res.evaluate(self._d_, param = 'nm_va')
        vm_va   = res.evaluate(self._d_, param = 'vm_va')

        self.assertEqual(np.sum(ans), anse)
        self.assertGreater(ansf, np.sum(ansb))
        self.assertEqual(np.sum(ansc), ansg)
        self.assertGreater(ansh, np.sum(ansd))

        # Force error
        with self.assertRaises(Exception):
            res.evaluate(self._d_, param = 'asda')

        # Test errant CE equation
        def bad():
            return 1

        new = OPC(num_bins = 3, ce = bad)

    def test_set_bins(self):
        res = OPC(dmin = 0.5, dmax = 2.5, num_bins = 1)

        self.assertEqual(round(res.bins[0, 2] - res.bins[0, 0], 2), 2.0)
        self.assertEqual(res.bins[0, 0], 0.5)
        self.assertEqual(round(res.bins[0, 2], 2), 2.5)

    def test_bins(self):
        res = OPC(dmin = 0.5, dmax = 2.5, num_bins = 1)

        l, h, w   = res.boxes(self._d_)

        self.assertEqual(round(w[0], 2), 2.0)
        self.assertEqual(round(l[0], 2), 0.5)
