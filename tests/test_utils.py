import unittest
import opcsim
import numpy as np
import random

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_make_bins(self):
        dmin = 0.5
        dmax = 10.
        n_bins = 2

        bins = opcsim.utils.make_bins(dmin, dmax, n_bins)

        self.assertEqual(bins.shape[0], n_bins)
        self.assertEqual(bins[0, 0], dmin)
        self.assertEqual(bins[-1, -1], dmax)

        ba = opcsim.utils.make_bins(dmin, dmax, n_bins, base=None)

        # Make sure the bin midpoints are correct
        self.assertEqual(np.mean([ba[0,0], ba[0,2]]), ba[0,1])
        self.assertGreater(ba[0,1], bins[0,1])

    def test_midpoints(self):
        dmin = 0.5
        dmax = 17.

        arr = np.array([[dmin, 2.5], [2.5, 10], [10, dmax]])

        bins = opcsim.utils.midpoints(arr)

        self.assertEqual(bins[0, 0], dmin)
        self.assertEqual(bins[-1, -1], dmax)

        with self.assertRaises(ValueError):
            bins = opcsim.utils.midpoints(np.array([[1, 2, 3, 4]]))
