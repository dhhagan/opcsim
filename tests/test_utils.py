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
    
    def test_k_eff(self):

        # calculate easy weighting
        k1 = opcsim.utils.k_eff(kappas=[0.5, 1.5], weights=[1, 0])
        self.assertEqual(k1, 0.5)

        # average
        k1 = opcsim.utils.k_eff(kappas=[0, 1], weights=[0.5, 0.5])
        self.assertEqual(k1, 0.5)

        # compute based on diameters
        k1 = opcsim.utils.k_eff(kappas=[0, 1], diams=[1., 2.])
        self.assertEqual(k1, 8./9.)
    
    def test_ri_eff(self):
        
        # calculate easy weighting
        r1 = opcsim.utils.ri_eff(
            species=[complex(1.5, 0), complex(0, 0)], weights=[1, 0])
        self.assertEqual(r1, complex(1.5, 0))

        # calculate everage
        r1 = opcsim.utils.ri_eff(
            species=[complex(1.5, 0), complex(2.0, 1.0)], weights=[0.5, 0.5])
        self.assertEqual(r1, complex(1.75, 0.5))

    def power_law_fit(self):
        xs = np.linspace(1, 10, 10)
        a, b = 2, 0.1
        y = [x*np.random.uniform(0.99, 1.01) for x in a*np.power(xs, b)]

        fitted = opcsim.utils.power_law_fit(xs, cscat=y)
        self.assertAlmostEqual(y.sum(), fitted.sum())

        # test with weights
        fitted = opcsim.utils.power_law_fit(xs, cscat=y, fit_kws=dict(sigma=np.power(y, 10)))
        self.assertAlmostEqual(y.sum(), fitted.sum())

    
    def test_squash_dips(self):
        xs = np.array([1, 2, 3, 4, 3, 6, 7, 8, 9, 10])
        s = opcsim.utils.squash_dips(xs)

        self.assertTrue((np.diff(s) < 0).any() == False)

