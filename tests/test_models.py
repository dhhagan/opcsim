import unittest
import opcsim
import pandas as pd
import numpy as np

from opcsim.distributions import *
from opcsim.models import *

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_opc_model(self):
        n_bins = 10
        dmin = 0.3
        dmax = 10.0
        wl = 0.658
        theta = (32.0, 88.0)

        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)

        self.assertIsInstance(opc, opcsim.OPC)
        self.assertEqual(opc.n_bins, 10)
        self.assertEqual(opc.dmin, dmin)
        self.assertEqual(opc.dmax, dmax)
        self.assertEqual(opc.wl, wl)
        self.assertEqual(opc.theta, theta)
        self.assertEqual(opc.label, None)
        self.assertEqual(opc.calibration_function, None)
        self.assertEqual(len(opc.bin_boundaries), n_bins+1)
        self.assertTrue((opc.dlogdp == np.log10(
            opc.bins[:, -1]) - np.log10(opc.bins[:, 0])).all())
        self.assertTrue((opc.ddp == opc.bins[:, -1] - opc.bins[:, 0]).all())

        # test for bins manually set as list of values (boundaries)
        bins = [0.38, 0.54, 0.78, 1.05, 1.5, 2.5, 3.5, 5.]
        opc = opcsim.OPC(wl=wl, bins=bins, theta=theta)

        # test for bins manually set as 3xn array
        self.assertIsInstance(opc, opcsim.OPC)
        self.assertEqual(opc.n_bins, len(bins)-1)
        self.assertEqual(opc.dmin, bins[0])
        self.assertEqual(opc.dmax, bins[-1])
        self.assertEqual(opc.wl, wl)
        self.assertEqual(opc.theta, theta)
        self.assertEqual(opc.bins.shape[0], opc.n_bins)
    
    def test_opc_calibrate(self):
        n_bins = 10
        dmin = 0.3
        dmax = 10.0
        wl = 0.658
        theta = (32.0, 88.0)

        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)
        
        # calibrate the OPC for PSL's
        self.assertIsNone(opc.calibration_function)
        opc.calibrate(material="psl")
        self.assertIsNotNone(opc.calibration_function)

        # test the calibration for a spec'd material
        # create an opc based on the above

        # calibrate the OPC for a random RI
        opc.calibrate(material=complex(1.5, 0))
        self.assertIsNotNone(opc.calibration_function)

        # try fitting the data this time...
        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)

        # calibrate the OPC for PSL's
        self.assertIsNone(opc.calibration_function)
        opc.calibrate(material="psl", method="linear")
        self.assertIsNotNone(opc.calibration_function)

        # try for an integer
        opc.calibrate(material=1.5)
        self.assertIsNotNone(opc.calibration_function)

        # calibrate the OPC for PSL's
        opc.calibrate(material="psl", method="piecewise")
        self.assertIsNotNone(opc.calibration_function)

        # try for a bad string
        with self.assertRaises(ValueError):
            opc.calibrate(material="random_thing")

        # non-existent calibration method
        with self.assertRaises(ValueError):
            opc.calibrate("psl", method="wrong")
    
    def test_opc_evaluate(self):
        n_bins = 10
        dmin = 0.3
        dmax = 10.0
        wl = 0.658
        theta = (32.0, 88.0)

        # build a distribution
        d = opcsim.AerosolDistribution()
        d.add_mode(n=1e3, gm=0.4, gsd=1.5, rho=1.6, refr=complex(1.5, 0))

        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)

        with self.assertRaises(Exception):
            h = opc.histogram(d)

        # calibrate the OPC for PSL's
        self.assertIsNone(opc.calibration_function)
        opc.calibrate(material="psl")
        self.assertIsNotNone(opc.calibration_function)

        # test the histogram
        h = opc.evaluate(d)

    def test_opc_histogram(self):
        n_bins = 10
        dmin = 0.3
        dmax = 10.0
        wl = 0.658
        theta = (32.0, 88.0)

        # build a distribution
        d = opcsim.AerosolDistribution()
        d.add_mode(n=1e3, gm=0.4, gsd=1.5, rho=1.6, refr=complex(1.5, 0))

        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)

        opc.calibrate(material="psl")
        self.assertIsNotNone(opc.calibration_function)

        # create the histogram
        h = opc.histogram(d, weight='number')

        self.assertEqual(len(h), opc.n_bins)

        # calculate surface area, volume, and mass distributions
        h = opc.histogram(d, weight='surface')
        h = opc.histogram(d, weight='volume')
        h = opc.histogram(d, weight='mass', rho=1.65)

        # test dN/dDp
        h = opc.histogram(d, weight="number", base=None)

        # force error
        with self.assertRaises(ValueError):
            h = opc.histogram(d, weight="unknown")

    def test_opc_integrate(self):
        n_bins = 10
        dmin = 0.3
        dmax = 10.0
        wl = 0.658
        theta = (32.0, 88.0)

        # build a distribution
        d = opcsim.AerosolDistribution()
        d.add_mode(n=1e3, gm=0.4, gsd=1.5, rho=1.6, refr=complex(1.5, 0))

        # create an opc based on the above
        opc = opcsim.OPC(wl=wl, n_bins=n_bins, dmin=dmin,
                         dmax=dmax, theta=theta)
        
        opc.calibrate(material="psl")
        self.assertIsNotNone(opc.calibration_function)

        # integrate in number space
        n1 = opc.integrate(d, dmin=0, dmax=1., weight="number")
        n2 = opc.integrate(d, dmin=0., dmax=2.5, weight="number")
        n3 = opc.integrate(d, dmin=0., dmax=10., weight="number")

        self.assertGreaterEqual(n2, n1)
        self.assertGreaterEqual(n3, n1)
        self.assertGreaterEqual(n3, n2)

        # force a value error
        with self.assertRaises(ValueError):
            n1 = opc.integrate(d, dmin=0., dmax=1., weight="bad-weight")
    
        # test in-between bounds
        n1 = opc.integrate(d, dmin=.6, dmax=1., weight="number")
        n2 = opc.integrate(d, dmin=.61, dmax=.62, weight="number")

        self.assertGreater(n1, n2)

        # integrate surface area, volume, and mass distributions
        n1 = opc.integrate(d, dmin=0., dmax=1., weight="surface")
        n2 = opc.integrate(d, dmin=0., dmax=1., weight="volume")
        n3 = opc.integrate(d, dmin=0., dmax=1., weight="mass", rho=1.5)

    def test_nephelometer(self):
        neph = opcsim.Nephelometer(wl=0.658, theta=(7., 173.))

        self.assertIsNone(neph.pm1_ratio)
        self.assertIsNone(neph.pm25_ratio)
        self.assertIsNone(neph.pm10_ratio)

        # calibrate the device to a distribution
        d = opcsim.AerosolDistribution()
        d.add_mode(n=1000, gm=.2, gsd=1.5, kappa=0.53, refr=complex(1.592, 0), rho=1.77)

        neph.calibrate(d, rh=0.)

        self.assertIsNotNone(neph.pm1_ratio)
        self.assertIsNotNone(neph.pm25_ratio)
        self.assertIsNotNone(neph.pm10_ratio)

        # test evaluate functionality
        vals = neph.evaluate(d, rh=0.)
        vals2 = neph.evaluate(d, rh=95.)

        self.assertGreaterEqual(vals2[1], vals[1])
        self.assertGreaterEqual(vals2[2], vals[2])
        self.assertGreaterEqual(vals2[3], vals[3])

