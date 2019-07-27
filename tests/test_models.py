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
        opc.calibrate(material="psl", method="smooth")
        self.assertIsNotNone(opc.calibration_function)

        # test the calibration for a spec'd material
        # create an opc based on the above

        # calibrate the OPC for a random RI
        opc.calibrate(material=complex(1.5, 0), method="smooth")
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
        opc.calibrate(material=1.5, method="smooth")
        self.assertIsNotNone(opc.calibration_function)

        # try for a bad string
        with self.assertRaises(ValueError):
            opc.calibrate(material="random_thing", method="smooth")

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
        opc.calibrate(material="psl", method="smooth")
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

        opc.calibrate(material="psl", method="smooth")
        self.assertIsNotNone(opc.calibration_function)

        # create the histogram
        lb, h, ddp = opc.histogram(d, weight='number')

        self.assertTrue((lb == opc.bins[:, 0]).all())
        self.assertTrue((ddp == opc.ddp).all())
        self.assertEqual(len(h), opc.n_bins)

        # calculate surface area, volume, and mass distributions
        lb, h, ddp = opc.histogram(d, weight='surface')
        lb, h, ddp = opc.histogram(d, weight='volume')
        lb, h, ddp = opc.histogram(d, weight='mass', rho=1.65)

        # test dN/dDp
        lb, h, ddp = opc.histogram(d, weight="number", base=None)

        # force error
        with self.assertRaises(ValueError):
            lb, h, ddp = opc.histogram(d, weight="unknown")

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
        
        opc.calibrate(material="psl", method="smooth")
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
