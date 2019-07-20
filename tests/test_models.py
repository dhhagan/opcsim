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

    # def test_opc_model(self):
    #     n_bins = 2
    #     dmin = 0.5
    #     dmax = 2.5

    #     opc = opcsim.models.OPC(n_bins=n_bins, dmin=dmin, dmax=dmax)

    #     self.assertEqual(dmin, opc.dmin)
    #     self.assertEqual(dmax, opc.dmax)
    #     self.assertEqual(n_bins, opc.n_bins)

    # def test_opc_with_bins(self):
    #     bins = np.array([[0.5, 2.5], [2.5, 10.]])

    #     with self.assertRaises(Exception):
    #         opc = opcsim.models.OPC(bins=bins)

    #     bins = opcsim.utils.midpoints(bins)

    #     opc = opcsim.models.OPC(bins=bins)

    #     self.assertEqual(opc.n_bins, bins.shape[0])

    #     self.assertIsNotNone(opc.dlogdp)

    # def test_opc_evaluate(self):
    #     opc = opcsim.models.OPC()
    #     d = opcsim.load_distribution("Urban")

    #     # Test dNdDp
    #     dNdDp = opc.evaluate(d, weight='number', base=None)
    #     dNdlnDp = opc.evaluate(d, weight='number', base='log')
    #     dNdlogDp = opc.evaluate(d, weight='number', base='log10')

    #     self.assertEqual(dNdDp.shape[0], opc.n_bins)
    #     self.assertEqual(dNdlnDp.shape[0], opc.n_bins)
    #     self.assertEqual(dNdlogDp.shape[0], opc.n_bins)

    #     # Test bad distribution
    #     with self.assertRaises(Exception):
    #         opc.evaluate(1)

    # def test_histogram(self):
    #     opc = opcsim.models.OPC()
    #     d = opcsim.load_distribution("Urban")

    #     bb, h, bw = opc.histogram(d)

    #     self.assertEqual(len(bb), len(h))
    #     self.assertEqual(len(bb), len(bw))

    # def test_opc_number(self):
    #     opc = opcsim.models.OPC()
    #     d = opcsim.load_distribution("Urban")

    #     vm = opc.number(d)
    #     va = opc.number(d, measured=False)

    #     # They should be equal because we have 100% efficiency
    #     self.assertEqual(vm, va)
    #     self.assertEqual(len(vm), opc.n_bins)

    # def test_opc_sa(self):
    #     opc = opcsim.models.OPC()
    #     d = opcsim.load_distribution("Urban")

    #     vm = opc.surface_area(d, dmax=5)
    #     va = d.cdf(dmin=0, dmax=5, weight='surface')

    #     self.assertNotEqual(vm, va)
    #     self.assertEqual(len(vm), opc.n_bins)

    # def test_opc_volume(self):
    #     opc = opcsim.models.OPC()
    #     d = opcsim.load_distribution("Urban")

    #     vm = opc.volume(d)
    #     va = d.cdf(dmin=0, dmax=5, weight='volume')

    #     self.assertNotEqual(vm, va)
    #     self.assertEqual(len(vm), opc.n_bins)

    # def test_counting_efficiency_simple(self):
    #     # Function to return 50% CE at all Dp
    #     eff50 = lambda dp: dp*0 + 0.5

    #     m100 = opcsim.models.OPC(n_bins=3)
    #     m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

    #     urban = opcsim.load_distribution("Urban")

    #     # Test with 100% efficiency
    #     n_100 = m100.evaluate(urban, weight='number', method='simple', base='log10')

    #     # Test with 50% efficiency
    #     n_50 = m50.evaluate(urban, weight='number', method='simple', base='log10')

    #     self.assertTrue(len(n_100) == 3)

    #     for pair in zip(n_50*2, n_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    #     # Test with 100% efficiency
    #     v_100 = m100.evaluate(urban, weight='volume', method='simple', base='log10')

    #     # Test with 50% efficiency
    #     v_50 = m50.evaluate(urban, weight='volume', method='simple', base='log10')

    #     self.assertTrue(len(v_100) == 3)

    #     for pair in zip(v_50*2, v_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    # def test_efficiency_subint_number_log10(self):
    #     # Function to return 50% CE at all Dp
    #     eff50 = lambda dp: dp*0 + 0.5

    #     m100 = opcsim.models.OPC(n_bins=3)
    #     m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

    #     urban = opcsim.load_distribution("Urban")

    #     # Test with 100% efficiency
    #     n_100 = m100.evaluate(urban, weight='number', method='subint', base='log10')

    #     # Test with 50% efficiency
    #     n_50 = m50.evaluate(urban, weight='number', method='subint', base='log10')

    #     self.assertTrue(len(n_100) == 3)

    #     # Make sure that the 50% eff * 2 = 100% efficiency
    #     for pair in zip(n_50*2, n_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    # def test_efficiency_subint_number_log(self):
    #     # Function to return 50% CE at all Dp
    #     eff50 = lambda dp: dp*0 + 0.5

    #     m100 = opcsim.models.OPC(n_bins=3)
    #     m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

    #     urban = opcsim.load_distribution("Urban")

    #     # Test with 100% efficiency
    #     n_100 = m100.evaluate(urban, weight='number', method='subint', base='log')

    #     # Test with 50% efficiency
    #     n_50 = m50.evaluate(urban, weight='number', method='subint', base='log')

    #     self.assertTrue(len(n_100) == 3)

    #     # Make sure that the 50% eff * 2 = 100% efficiency
    #     for pair in zip(n_50*2, n_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    # def test_efficiency_subint_number_none(self):
    #     # Function to return 50% CE at all Dp
    #     eff50 = lambda dp: dp*0 + 0.5

    #     m100 = opcsim.models.OPC(n_bins=3)
    #     m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

    #     urban = opcsim.load_distribution("Urban")

    #     # Test with 100% efficiency
    #     n_100 = m100.evaluate(urban, weight='number', method='subint', base='none')

    #     # Test with 50% efficiency
    #     n_50 = m50.evaluate(urban, weight='number', method='subint', base='none')

    #     self.assertTrue(len(n_100) == 3)

    #     # Make sure that the 50% eff * 2 = 100% efficiency
    #     for pair in zip(n_50*2, n_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    # def test_efficiency_subint_volume(self):
    #     urban = opcsim.load_distribution("Urban")

    #     eff50 = lambda dp: dp*0 + 0.5

    #     m100 = opcsim.models.OPC(n_bins=3)
    #     m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

    #     # Test with 100% efficiency
    #     v_100 = m100.evaluate(urban, weight='volume', method='subint', base='log10')

    #     # Test with 50% efficiency
    #     v_50 = m50.evaluate(urban, weight='volume', method='subint', base='log10')

    #     self.assertTrue(len(v_100) == 3)

    #     for pair in zip(v_50*2, v_100):
    #         self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    # def test_false_model_method(self):
    #     model = opcsim.models.OPC(n_bins=3)

    #     urban = opcsim.load_distribution("Urban")

    #     with self.assertRaises(Exception):
    #         model.evaluate(urban, method='fake')

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

        # test raise for no bins set
        with self.assertRaises(Exception):
            opc = opcsim.OPC(wl=wl)

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
