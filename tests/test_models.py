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
        n_bins = 2
        dmin = 0.5
        dmax = 2.5

        opc = opcsim.models.OPC(n_bins=n_bins, dmin=dmin, dmax=dmax)

        self.assertEqual(dmin, opc.dmin)
        self.assertEqual(dmax, opc.dmax)
        self.assertEqual(n_bins, opc.n_bins)

    def test_opc_with_bins(self):
        bins = np.array([[0.5, 2.5], [2.5, 10.]])

        with self.assertRaises(Exception):
            opc = opcsim.models.OPC(bins=bins)

        bins = opcsim.utils.midpoints(bins)

        opc = opcsim.models.OPC(bins=bins)

        self.assertEqual(opc.n_bins, bins.shape[0])

        self.assertIsNotNone(opc.dlogdp)

    def test_opc_evaluate(self):
        opc = opcsim.models.OPC()
        d = opcsim.load_distribution("Urban")

        # Test dNdDp
        dNdDp = opc.evaluate(d, weight='number', base=None)
        dNdlnDp = opc.evaluate(d, weight='number', base='log')
        dNdlogDp = opc.evaluate(d, weight='number', base='log10')

        self.assertEqual(dNdDp.shape[0], opc.n_bins)
        self.assertEqual(dNdlnDp.shape[0], opc.n_bins)
        self.assertEqual(dNdlogDp.shape[0], opc.n_bins)

        # Test bad distribution
        with self.assertRaises(Exception):
            opc.evaluate(1)

    def test_histogram(self):
        opc = opcsim.models.OPC()
        d = opcsim.load_distribution("Urban")

        bb, h, bw = opc.histogram(d)

        self.assertEqual(len(bb), len(h))
        self.assertEqual(len(bb), len(bw))

    def test_opc_number(self):
        opc = opcsim.models.OPC()
        d = opcsim.load_distribution("Urban")

        vm = opc.number(d)
        va = opc.number(d, measured=False)

        # They should be equal because we have 100% efficiency
        self.assertEqual(vm, va)
        self.assertEqual(len(vm), opc.n_bins)

    def test_opc_sa(self):
        opc = opcsim.models.OPC()
        d = opcsim.load_distribution("Urban")

        vm = opc.surface_area(d)
        va = opc.surface_area(d, measured=False)

        self.assertNotEqual(vm, va)
        self.assertEqual(len(vm), opc.n_bins)

    def test_opc_volume(self):
        opc = opcsim.models.OPC()
        d = opcsim.load_distribution("Urban")

        vm = opc.volume(d)
        va = opc.volume(d, measured=False)

        self.assertNotEqual(vm, va)
        self.assertEqual(len(vm), opc.n_bins)

    def test_counting_efficiency_simple(self):
        # Function to return 50% CE at all Dp
        eff50 = lambda dp: dp*0 + 0.5

        m100 = opcsim.models.OPC(n_bins=3)
        m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

        urban = opcsim.load_distribution("Urban")

        # Test with 100% efficiency
        n_100 = m100.evaluate(
                        urban,
                        weight='number',
                        method='simple',
                        base='log10')

        # Test with 50% efficiency
        n_50 = m50.evaluate(
                        urban,
                        weight='number',
                        method='simple',
                        base='log10')

        self.assertTrue(len(n_100) == 3)

        for pair in zip(n_50*2, n_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

        # Test with 100% efficiency
        v_100 = m100.evaluate(
                        urban,
                        weight='volume',
                        method='simple',
                        base='log10')

        # Test with 50% efficiency
        v_50 = m50.evaluate(
                        urban,
                        weight='volume',
                        method='simple',
                        base='log10')

        self.assertTrue(len(v_100) == 3)

        for pair in zip(v_50*2, v_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    def test_efficiency_subint(self):
        # Function to return 50% CE at all Dp
        eff50 = lambda dp: dp*0 + 0.5

        m100 = opcsim.models.OPC(n_bins=3)
        m50 = opcsim.models.OPC(n_bins=3, ce=eff50)

        urban = opcsim.load_distribution("Urban")

        # Test with 100% efficiency
        n_100 = m100.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='log10')

        # Test with 50% efficiency
        n_50 = m50.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='log10')

        self.assertTrue(len(n_100) == 3)

        for pair in zip(n_50*2, n_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

        # Test with 100% efficiency
        v_100 = m100.evaluate(
                        urban,
                        weight='volume',
                        method='subint',
                        base='log10')

        # Test with 50% efficiency
        v_50 = m50.evaluate(
                        urban,
                        weight='volume',
                        method='subint',
                        base='log10')

        self.assertTrue(len(n_100) == 3)

        for pair in zip(v_50*2, v_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

        # Test with 100% efficiency
        n_100 = m100.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='log')

        # Test with 50% efficiency
        n_50 = m50.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='log')

        self.assertTrue(len(n_100) == 3)

        for pair in zip(n_50*2, n_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

        # Test with 100% efficiency
        n_100 = m100.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='none')

        # Test with 50% efficiency
        n_50 = m50.evaluate(
                        urban,
                        weight='number',
                        method='subint',
                        base='none')

        self.assertTrue(len(n_100) == 3)

        for pair in zip(n_50*2, n_100):
            self.assertEqual(round(pair[0], 2), round(pair[1], 2))

    def test_false_model_method(self):
        model = opcsim.models.OPC(n_bins=3)

        urban = opcsim.load_distribution("Urban")

        with self.assertRaises(Exception):
            model.evaluate(urban, method='fake')
