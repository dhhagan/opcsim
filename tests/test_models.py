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

        self.assertNotEqual(vm, va)
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
