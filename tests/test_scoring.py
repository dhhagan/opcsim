#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

    def test_nv_score(self):
        eta_50 = lambda dp: dp*0 + 0.5
        eta_100 = lambda dp: dp*0 + 1

        urban = opcsim.load_distribution("Urban")

        opc50 = opcsim.models.OPC(n_bins=5, ce=eta_50)
        opc100 = opcsim.models.OPC(n_bins=5, ce=eta_100)

        nv50 = opcsim.metrics.nv_score(opc50, urban)
        nv100 = opcsim.metrics.nv_score(opc100, urban)

        self.assertNotEqual(nv50, nv100)

        # Make sure if we integrate to 1 micron, we get a higher N/V than
        # if we integrate to 2.5 microns.
        nv_1 = opcsim.metrics.nv_score(opc100, urban, dmax=1)
        nv_25 = opcsim.metrics.nv_score(opc100, urban, dmax=2.5)

        self.assertLess(nv_25, nv_1)

    def test_two_opcs(self):
        urban = opcsim.load_distribution("Urban")

        opc1 = opcsim.models.OPC(n_bins=5, dmin=1, dmax=2.5)
        opc2 = opcsim.models.OPC(n_bins=5, dmin=0.5, dmax=2.5)

        nv_1 = opcsim.metrics.nv_score(opc1, urban)
        nv_2 = opcsim.metrics.nv_score(opc2, urban)

        self.assertLess(nv_1, nv_2)

    def test_nv_score(self):
        eta_50 = lambda dp: dp*0 + 0.5
        eta_100 = lambda dp: dp*0 + 1

        urban = opcsim.load_distribution("Urban")

        opc50 = opcsim.models.OPC(n_bins=5, ce=eta_50)
        opc100 = opcsim.models.OPC(n_bins=5, ce=eta_100)

        vv50 = opcsim.metrics.vv_score(opc50, urban)
        vv100 = opcsim.metrics.vv_score(opc100, urban)

        self.assertNotEqual(vv50, vv100)
        self.assertLess(vv50, vv100)

    def test_vv_2_sensors(self):
        urban = opcsim.load_distribution("Urban")

        opc1 = opcsim.models.OPC(n_bins=2, dmin=0.5, dmax=2.5)
        opc2 = opcsim.models.OPC(n_bins=20, dmin=0.5, dmax=2.5)

        # with more bins across the same range, more bins should get V/V closer
        # to 1
        vv1 = opcsim.metrics.vv_score(opc1, urban, dmin=0.5, dmax=2.5)
        vv2 = opcsim.metrics.vv_score(opc2, urban, dmin=0.5, dmax=2.5)

        check1 = abs(1 - vv1)
        check2 = abs(1 - vv2)

        self.assertGreater(check1, check2)
