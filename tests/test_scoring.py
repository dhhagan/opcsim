#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import opcsim
import pandas as pd
import numpy as np

from opcsim.distributions import *
from opcsim.models import *
from opcsim.metrics import compute_bin_assessment

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_compute_bin_assessments(self):
        opc = OPC(wl=0.658, n_bins=5)

        opc.calibrate("psl")

        d = AerosolDistribution("test")
        d.add_mode(n=10000, gm=0.1, gsd=2, label="mode 1", kappa=0.53, refr=complex(1.55, 0))

        rv = compute_bin_assessment(opc, refr=complex(1.55, 0), kappa=0.53)

        self.assertTrue(isinstance(rv, pd.DataFrame))
