import unittest
import opcsim
import pandas as pd
import os
import matplotlib as mpl

from opcsim.distributions import *
from opcsim.models import *
from opcsim.plots import *

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_histplot(self):
        opc = opcsim.OPC()
        d = opcsim.load_distribution("Urban")

        ax = opcsim.plots.histplot(opc.evaluate(d), opc.bins)

        self.assertIsNotNone(ax)
