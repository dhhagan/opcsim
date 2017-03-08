import unittest
import opcsim
import pandas as pd
import os

from opcsim.distributions import *
from opcsim.models import *
from opcsim.viz import distplot

"""
_mode_          = (1000, 0.05, 1.25)
_mode2_         = (1000, 0.1, 1.5)

_mode_label     = 'ex_mode'
_mode_label2    = 'ex_mode2'

def fake_dist():
    dist = AerosolDistribution()

    dist.add_mode(_mode_, _mode_label)

    return dist


class SetupTestCase(unittest.TestCase):
    def setUp(self):
        self._d_    = fake_dist()

    def tearDown(self):
        pass

    def test_distplot(self):
        plt = distplot(self._d_)

        # Test bad weight
        with self.assertRaises(Exception):
            distplot(self._d_, weight = 'error')

        plt = distplot(self._d_, with_modes = True)

    def test_multi_plot_distplot(self):
        plt = distplot(self._d_, weight = ['number', 'surface'])
"""
