#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

    def test_pdfplot(self):
        opc = opcsim.OPC()
        d = opcsim.load_distribution("Urban")

        ax = opcsim.plots.pdfplot(d)

        self.assertIsNotNone(ax)

        # Try with a different weight
        ax = opcsim.plots.pdfplot(d, weight='volume')

        self.assertIsNotNone(ax)

        with self.assertRaises(Exception):
            ax = opcsim.plots.pdfplot(1)

        # Test invalid weight
        with self.assertRaises(ValueError):
            ax = opcsim.plots.pdfplot(d, weight='mass2')

        # Test with_modes
        ax = opcsim.plots.pdfplot(d, with_modes=True)
        self.assertIsNotNone(ax)

    def test_pdf_plot_with_fill(self):
        opc = opcsim.OPC()
        d = opcsim.load_distribution("Urban")

        ax = opcsim.plots.pdfplot(d, fill=True)

        self.assertIsNotNone(ax)

        # Try with a different weight
        ax = opcsim.plots.pdfplot(d, weight='volume')

    def test_cdfplot(self):
        opc = opcsim.OPC()
        d = opcsim.load_distribution("Urban")

        ax = opcsim.plots.cdfplot(d)

        self.assertIsNotNone(ax)

        # Try with a different weight
        ax = opcsim.plots.cdfplot(d, weight='volume')

        self.assertIsNotNone(ax)

        with self.assertRaises(Exception):
            ax = opcsim.plots.cdfplot(1)

        # Test invalid weight
        with self.assertRaises(ValueError):
            ax = opcsim.plots.cdfplot(d, weight='mass2')
