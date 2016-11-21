import unittest
import opcsim
import pandas as pd
import os

from opcsim.equations import *

_mode_ = (1000, 0.05, 1.25)
_mode_label = 'ex_mode'

class SetupTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_dn_ddp(self):
        val = round(opcsim.equations.pdf.dn_ddp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 5.743)

        # Make sure it works on an array
        vals = opcsim.equations.pdf.dn_ddp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 5.743)
        self.assertEqual(round(vals[1], 3), 357.566)

    def test_dn_dlndp(self):
        val = round(opcsim.equations.pdf.dn_dlndp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 1.436)

        # Test on an array
        vals = opcsim.equations.pdf.dn_dlndp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 1.436)
        self.assertEqual(round(vals[1], 3), 178.783)

    def test_dn_dlogdp(self):
        val = round(opcsim.equations.pdf.dn_dlogdp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 3.306)

        vals = opcsim.equations.pdf.dn_dlogdp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 3.306)
        self.assertEqual(round(vals[1], 3), 411.663)

    def test_ds_ddp(self):
        val     = round(opcsim.equations.pdf.ds_ddp(0.25, 100, 0.5, 1.25), 3)

        check   = round(opcsim.equations.pdf.dn_ddp(0.25, 100, 0.5, 1.25) * np.pi * (.25 ** 2), 3)

        self.assertEqual(val, check)

        # Make sure it works for arrays as well
        vals    = opcsim.equations.pdf.ds_ddp(np.array([0.25, 0.5]), 100, 0.5, 1.25)
        checks  = opcsim.equations.pdf.dn_ddp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), round(checks[0] * np.pi * 0.25 ** 2 , 3))
        self.assertEqual(round(vals[1], 3), round(checks[1] * np.pi * 0.5 ** 2 , 3))

    def test_ds_dlndp(self):
        val     = round(opcsim.equations.pdf.ds_dlndp(0.25, 100, 0.5, 1.25), 3)

        check   = round(opcsim.equations.pdf.dn_dlndp(0.25, 100, 0.5, 1.25) * np.pi * (.25 ** 2), 3)

        self.assertEqual(val, check)

        # Make sure it works for arrays as well
        vals    = opcsim.equations.pdf.ds_dlndp(np.array([0.25, 0.5]), 100, 0.5, 1.25)
        checks  = opcsim.equations.pdf.dn_dlndp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), round(checks[0] * np.pi * 0.25 ** 2 , 3))
        self.assertEqual(round(vals[1], 3), round(checks[1] * np.pi * 0.5 ** 2 , 3))

    def test_ds_dlogdp(self):
        val     = round(opcsim.equations.pdf.ds_dlogdp(0.25, 100, 0.5, 1.25), 3)

        check   = round(opcsim.equations.pdf.ds_ddp(0.25, 100, 0.5, 1.25) * np.log(10) * .25, 3)

        self.assertEqual(val, check)

        # Make sure it works for arrays as well
        vals    = opcsim.equations.pdf.ds_dlogdp(np.array([0.25, 0.5]), 100, 0.5, 1.25)
        checks  = opcsim.equations.pdf.ds_ddp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), round(checks[0] * np.log(10) * 0.25, 3))
        self.assertEqual(round(vals[1], 3), round(checks[1] * np.log(10) * 0.5, 3))

    def test_dv_ddp(self):
        val     = round(opcsim.equations.pdf.dv_ddp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 0.047)

        vals    = opcsim.equations.pdf.dv_ddp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 0.047)
        self.assertEqual(round(vals[1], 3), 23.403)

    def test_dv_dlndp(self):
        val     = round(opcsim.equations.pdf.dv_dlndp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 0.012)

        vals    = opcsim.equations.pdf.dv_dlndp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 0.012)
        self.assertEqual(round(vals[1], 3), 11.701)

    def test_dv_dlogdp(self):
        val     = round(opcsim.equations.pdf.dv_dlogdp(0.25, 100, 0.5, 1.25), 3)

        self.assertEqual(val, 0.027)

        vals    = opcsim.equations.pdf.dv_dlogdp(np.array([0.25, 0.5]), 100, 0.5, 1.25)

        self.assertEqual(round(vals[0], 3), 0.027)
        self.assertEqual(round(vals[1], 3), 26.943)

    def test_cdf_n(self):
        N   = 1000.0
        val = opcsim.equations.cdf.Nt(10, N, 0.4, 1.25)

        # Make sure integrating all the way out returns the total number of
        # particles
        self.assertEqual(N, val)

    def test_cdf_s(self):
        N   = 1000.0
        Dp  = 0.25
        S   = N * np.pi * Dp ** 2

        val = opcsim.equations.cdf.St(10, N, Dp, 1.25)

        self.assertGreaterEqual(val, S)

    def test_cdf_v(self):
        N   = 1000.0
        Dp  = 0.25
        V   = N * np.pi / 6. * Dp ** 3

        val = opcsim.equations.cdf.Vt(10, N, Dp, 1.25)

        self.assertGreaterEqual(val, V)
