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

    def test_coef_pi_tau(self):
        pi, tau = opcsim.mie.coef_pi_tau(theta=30.0, x=.5)
        self.assertAlmostEqual(pi[0], 1.)
        self.assertAlmostEqual(0.86, round(tau[0], 2), places=1)

    def test_coef_ab(self):
        a, b = opcsim.mie.coef_ab(refr=complex(1.5, 0), x=.5)
        self.assertAlmostEqual(a[0].real, 6.06e-4, places=2)
        self.assertAlmostEqual(b[0].real, 7.5e-7, places=2)
    
    def test_s1s2(self):
        s1, s2 = opcsim.mie.s1s2(refr=complex(1.5, 0), x=.5, theta=30.)
        self.assertAlmostEqual(s1.real, 9.1e-4, places=2)
        self.assertAlmostEqual(s2.real, 1.8e-4, places=2)
    
    def test_cscat(self):
        rv = opcsim.mie.cscat(dp=0.5, wl=0.658, refr=complex(
            1.9, .5), theta1=32., theta2=88.)
        self.assertTrue(type(rv), float)