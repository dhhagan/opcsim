
.. _nephelometer_tutorial:


Using OPCSIM to Simulate a Nephelometer
=======================================

.. code:: ipython3

    # Make imports
    import opcsim
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticks
    import seaborn as sns
    
    %matplotlib inline
    
    # turn off warnings temporarily
    import warnings
    warnings.simplefilter('ignore')
    
    # Let's set some default seaborn settings
    sns.set(context='notebook', style='ticks', palette='dark', font_scale=1.75, 
            rc={'figure.figsize': (12,6), **opcsim.plots.rc_log})

.. code:: ipython3

    class Nephelometer(object):
        def __init__(self, wl, theta=(7., 173.), **kwargs):
            self.wl = wl
            self.theta = theta
            self.calibration_function = None
            self.pm1_ratio = None
            self.pm25_ratio = None
            self.pm10_ratio = None
        
        def calibrate(self, distribution, **kwargs):
            """Set the calibration function that relates the total scattered light
            to mass.
            """
            n_bins = kwargs.pop("n_bins", 100)
            
            # compute PM1, PM25, and PM10 for the distribution
            pm1 = distribution.cdf(dmin=0., dmax=1., weight="mass")
            pm25 = distribution.cdf(dmin=0., dmax=2.5, weight="mass")
            pm10 = distribution.cdf(dmin=0., dmax=10., weight="mass")
            
            print (pm1, pm25, pm10)
            
            # compute the total scattered light the Nephelometer "see's"
            total_cscat = self._sum_across_distribution(distribution, n_bins=n_bins, rh=0.)
                
            # set the ratios
            self.pm1_ratio = total_cscat / pm1
            self.pm25_ratio = total_cscat / pm25
            self.pm10_ratio = total_cscat / pm10
            
            return
        
        def _sum_across_distribution(self, distribution, n_bins=100, rh=0., **kwargs):
            """Return the total Cscat value when summed across an entire distribution.
            """
            total_cscat = 0.
            for m in distribution.modes:
                gm = m["GM"]
                gsd = m["GSD"]
                refr = m["refr"]
                
                # alter GM and RI per RH conditions
                gm = opcsim.utils.k_kohler(diam_dry=gm, kappa=m["kappa"], rh=rh)
                
                pct_dry = m["GM"]**3 / gm**3
                
                refr = opcsim.utils.ri_eff(species=[refr, complex(1.333, 0)], weights=[pct_dry, 1-pct_dry])
    
                # compute the range over which to make calculations - we use Dpg/sigma^4 to Dpg*sigma^4
                bounds = np.logspace(start=np.log10(gm/(gsd**4)), stop=np.log10(gm*(gsd**4)), num=n_bins)
                
                # compute the midpoints for each bin
                midpoints = np.mean([bounds[:-1], bounds[1:]], axis=0)
                
                # compute the number of particles in each bin
                n = np.array([distribution.cdf(dmin=a, dmax=b, mode=m["label"], rh=rh) for a, b in zip(bounds[:-1], bounds[1:])])
                
                # compute the mean Cscat value for each bin
                cscat = [opcsim.mie.cscat(dp, wl=self.wl, refr=refr, theta1=self.theta[0], theta2=self.theta[1]) for dp in midpoints]
                
                # take the product of the two - results in total cm2
                signal = (n*cscat).sum()
                
                # add to the running total
                total_cscat += signal
                
            return total_cscat
        
        def evaluate(self, distribution, rh=0., **kwargs):
            """
            """
            total_cscat = self._sum_across_distribution(distribution, rh=rh, **kwargs)
            
            pm1 = total_cscat / self.pm1_ratio
            pm25 = total_cscat / self.pm25_ratio
            pm10 = total_cscat / self.pm10_ratio
            
            return total_cscat, pm1, pm25, pm10
        
    d = opcsim.AerosolDistribution("Cal1")
    d.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0.53, refr=complex(1.521, 0), rho=1.77)
    
    tmp = Nephelometer(wl=0.658)
    
    tmp.calibrate(d)


.. parsed-literal::

    40.160330923901654 44.716004553849146 44.75994510431068


.. code:: ipython3

    tmp.evaluate(d)




.. parsed-literal::

    (1.5146881014163737e-06,
     40.160330923901654,
     44.716004553849146,
     44.75994510431068)



.. code:: ipython3

    d2 = opcsim.AerosolDistribution("Cal2")
    d2.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0, refr=complex(1.592, 0), rho=1.)
    
    tmp.evaluate(d2)




.. parsed-literal::

    (1.748375082925014e-06,
     46.35629067377946,
     51.61481634192903,
     51.665536066659904)



.. code:: ipython3

    d.cdf(dmin=0., dmax=1., weight="mass")




.. parsed-literal::

    40.160330923901654



