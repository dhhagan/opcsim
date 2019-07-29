
.. _nephelometer_tutorial:


Using OPCSIM to Simulate a Nephelometer
=======================================

This section of the tutorial will walk you through how we model
Nephelometers, how you can build/model a Nephelometer, and how we can
evaluate Nephelometers across a wide range of conditions using this
tool.

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

Nephelometer Representation
---------------------------

In OPCSIM, we define a Nephelometer using two parameters: the wavelength
of light used in the device and its viewing angle. Unlike photometers
and some optical particle counters, most low-cost commercial
nephelometers gather light across as wide a range of angles as possible.
This minimizes some of the uncertainty associated with the Mie resonance
and allows manufacturers to use cheap photo-detectors while still
gathering enough signal to distinguish from noise.

To build a Nephelometer, simply initialize using the
``opcsim.Nephelometer`` class:

.. code:: ipython3

    # init a nephelometer with a 658 nm laser, gathering light from between 7-173 degrees
    neph = opcsim.Nephelometer(wl=0.658, theta=(7., 173))
    
    neph




.. parsed-literal::

    <opcsim.models.Nephelometer at 0x110594160>



Calibration
-----------

Nephelometers gather the total scattered light from many anglees across
an entire aerosol distribution. Typically, users of low-cost
nephelometers co-locate their device with a reference device of higher
(or known) quality and simply compare the output signal from the
nephelometer to the integrated mass value (i.e. :math:`PM_1`,
:math:`PM_{2.5}`, or :math:`PM_{10}`) from the reference device. To keep
things as simple and realistic as possible, we follow this approach.

To calibrate a nephelometer in OPCSIM, you provide an aerosol
distribution to the ``calibrate`` method - the actual mass values for
:math:`PM_1`, :math:`PM_{2.5}`, and :math:`PM_{10}` are calculated
exactly and the total scattered light is computed as well. The ratio
between the total scattered light and each of the mass loadings are
stored as calibration factors and are used again when evaluating
previously unseen distributions.

To calibrate our nephelometer above to a synthetic distribution of
Ammonium Sulfate:

.. code:: ipython3

    d1 = opcsim.AerosolDistribution("AmmSulf")
    
    d1.add_mode(n=1e4, gm=125e-3, gsd=1.5, refr=complex(1.521, 0), kappa=0.53, rho=1.77)
    
    # calibrate the nephelometer at 0% RH
    neph.calibrate(d1, rh=0.)

We can explore the calibration factors that were just determined - the
units are a bit arbitrary, since we don’t consider the intensity/power
of the laser as we assume it is constant. Thus, these units are
something like :math:`cm^2/(\mu g/ m^3)`

.. code:: ipython3

    neph.pm1_ratio




.. parsed-literal::

    1.1744058563022677e-08



Similarly, we get ratio’s for :math:`PM_{2.5}` and :math:`PM_{10}`:

.. code:: ipython3

    neph.pm25_ratio




.. parsed-literal::

    1.174352137965417e-08



.. code:: ipython3

    neph.pm10_ratio




.. parsed-literal::

    1.1743521375694491e-08



Evaluating a Nephelometer for New Aerosol Distributions
-------------------------------------------------------

The entire point of this tool is to be able to simulate what would
happen under different circumstances. To do so, we use the ``evaluate``
method, which takes an AerosolDistribution as an argument (as well as an
optional relative humidity) and returns the total scattered light,
:math:`PM_1`, :math:`PM_{2.5}`, and :math:`PM_{10}`.

.. code:: ipython3

    # evaluate the same distribution we used to calibrate
    neph.evaluate(d1, rh=0.)




.. parsed-literal::

    (4.4544608528839e-07, 37.92948433439533, 37.93121934108556, 37.9312193538752)



.. code:: ipython3

    # evaluate the same distribution we used to calibrate, but at a higher RH
    neph.evaluate(d1, rh=85.0)




.. parsed-literal::

    (2.0830803505720823e-06,
     177.3731235580573,
     177.38123712884374,
     177.3812371886531)



What if we went ahead and tried to evaluate on a totally unseen
distribution? Let’s go ahead and evaluate on an **urban** distribution:

.. code:: ipython3

    d2 = opcsim.load_distribution("urban")
    
    d2




.. parsed-literal::

    AerosolDistribution: urban



First, let’s determine the actual :math:`PM_1`, :math:`PM_{2.5}`, and
:math:`PM_{10}` loadings for this distribution:

.. code:: ipython3

    print ("PM1 = {:.2f} ug/m3".format(d2.cdf(dmin=0., dmax=1., weight='mass', rho=1.65)))
    print ("PM2.5 = {:.2f} ug/m3".format(d2.cdf(dmin=0., dmax=2.5, weight='mass', rho=1.65)))
    print ("PM10 = {:.2f} ug/m3".format(d2.cdf(dmin=0., dmax=10., weight='mass', rho=1.65)))


.. parsed-literal::

    PM1 = 8.97 ug/m3
    PM2.5 = 9.00 ug/m3
    PM10 = 9.00 ug/m3


Next, let’s evaluate the Nephelometer:

.. code:: ipython3

    neph.evaluate(d2, rh=0.)




.. parsed-literal::

    (1.7166105544467465e-07,
     14.61684259521375,
     14.617511212784953,
     14.617511217713684)



So, we’re off by about a factor of 2, in part due to differences in
assumed density and in part due to the fact the urban distribution
scatters less light per unit mass than our calibration aerosol.

