
.. _model_tutorial:


Using OPCSIM to Build and Model an Optical Particle Counter (OPC)
=================================================================

The following tutorial will show you how a low-cost OPC is represented
when using the opcsim software. You will learn how to build a model OPC,
how to mimic a calibration for specific aerosols, and how to evaluate
the OPC against a simulated aerosol distribution. Visualization tools
will also be discussed.

First, we import the python libraries we need and set the styles used
for plotting throughout this tutorial.

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
    
    print (opcsim.__version__)


.. parsed-literal::

    0.3.0+8.g007204b.dirty


The OPC Model
=============

The ``opcsim.OPC`` class provides a simple way to model most of the
functionality of low-cost optical particle counters. This model is based
on just a few instrument parameters which should be defined in the
manufacturer’s speficication sheet.

An OPC is defined by just the laser wavelength (``wl``), the exact
``bins`` the OPC uses as its output, and the angles for which the
scattered light is collected between (``theta``). Other work has
considered much more detail that this, including the electrical
properties and characteristics of the photodetector, etc. However, this
typically requires knowledge or information that is not commonly
available in the datasheet of a low-cost particle counter, and thus we
try to provide a method that does not rely on this information.

Within the software itself, there is flexibility in how exactly you
define the bins. In the end, you will end up with a 3xn array of values,
where n is the number of bins. Each bin is then defined by its left bin
boundary, midpoint, and right bin boundary.

To simulate an OPC using the ``opcsim.OPC`` class, we initiate as
follows:

.. code:: ipython3

    opc = opcsim.OPC(wl=0.658, n_bins=5, theta=(32., 88.))
    
    opc




.. parsed-literal::

    <class 'opcsim.models.OPC'>



When initiated with no arguments, the default arguments are used and can
be looked up in the API documentation. This sets
:math:`dmin=0.5\;\mu m`, :math:`dmax=2.5\;\mu m`, and
:math:`n_{bins}=5`. We can view the number of bins using the
``OPC.n_bins`` attribute.

.. code:: ipython3

    opc.n_bins




.. parsed-literal::

    5



We can also view the bin boundaries and midpoint diameters using the
``OPC.bins`` attribute. Here, we receieve a **3xn** array where the
first entry is the left bin boundary, the middle is the midpoint
diameter, and the last entry is the right bin boundary.

.. code:: ipython3

    opc.bins




.. parsed-literal::

    array([[0.5       , 0.58730947, 0.68986483],
           [0.68986483, 0.8103283 , 0.95182697],
           [0.95182697, 1.11803399, 1.3132639 ],
           [1.3132639 , 1.54258466, 1.81194916],
           [1.81194916, 2.12834981, 2.5       ]])



Building more specific OPC’s
----------------------------

We can build more complex - or specified - models by increasing the
number of bins in a couple of ways: (1) we can change the minimum or
maximum cutoffs, or the total number of bins:

.. code:: ipython3

    opc_10bins = opcsim.OPC(wl=0.658, n_bins=10, dmin=0.38, dmax=17.5)
    
    opc_10bins.bins




.. parsed-literal::

    array([[ 0.38      ,  0.46019969,  0.55732566],
           [ 0.55732566,  0.67495025,  0.81739972],
           [ 0.81739972,  0.98991341,  1.19883643],
           [ 1.19883643,  1.45185303,  1.75826924],
           [ 1.75826924,  2.12935514,  2.57875939],
           [ 2.57875939,  3.12301123,  3.78212839],
           [ 3.78212839,  4.58035343,  5.54704531],
           [ 5.54704531,  6.71775925,  8.13555449],
           [ 8.13555449,  9.85257798, 11.93198238],
           [11.93198238, 14.45024885, 17.5       ]])



If we are trying to model a specific OPC that has pre-defined bins, we
can also do that with the help of some utility methods. The bins
argument in the OPC class requires a **3xn** array as seen above. Often,
you may only have the bin boundary information and not the midpoints. We
can define the bins as an array of the bin boundaries, where there are 1
more entries than total number of bins, as follows:

.. code:: ipython3

    opc = opcsim.OPC(wl=0.658, bins=[0.38, 0.54, 0.78, 1.05, 1.5, 2.5])
    
    opc.bins




.. parsed-literal::

    array([[0.38      , 0.45299007, 0.54      ],
           [0.54      , 0.64899923, 0.78      ],
           [0.78      , 0.90498619, 1.05      ],
           [1.05      , 1.25499004, 1.5       ],
           [1.5       , 1.93649167, 2.5       ]])



This more or less covers how we can build various OPC’s - if you are
still unsure of what to do in a specific case, feel free to post
questions on the GitHub repository under the ‘Issues’ tab.

OPC Calibration
===============

Probably the most important part of our representation is how we
calibrate our sensor. Unlike nephelometers of photometers, OPC’s scatter
the light off of, and count, each and every particle (unless there are
too many and coincidence becomes an issue, but that’s a whole different
problem). The light scattered off a particle is proportional to the
scattering cross-section of the particle in question. To represent this
in our model, we perform Mie Theory calculations for each and every
particle in the aerosol distribution, calculating the scattering
cross-section (:math:`C_{scat}`) assuming the particle is spherical and
has some refractive index as defined when you built the distribution.

To determine which particle size bin a certain :math:`C_{scat}` value
belongs to, we must calibrate the sensor. There are many ways to do
this, and they are documented throughout the literature, especially for
expensive, nice OPC’s. Our goal is to come up with a functional
relationship that returns an OPC particle size bin for any given
:math:`C_{scat}` value. We provide two simple methods as built-in
options within the software, but users can elect to define their own.

First, let’s take a look at what a simple scattering pattern looks like
for an OPC with a wavelength of 658 nm and a collection angle of 32-88
degrees, assuming we are using PSL’s.

.. code:: ipython3

    opc = opcsim.OPC(wl=0.658, theta=(32.0, 88.0), n_bins=10)
    
    opc




.. parsed-literal::

    <class 'opcsim.models.OPC'>



.. code:: ipython3

    # generate an array of particle diameters
    dp = np.logspace(-1, 1.25, 250)
    vals = [opcsim.mie.cscat(x, wl=0.658, refr=complex(1.59, 0), theta1=32., theta2=88.) for x in dp]
    
    fig, ax = plt.subplots(1, figsize=(6, 6))
    
    ax.plot(dp, vals, lw=3)
    
    ax.semilogx()
    ax.semilogy()
    
    ax.set_title("Calibration Curve for PSLs", y=1.02)
    ax.set_xlabel("$D_p\;[\mu m]$")
    ax.set_ylabel("$C_{scat}$")
    
    sns.despine(offset=5)



.. image:: model_files/model_16_0.png


As we can see in the figure above, there are :math:`C_{scat}` values
that are not monotonically increasing as the particle diameter
increases, which may make it difficult to assign that value to 1 bin
specifically. This is a very well known and documented issue with OPC’s
and should play a role in choosing bin boundaries. Additionally, it only
gets more complicated as the particle optical properties change!

Now, we will examine the two simple approaches to modeling our ways out
of this.

Method 1: ``smooth``
--------------------

If we were to plot the bin boundaries for a given OPC on the figure
above, we would likely see periods where the :math:`C_{scat}` value
decreases from one bin boundary to the next, which makes our life
difficult. Thus, we can smooth it out by simply interpolating across any
of these ‘down’ periods. Note that this can make it quite difficult for
particles to get assigned to certain bins, especially if the
:math:`C_{scat}` values plateau - this can be somewhat solved by
combining these bins into one bin to reduce the uncertainty in bin
assignment, though it will make the correct particle sizing in that bin
less precise as well.

We can examine how this works by showing a quick example. Let’s grab the
:math:`C_{scat}` values for the OPC above at each bin boundary:

.. code:: ipython3

    opc.bin_boundaries




.. parsed-literal::

    array([0.5       , 0.58730947, 0.68986483, 0.8103283 , 0.95182697,
           1.11803399, 1.3132639 , 1.54258466, 1.81194916, 2.12834981,
           2.5       ])



.. code:: ipython3

    v = np.array([opcsim.mie.cscat(x, wl=0.658, theta1=32.,
                                          theta2=88., refr=complex(1.59, 0)) for x in opc.bin_boundaries])
    
    # print out the cscat values for each bin boundary
    v




.. parsed-literal::

    array([2.86492750e-09, 4.48206290e-09, 5.25272279e-09, 6.02393845e-09,
           6.64256044e-09, 9.09772227e-09, 1.25785156e-08, 1.25676003e-08,
           1.40492264e-08, 1.99683184e-08, 4.17442631e-08])



Let’s check to see if there are any points where the value actually
decreases as the particle diameter increases.

.. code:: ipython3

    np.diff(v) < 0




.. parsed-literal::

    array([False, False, False, False, False, False,  True, False, False,
           False])



Yup! We got one. So, let’s apply the ``opcsim.utils.squash_dips``
function to remove them:

.. code:: ipython3

    opcsim.utils.squash_dips(v)




.. parsed-literal::

    array([2.86492750e-09, 4.48206290e-09, 5.25272279e-09, 6.02393845e-09,
           6.64256044e-09, 9.09772227e-09, 1.08326613e-08, 1.25676003e-08,
           1.40492264e-08, 1.99683184e-08, 4.17442631e-08])



We can see above that the offending value has now been modified.

Method 2: ``linear``
--------------------

The second approach to generating our map of :math:`C_{scat}` to
:math:`D_p` values is to fit a line to the data in log-log space. Here,
we begin by plotting the bin boundaries on top of the PSL scattering
line for an OPC that models the Alphasense OPC-N2.

.. code:: ipython3

    from scipy.interpolate import interp1d
    
    fig, ax = plt.subplots(1, figsize=(8, 8))
    
    ax.plot(dp, vals, lw=3)
    
    ax.set_xlabel("$D_p$, $\mu m$")
    ax.set_ylabel("$C_{scat}$, $cm^2$/particle")
    
    # draw some bins
    f = interp1d(x=dp, y=vals)
    
    bb = np.array([0.38, 0.54, 0.78, 1.05, 1.34, 1.59, 2.07, 3., 4., 
                   5., 6.5, 8., 10., 12., 14., 16., 17.5])
    
    ax.scatter(bb, f(bb), color='r', s=50, label="OPC-N2 Bin Boundaries")
    ax.legend()
    
    ax.semilogx()
    ax.semilogy()
    
    sns.despine(offset=5)



.. image:: model_files/model_27_0.png


We can fit a line to these data points, which will generate an easy
mapping between the two values - however, it will likely undersize
particles in some regions and oversize them in others. You can easily
create complicated polynomial or machine-learning approaches to matching
these values; however, it will all get thrown into chaos when we begin
to incorporate particles of various optical properties and morphologies.

So, how do we calibrate our sensor?
-----------------------------------

It is recommended you begin with one of the two approaches above, which
are both built into the ``OPC.calibrate`` function. Once you get a hang
of those, you can go ahead and get more sophisticated.

To calibrate our sensor, we simply define the material we want to use to
calibrate. We can define the material in one of two ways:

1. you can input a string with the material name
2. you can input the materials complex refractive index at the
   wavelength of your device

The material lookup table is quite limited and are mostly values taken
as close to 658 nm as possible. If you would like to add to this table,
please open an issue on GitHub.

Current Options
~~~~~~~~~~~~~~~

==================== ========================= ================
option               Material                  Refractive Index
==================== ========================= ================
``psl``              Polystyrene Latex Spheres 1.592 + 0i
``ammonium_sulfate`` Ammonium Sulfate          1.521 + 0i
``sodium_chloride``  Sodium Chloride           1.5405 + 0i
``sodium_nitrate``   Sodium Nitrate            1.448 + 0i
``black_carbon``     Black Carbon              1.95 + 0.79i
``sulfuric_acid``    Sulfuric Acid             1.427 + 0i
``soa``              Secondary Organic Aerosol 1.4 + 0.002i
``h20``              Water                     1.333 + 0i
==================== ========================= ================

To use the method, we simply state the material and the method. We can
also get more complicated (check the API docs) if you want to send
custom arguments to the fitting function, etc. A simple example would
be:

.. code:: ipython3

    opc.calibrate("psl", method="smooth")
    
    opc.calibration_function




.. parsed-literal::

    functools.partial(<bound method OPC._digitize_opc_bins of <class 'opcsim.models.OPC'>>, cscat_boundaries=array([2.86492750e-09, 4.48206290e-09, 5.25272279e-09, 6.02393845e-09,
           6.64256044e-09, 9.09772227e-09, 1.08326613e-08, 1.25676003e-08,
           1.40492264e-08, 1.99683184e-08, 4.17442631e-08]))



Alternatively, we could calibrate using the linear fit:

.. code:: ipython3

    opc.calibrate("psl", method="linear")
    
    opc.calibration_function




.. parsed-literal::

    functools.partial(<bound method OPC._digitize_opc_bins of <class 'opcsim.models.OPC'>>, cscat_boundaries=array([3.12354813e-09, 3.85694293e-09, 4.76253547e-09, 5.88075700e-09,
           7.26153183e-09, 8.96650627e-09, 1.10718009e-08, 1.36714091e-08,
           1.68813933e-08, 2.08450672e-08, 2.57393936e-08]))



When the calibration is executed, it saves the calibration function,
which is a digitizer which simply uses the :math:`C_{scat}` values we
computed at each bin boundary during the calibration to return the OPC
bin for any given :math:`C_{scat}` value. If a :math:`C_{scat}` value is
either too high or too low (above or below the boundaries we computed),
then it will not assign that particle to a bin."

Evaluating an OPC for a given Aerosol Distribution
==================================================

The entire point of this software is to model how various OPC’s and
other particle sensors “see” aerosols in the real world. The Aerosol
Distribution tutorial showed how we can model realistic aerosol
distributions. This section will review how the OPC’s we just built
“see” these distributions. Generally speaking, the process works as
follows:

1. For every particle in the aerosol distribution, we calculate the
   scattering cross-section
2. For each of the scattering cross-section’s we just calculated, we
   assign them to a bin of the OPC based on the calibration curve we
   generated
3. We end up with a histogram in the form of the cumulative number of
   particles in each OPC bin, just as we would if we were using a
   commercially available OPC

The base method we use to obtain these numbers is ``OPC.evaluate``. This
method returns the number of particles the OPC “sees” in each bin for a
given aerosol distribution. One component that hasn’t yet been discussed
is how the OPC reacts as relative humidity changes. Each of these
methods takes an argument for relative humidity (``rh``) which will
automatically calculate the change in particle size due to hygroscopic
growth per :math:`\kappa`-kohler theory; the refractive index and
density will also change accordingly based on the amount of growth that
takes place.

To evaluate a distribution, we need to first define the
AerosolDistribution and the OPC.

-  Show examples for a few different scenarios, including growth due to
   RH
-  show examples for different calibration approaches

.. code:: ipython3

    # build a single mode of ammonium sulfate
    d = opcsim.AerosolDistribution("Amm. Sulfate")
    
    d.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0.53, rho=1.77, refr=complex(1.521, 0))
    
    d




.. parsed-literal::

    AerosolDistribution: Amm. Sulfate



.. code:: ipython3

    # build an OPC
    opc = opcsim.OPC(wl=0.658, n_bins=24, dmin=0.35, dmax=40.)
    
    opc




.. parsed-literal::

    <class 'opcsim.models.OPC'>



Let’s go ahead and calibrate the OPC using ammonium sulfate:

.. code:: ipython3

    opc.calibrate("ammonium_sulfate", method="smooth")
    
    # ax = opcsim.plots.calibration_plot(opc)

Now, let’s go ahead and evaluate the OPC for the previously defined
distribution. This will return the number of particles in each bin:

.. code:: ipython3

    vals = opc.evaluate(d, rh=0.0)
    
    vals




.. parsed-literal::

    array([111.,  71.,  34.,  22.,   7.,   4.,   0.,   0.,   0.,   0.,   0.,
             0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
             0.,   0.])



Interating across the Particle Size Distribution
------------------------------------------------

Another important task for evaluation OPC’s is comparing the integrated
values across some pre-defined particle size range. We can accomplish
this by using the ``OPC.integrate`` method to calculate either the total
number of particles, total surface area, total volume, or total mass
betweeen any two diameters.

Example:

.. code:: ipython3

    opc.integrate(d, dmin=0., dmax=2.5, weight='number', rh=0.)




.. parsed-literal::

    249.0



.. code:: ipython3

    opc.integrate(d, dmin=0., dmax=2.5, weight='number', rh=85.)




.. parsed-literal::

    562.0



Visualizing the OPC Histogram
-----------------------------

Last, it is fairly useful and important for us to be able to easily
visualize our results. There are functions built-in to compute the
histogram to use with ``opcsim.plots.histplot`` to easily plot the
particle size distribution.

Let’s go ahead and visualize our results from earlier:

.. code:: ipython3

    # compute the histogram
    lb, data1, ddp = opc.histogram(d, weight="number", base="log10", rh=0.0)
    lb, data2, ddp = opc.histogram(d, weight="number", base="log10", rh=95.0)

.. code:: ipython3

    fig, ax = plt.subplots(1, figsize=(8, 6))
    
    ax = opcsim.plots.pdfplot(d, ax=ax, weight="number")
    ax = opcsim.plots.histplot(data=data1, bins=opc.bins, ax=ax, label="RH=0%", plot_kws=dict(linewidth=2, fill=True, alpha=.5))
    ax = opcsim.plots.histplot(data=data2, bins=opc.bins, ax=ax, label="RH=95%", plot_kws=dict(linewidth=2))
    
    ax.set_ylim(0, None)
    ax.set_xlim(0.01, None)
    ax.legend()
    
    sns.despine()



.. image:: model_files/model_46_0.png


Effects of Aerosol Optical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the key (and cool) features of OPCSIM is it allows us to explore
how an OPC would react to aerosols of different optical properties,
namely changes in the refractive index.

Here, we calibrate our OPC to PSL’s (commonly done in the lab) and then
try to evaluate our OPC on an Urban distribution.

First, lets build our OPC and calibrate it with PSL’s
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    opc = opcsim.OPC(
        wl=658e-3, 
        bins=[0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0, 2.5, 5., 7.5, 10.], 
        theta=(32., 88.))
    
    opc.calibrate(material="psl", method="smooth")
    
    opc.calibration_function




.. parsed-literal::

    functools.partial(<bound method OPC._digitize_opc_bins of <class 'opcsim.models.OPC'>>, cscat_boundaries=array([3.91693267e-10, 9.62064230e-10, 1.73340831e-09, 2.23944162e-09,
           2.86492750e-09, 3.94881459e-09, 4.56679048e-09, 5.38850764e-09,
           5.90727633e-09, 6.49417520e-09, 8.75476317e-09, 9.97510289e-09,
           1.11954426e-08, 2.27235414e-08, 4.17442631e-08, 1.10415287e-07,
           2.54973417e-07, 4.34401239e-07]))



.. _next-lets-build-our-distribution-of-urban-aerosols-using-ri-values-from-ebert-et-al-20041:

Next, let’s build our Distribution of Urban Aerosols using RI values from `Ebert, et al (2004) <https://www.sciencedirect.com/science/article/pii/S1352231004008027>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    d1 = opcsim.AerosolDistribution("PSL")
    d2 = opcsim.AerosolDistribution("Urban Upper Bound Estimate")
    
    d1.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0., refr=complex(1.591, 0))
    # d1.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0.2, refr=complex(1.6, 0.034))
    d2.add_mode(n=1e3, gm=250e-3, gsd=1.65, kappa=0.2, refr=complex(1.73, 0.086))

Now, let’s compute the histogram for a few different scenarios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    lb, data1, ddp = opc.histogram(d1, weight="number", base="log10", rh=0.0)
    lb, data2, ddp = opc.histogram(d2, weight="number", base="log10", rh=.0)

Finally, let’s plot the results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    fig, ax = plt.subplots(1, figsize=(8, 6))
    
    ax = opcsim.plots.pdfplot(d1, ax=ax, weight="number", label="Actual Distribution")
    
    # ax.plot(opc.midpoints, data1)
    ax = opcsim.plots.histplot(data=data1, bins=opc.bins, ax=ax, label="PSL", plot_kws=dict(linewidth=0, fill=True, alpha=.35))
    ax = opcsim.plots.histplot(data=data2, bins=opc.bins, ax=ax, label="Urban Aerosol", plot_kws=dict(linewidth=3))
    
    ax.set_ylim(0, None)
    ax.set_xlim(0.01, None)
    ax.legend(bbox_to_anchor=(.85, 1))
    
    sns.despine()



.. image:: model_files/model_54_0.png


.. code:: ipython3

    opc.evaluate(d2, rh=0)




.. parsed-literal::

    array([122., 123.,  37.,  59.,  59.,   0.,   0.,   0.,   0.,   0.,   0.,
             0.,   0.,   0.,   0.,   0.,   0.])



