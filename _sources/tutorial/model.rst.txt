
.. _model_tutorial:


Using OPCSIM to Build a Simulated OPC
=====================================

The following tutorial will show you how we represent a low-cost optical
particle counter in the opcsim software. You will learn how to build a
model OPC and evaluate it against a simulated aerosol distribution.
Visualization tools will also be discussed.

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

The OPC Model
=============

The ``opcsim.OPC`` class provides an easy-to-use interface for
simulating low-cost optical particle counters. This simple model is
based on just a few instrument parameters. We define an OPC with the
following statements:

-  An OPC has ``n_bins`` histogram bins
-  An OPC has a lower cutoff at some diameter (``dmin``) and upper
   cutoff at some diameter (``dmax``)
-  An OPC has a counting efficiency that can either be constant or vary
   with particle diameter as defined by the function ``ce``

To simulate an OPC using the ``opcsim.OPC`` class, we initiate as
follows:

.. code:: ipython3

    opc = opcsim.OPC()

When initiated with no arguments, the default arguments are used. This
sets :math:`dmin=0.5\;\mu m`, :math:`dmax=2.5\;\mu m`, and
:math:`n_{bins}=1`. We can view the number of bins using the
``OPC.n_bins`` attribute.

.. code:: ipython3

    opc.n_bins




.. parsed-literal::

    1



We can also view the bin boundaries and midpoint diameters using the
``OPC.bins`` attribute. Here, we receieve a ``3xn`` array where the
first entry is the left bin boundary, the middle is the midpoint
diameter, and the last entry is the right bin boundary.

.. code:: ipython3

    opc.bins




.. parsed-literal::

    array([[0.5       , 1.11803399, 2.5       ]])



Building an OPC
---------------

We can build more complex models by increasing the number of bins in a
couple of ways: (1) we can change the minimum or maximum cutoffs, or the
total number of bins:

.. code:: ipython3

    opc_10bins = opcsim.OPC(n_bins=10, dmin=0.3, dmax=17)
    
    opc_10bins.bins




.. parsed-literal::

    array([[ 0.3       ,  0.36710275,  0.44921476],
           [ 0.44921476,  0.54969325,  0.67264635],
           [ 0.67264635,  0.82310108,  1.0072089 ],
           [ 1.0072089 ,  1.23249719,  1.50817703],
           [ 1.50817703,  1.84551978,  2.25831796],
           [ 2.25831796,  2.76344911,  3.38156589],
           [ 3.38156589,  4.13794047,  5.06349775],
           [ 5.06349775,  6.19607983,  7.58199316],
           [ 7.58199316,  9.2779018 , 11.35314422],
           [11.35314422, 13.89256822, 17.        ]])



If we are trying to mimic a specific OPC that has pre-defined bins, we
can also do that with the help of some utility methods. The bins
argument in the OPC class requires a ``3xn`` array as seen above. Often,
you may only have the bin boundary information and not the midpoints.
Typically, we use the logarithmic mean rather than the arithmetic mean,
though we have made both available through the ``opcsim.midpoints``
utility function.

For example, let’s calculate the bins for an OPC like the Dylos DC1100
Pro. This OPC has two bins (0.5-2.5, 2.5-10). How do we build the bins?

.. code:: ipython3

    bins = np.array([[0.5, 2.5], [2.5, 10]])
    
    bins = opcsim.midpoints(bins)
    
    bins




.. parsed-literal::

    array([[ 0.5       ,  1.11803399,  2.5       ],
           [ 2.5       ,  5.        , 10.        ]])



If we build bins from ‘scratch’ as above, when we initiate the OPC
model, we need to only include the bins as an argument:

.. code:: ipython3

    dylos = opcsim.OPC(bins=bins)

Define OPC Counting Efficiency
------------------------------

The last argument of interest to the OPC model is the counting
efficiency (``ce``). The counting efficiency argument must be a callable
function that accepts the particle diameter and returns a float. By
default, counting efficiency is set to return :math:`\eta=1` at all
diameters. You can provide any function you want.

Let’s define some counting efficiency functions that we can then
incorporate into various simulated OPCs:

.. code:: ipython3

    # Define efficiency based on an exponential function
    η_exp = lambda dp: 1 - np.exp(-5*dp)
    
    # Define efficiency based on a tanh function
    η_tanh = lambda dp: np.tanh(2*dp)
    
    # Define a function that rises linearly from 100nm to 1um, and then stays at 1
    η_linear = lambda dp: [np.piecewise(i, [i < 1., i >= 1.], [i, 1]) for i in dp]

Let’s go ahead and visualize these functions really quick to get a
better idea

.. code:: ipython3

    # Create an array of diameters
    diams = np.logspace(-2,1, 50)
    
    fig, ax = plt.subplots(1)
    
    ax.plot(diams, η_exp(diams), marker='o', label="$\eta=1-exp(-5D_p)$")
    ax.plot(diams, η_tanh(diams), marker='*', label="$\eta=tanh(2D_p)$")
    ax.plot(diams, η_linear(diams), marker='^', label="$\eta=linear$")
    
    ax.semilogx()
    
    sns.despine(offset=5)
    
    ax.set_xlabel("Diameter")
    ax.set_ylabel("Counting Efficiency")
    
    # Move the legend
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    
    ax.xaxis.set_major_formatter(mticks.FormatStrFormatter("%.3g"))



.. image:: model_files/model_18_0.png


Now that we have a better understanding of what the counting efficiency
function looks like (and how you can define your own), let’s go ahead
and show how to build an OPC that uses one of these functions.

Let’s go ahead and build a 10-bin OPC that uses the tanh counting
efficiency from above:

.. code:: ipython3

    opc_tanh = opcsim.OPC(n_bins=10, ce=η_tanh)

That more or less covers how we build an OPC. Next, how do we determine
what an OPC “sees” given an aerosol distribution?

Evaluate the OPC for a Given ``AerosolDistribution``
====================================================

To evaluate the OPC, we need to determine how many particles the OPC
‘sees’ in each size bin. Once we have this value, we can convert to
surface area, volume, or mass in order to compare to the true amount of
mass present in the underlying aerosol distribution.

There are two methods we use to do this:

1. ``simple`` method

   The simple method means we evaluate the PDF of the aerosol
   distribution at each bin midpoint. Depending on the ``weight`` and
   ``base`` we are evaluating at, this returns the
   :math:`d[weight]/d[base]D_p` value at the given bin. We take into
   account the counting efficiency by multiplying this value by the
   ``ce`` function evaluated at the midpoint diameter for each bin.
   Mathematically, this would be represented as:

   .. math:: \frac{d[weight]}{d[base]D_p}=\sum_{i=1}^{n_{bins}}PDF(D_{p,midpoint})*CE(D_{p,midpoint})

2. ``subint`` method

   The subintegration method takes a more continuous approach; the total
   number of particles in each bin is calculated by integrating the
   product of the CDF and the counting efficiency function within each
   individual bin. This provides a more “accurate” result. Essentially,
   if you assume the OPC has 100% counting efficiency, this would return
   the actual number of particles present in the given bin.

We assume that an OPC “sees” particle number concentration, and not some
correlation to particle volume. Thus, each evaluation is completed by
first evaluating the aerosol distribution in number-weighted space, and
then converting to number, surface area, or volume by multiplying by the
respective multiplier. The multiplier is determined at the bin midpoint,
which is important.

The ``opcsim`` library provides a few ways to obtain these values.

``opcsim.OPC.evaluate``
-----------------------

The ``opcsim.OPC.evaluate`` method returns an array of values where each
value is the number of {particles, surface area, volume} in each bin. It
will return data in the format :math:`d[weight]/d[base]D_p` where the
default is to return :math:`dN/dlogD_p` (``weight='number'``,
``base='log10'``). It can also be evaluated with either the ``simple``
evaluation method or the ``subint`` evaluation method depending on the
``method`` keyword argument provided.

For example, to evaluate a 5-bin OPC and return :math:`dN/dlogD_p`
values for each bin using the default ``subint`` method, we would do the
following:

.. code:: ipython3

    # Build a 5-bin OPC
    opc = opcsim.OPC(n_bins=5, dmin=0.3, dmax=2.5)
    
    # load the urban distribution
    urban = opcsim.load_distribution("Urban")
    
    # evaluate the number-weighted distribution
    opc.evaluate(distribution=urban)




.. parsed-literal::

    array([3.32717067e+02, 4.44738784e+01, 2.75920424e+00, 7.85362827e-02,
           1.01796109e-03])



To compare to the ``simple`` method, we can grab that data as well:

.. code:: ipython3

    opc.evaluate(urban, method='simple')




.. parsed-literal::

    array([3.04815785e+02, 3.57911500e+01, 1.87041608e+00, 4.33304183e-02,
           4.44549802e-04])



As you can see, they are similar, but not exactly the same. What if we
want to grab :math:`dV/dlogD_p`?

.. code:: ipython3

    opc.evaluate(urban, weight='volume')




.. parsed-literal::

    array([8.88552496e+00, 4.23842272e+00, 9.38370544e-01, 9.53129938e-02,
           4.40863558e-03])



``opcsim.OPC.number``
---------------------

Although the log-weighted values are ideal for visualization, when it
comes to evaluating the OPC performance, we want the actual number of
particles, surface area, or volume within each bin. To get this data, we
could either multiply the above results by the log difference of the
bins, or we can use one of the other methods made available.

The ``opcsim.OPC.number`` method returns the total number of particles
the OPC “sees” in each bin per a given distribution. You can also access
the “True” number of particles in each bin (i.e. the integrated CDF of
the underyling aerosol distribution) by changing the ``measured``
argument to be ``False``.

For example, let’s grab the total number of particles/cc in each bin of
the previous OPC per the Urban distribution:

.. code:: ipython3

    opc.number(urban)




.. parsed-literal::

    array([6.12744230e+01, 8.19047626e+00, 5.08145402e-01, 1.44635364e-02,
           1.87471533e-04])



``opcsim.OPC.surface_area``
---------------------------

Similar to the ``number`` method above, we can do the same for surface
area.

To get the surface area within each bin, we do the following:

.. code:: ipython3

    opc.surface_area(urban)




.. parsed-literal::

    array([2.64749631e+01, 8.26404722e+00, 1.19728940e+00, 7.95816850e-02,
           2.40880403e-03])



``opcsim.OPC.volume``
---------------------

Similar to the ``number`` and ``surface_area`` methods above, we can do
the same for volume.

To get the volume within each bin, we do the following:

.. code:: ipython3

    opc.volume(urban)




.. parsed-literal::

    array([1.63639160e+00, 7.80563826e-01, 1.72813839e-01, 1.75531984e-02,
           8.11910864e-04])



Plotting OPC Response to the Urban Distribution
-----------------------------------------------

Now that we know how to evaluate the response of an OPC to the urban
distribution, how can we easily visualize it? Well, we have the handy
function ``opcsim.plots.histplot`` to do that! All we need is the data
to plot (evaluated PDF) and the OPC bins.

Let’s go ahead and plot the response of a 10-bin OPC to the Urban
Aerosol Distribution:

.. code:: ipython3

    # Set the 10-bin OPC
    opc = opcsim.OPC(n_bins=10, dmin=0.3, dmax=2.5)
    
    # Load the urban distribution
    urban = opcsim.load_distribution("Urban")
    
    # Plot
    ax = opcsim.plots.histplot(opc.evaluate(urban), opc.bins)
    
    ax.set_ylabel("$dN/dlogD_p$")
    
    # Remove the spine
    sns.despine()



.. image:: model_files/model_34_0.png


Why don’t we go ahead and overlay the distribution itself:

.. code:: ipython3

    # Plot
    ax = opcsim.plots.histplot(opc.evaluate(urban), opc.bins)
    
    # Add the distribution to the plot
    ax = opcsim.plots.pdfplot(urban, ax=ax)
    
    sns.despine()



.. image:: model_files/model_36_0.png


The above plots are in number-space. The primary use of these low-cost
sensors is to estimate mass, so why don’t we go ahead and plot this in
volume space?

.. code:: ipython3

    # Plot
    ax = opcsim.plots.histplot(opc.evaluate(urban, weight='volume'), opc.bins)
    
    # Add the distribution to the plot
    ax = opcsim.plots.pdfplot(urban, weight='volume', ax=ax)
    
    ax.set_xlim(0.01, 10)
    
    sns.despine()



.. image:: model_files/model_38_0.png


Each of these plots uses the ``method='subint'`` integration method. How
does it change if we use the ``simple`` method instead?

.. code:: ipython3

    # Plot
    ax = opcsim.plots.histplot(opc.evaluate(urban, weight='volume'), opc.bins)
    ax = opcsim.plots.histplot(opc.evaluate(urban, weight='volume', method='simple'), opc.bins, ax=ax)
    
    # Add the distribution to the plot
    ax = opcsim.plots.pdfplot(urban, weight='volume', ax=ax)
    
    # Add a legend and set limits
    ax.legend(["Urban PDF", "subint", "simple"], bbox_to_anchor=(1.5, 1.05))
    ax.set_xlim(0.01, 10)
    
    sns.despine()



.. image:: model_files/model_40_0.png


So it doesn’t look too different from this picture, but it can have
reasonable impacts. That should be a fairly in depth introduction to
setting up, evaluating, and visualizing a simulated OPC.
