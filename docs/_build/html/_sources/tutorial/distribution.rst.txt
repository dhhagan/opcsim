
.. _distribution_tutorial:


The following tutorial will show you how an aerosol distribution is
represented in the opcsim software. You will learn how to import sample
datasets and how to create your own distribution from scratch.
Additional visualization tools are also discussed.

First, we import the python libraries we need and set the styles used
for plotting throughout this tutorial.

.. code:: ipython3

    # Make imports
    import opcsim
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    %matplotlib inline
    
    # turn off warnings temporarily
    import warnings
    warnings.simplefilter('ignore')
    
    # Let's set some default seaborn settings
    sns.set(context='notebook', style='ticks', palette='dark', font_scale=1.75, rc={'figure.figsize': (10,5)})

The Aerosol Distribution
========================

For the purpose of evaluating the performance of low-cost optical
particle counters, we are going to assume that every aerosol
distribution can be described as the sum of :math:`n` lognormal
distributions (S+P 8.54). Thus, it follows that:

.. math:: n_N^o(logD_p)=\sum_{i=1}^n \frac{N_i}{\sqrt{2\pi} * log\sigma_i}exp\Big(-\frac{(logD_p - logD_{pi})^2}{2log^2\sigma_i}\Big)

where :math:`N_i` is the number concentration, :math:`D_{pi}` is the
median particle diameter, and :math:`\sigma` is the standard deviation.
Thus, we need :math:`3n` parameters to describe a complete aerosol
distribution.

Using the ``opcsim.AerosolDistribution`` class, we can build our own
distributions by defining each aerosol mode as its own lognormal mode.

Ex: Initialize an Aerosol Distribution with a Single Mode
(:math:`N=1000`, :math:`D_{pg}=100\;nm`, :math:`\sigma=1.5`)

.. code:: ipython3

    # Initialize a distribution
    sample = opcsim.AerosolDistribution()
    
    # Add a mode with N=1000, GM=0.1, GSD=1.5
    sample.add_mode(n=1000, gm=0.1, gsd=1.5, label="Mode I")

Most aerosol distributions are composed of multiple lognormal modes.
Table 8.3 in Seinfeld and Pandis shows some sample parameters for
several different distributions. The urban aerosol distribution can be
described as follows:

+--------+-------------+------------------+-----------------------+
| Mode   | :math:`N`   | :math:`D_{pg}`   | :math:`log\sigma_i`   |
+========+=============+==================+=======================+
| I      | 7100        | 0.0117           | 0.232                 |
+--------+-------------+------------------+-----------------------+
| II     | 6320        | 0.0373           | 0.250                 |
+--------+-------------+------------------+-----------------------+
| III    | 960         | 0.151            | 0.204                 |
+--------+-------------+------------------+-----------------------+

How would we go about building this distribution? We can add as many
modes as we would like, just as we did above. Also, if you look at the
API documentation for the ``opcsim.AerosolDistribution`` class, we see
that you can add a label for the distribution as an argument upon
initiation of the class instance.

.. code:: ipython3

    urban = opcsim.AerosolDistribution("Urban")
    
    urban.add_mode(7100, 0.0117, 10**0.232, "Mode I")
    urban.add_mode(6320, 0.0373, 10**0.25, "Mode II")
    urban.add_mode(960, 0.151, 10**0.204, "Mode III")

To make things even easier, there are several "sample" distributions
included in the package and can be accessed via the
``opcsim.load_distribution`` function. The parameters for the
distributions are taken directly from Seinfeld and Pandis Table 8.3. To
read in the distribution, simply provide the name of the type of
distribution you would like to use. Options include:

-  Urban
-  Marine
-  Rural
-  Remote Continental
-  Free Troposphere
-  Polar
-  Desert

To read in the urban distribution, we would do the following:

.. code:: ipython3

    urban = opcsim.load_distribution("Urban")

Probability Distribution Function
=================================

Number Distribution
-------------------

Aerosol distributions are typically depicted using the probabliity
distribution function. In Seinfeld and Pandis, they refer to it as the
"Number Distribution Function". When plotted in log-space (i.e.
:math:`dN/dlogD_p`), the area under the curve is the aerosol number
concentration.

Mathematically, the PDF in number-space looks like the following:

.. math:: n_N^o(logD_p)=\frac{dN}{dlogD_p}=\frac{N_t}{\sqrt{2\pi} \; log\sigma_g}exp\Big(-\frac{(logD_p - logD_{pg})^2}{2log^2\sigma_g}\Big)

Surface Area Distribution
-------------------------

It is also quite useful to look at the surface area and volume
distributions. The surface area probability distribution can easily be
obtained by relating to the number probability distribution in the
following way:

.. math:: n_S^o(logD_p)=\pi D_p^2 n_N^o(logD_p)=\frac{dS}{dlogD_p}=\frac{\pi D_p^2 N_t}{\sqrt{2\pi} \; log\sigma_g}exp\Big(-\frac{(logD_p - logD_{pg})^2}{2log^2\sigma_g}\Big)

Volume Distribution
-------------------

Likewise, for the volume distribution, we get:

.. math:: n_V^o(logD_p)=\frac{\pi}{6} D_p^3 n_N^o(logD_p)=\frac{dV}{dlogD_p}=\frac{\pi D_p^3 N_t}{6\sqrt{2\pi} \; log\sigma_g}exp\Big(-\frac{(logD_p - logD_{pg})^2}{2log^2\sigma_g}\Big)

``opcsim`` provides the ``AerosolDistribution.pdf`` method to easily
calculate the distribution at any particle diameter. The arguments of
the function are the particle diameter (``dp``), the base (``none``,
``log``, or ``log10``), the weight (``number``, ``surface``, ``volume``,
or ``mass``), and an optional ``mode`` parameter in case you would like
to examine only one of the modes of the distribution at a time. The
default arguments are set to be the most useful/common ones (i.e.
``weight='number'``, ``base='log10'``). If calculating the mass-weighted
PDF, you can also provide an optional keyword argument ``rho``; the
default value for particle density is :math:`1\;gcm^{-3}`.

To calculate the number probability for the urban aerosol distribution
at 0.1 micron, we do the following:

.. code:: ipython3

    urban.pdf(0.1)




.. parsed-literal::

    3606.2139576648124



This gives us the number concentration probability at 1 micron in units
of :math:`particles\;cm^{-3}`. We can also calculate a whole range of
values by providing an array for the ``dp`` value:

.. code:: ipython3

    urban.pdf(np.array([0.1, 0.2, 0.3]))




.. parsed-literal::

    array([ 3606.21395766,  1712.82519467,   659.56432207])



To calculate the volume-weighted PDF at some particle diameter
(:math:`dV/dlogD_p`), we could do the following:

.. code:: ipython3

    urban.pdf(0.1, weight='volume')




.. parsed-literal::

    1.8882092127787917



This returns :math:`dV/dlogDp` at particle diameter
:math:`D_p=0.1\;\mu m` in units of :math:`\mu m^3 cm^{-3}`.

Visualizing the Distribution
============================

Visualizing the PDF for an aerosol distribution is extremely helpful.
The function ``opcsim.plots.pdfplot`` has been included to make this
simple.

To plot the pdf of an aerosol distribution, the only required input is
the ``opcsim.AerosolDistribution`` object. The function returns a
matplotlib axis object which makes it extremely easy to add to modify
the plot using normal matplotlib syntax.

Let's plot the urban distribution we built earlier.

.. code:: ipython3

    ax = opcsim.plots.pdfplot(urban)



.. image:: distribution_files/distribution_17_0.png


We can also go ahead and plot each individual mode along with the entire
distribution using the ``with_modes`` argument:

.. code:: ipython3

    ax = opcsim.plots.pdfplot(urban, with_modes=True)
    
    ax.legend(loc='best')
    
    sns.despine()



.. image:: distribution_files/distribution_19_0.png


Still staying in number space, we can go ahead and plot all of the
available sample distributions to get a feel for just how different they
are!

.. code:: ipython3

    fig, ax = plt.subplots(1, figsize=(14,7))
    
    # Iterate over every sample in the library
    for i, sample in enumerate(opcsim.distributions.DISTRIBUTION_DATA.keys()):
        # Load the sample dataset
        _sample = opcsim.load_distribution(sample)
        
        # if we've used more colors than we have available in this palette, change the linestyle
        ls = '-' if i < 6 else '--'
        
        opcsim.plots.pdfplot(_sample, ax=ax, plot_kws={'linestyle': ls})
        
    # Add a legend
    ax.legend(loc='best')
    
    # remove the spine
    sns.despine()



.. image:: distribution_files/distribution_21_0.png


Finally, we can also go ahead and look at one distribution in number,
surface area, and volume weighted views:

.. code:: ipython3

    fig, ax = plt.subplots(3, figsize=(12,9), sharex=True)
    
    opcsim.plots.pdfplot(urban, weight='number', ax=ax[0])
    opcsim.plots.pdfplot(urban, weight='surface', ax=ax[1])
    opcsim.plots.pdfplot(urban, weight='volume', ax=ax[2])
    
    #fig.subplots_adjust(hspace=0)
    
    ax[0].set_ylabel("Number")
    ax[1].set_ylabel("Surface Area")
    ax[2].set_ylabel("Volume")
    
    plt.tight_layout(h_pad=0)
    plt.show()



.. image:: distribution_files/distribution_23_0.png


Cumulative Distribution Function
================================

We can easily obtain the integrated value for number of particles, total
surface area, total volume, or total mass by integrating the correct
CDF.

Number CDF
----------

The total number of particles between two particle diameters can be
found by completing the following integration

.. math:: N_t=\int_{dmin}^{dmax}n_N(D_p)dD_p

Surface Area CDF
----------------

We can find the total particle surface area between two diameters using
the following integral:

.. math:: S_t=\pi \int_{dmin}^{dmax}D_p^2 n_N(D_p)dD_p

Volume CDF
----------

We can find the total particle volume between two diameters using the
following integral:

.. math:: V_t=\frac{\pi}{6} \int_{dmin}^{dmax}D_p^3 n_N(D_p)dD_p

To evaluate the CDF, we use the ``opcsim.AerosolDistribution.cdf``
method. For example, to evaluate the number of particles with diameter
less than :math:`D_p=2.5\;\mu m`, we do the following:

.. code:: ipython3

    urban.cdf(dmax=2.5)




.. parsed-literal::

    14379.999998896783



If we want to calculate the total number of particles within some size
range, we can add the ``dmin`` argument. For example, let's find the
total number of particles between 1 and 2.5 microns:

.. code:: ipython3

    urban.cdf(dmin=1, dmax=2.5)




.. parsed-literal::

    0.027425959342281203



What about the total volume of particles less than
:math:`D_p=1 \; \mu m`? (i.e. :math:`PM_1`)

.. code:: ipython3

    urban.cdf(dmax=1, weight='volume')




.. parsed-literal::

    5.4345309918199618



Last, how about the total mass of particles in the Urban distribution if
we set the particle density :math:`\rho=1.65\;gcm^{-3}`:

.. code:: ipython3

    urban.cdf(dmax=100, weight='mass', rho=1.65)




.. parsed-literal::

    9.0013585563083698



Although we wouldn't normally plot the CDF, we easily can to visualize
where most of the [number, surface area, mass] is within the
distribution using the ``opcsim.plots.cdfplot`` function:

.. code:: ipython3

    ax = opcsim.plots.cdfplot(urban)



.. image:: distribution_files/distribution_33_0.png


Lastly, we can plot the total volume CDF to get an idea of where the
mass is distributed:

.. code:: ipython3

    ax = opcsim.plots.cdfplot(urban, weight='mass', rho=1.65)



.. image:: distribution_files/distribution_35_0.png


