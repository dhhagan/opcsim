
.. _scoring_tutorial:

Evaluating the Results of a Simulated OPC
=========================================

To determine the efficacy of our simulated OPC’s, we must derive some
sort of metric to use to compare to other modeled OPC’s and other
distributions. Perusing through the current literature, you will find
two primary methods that encapsulate how this is done.

1. Number-Volume Correlation Method (``opcsim.metrics.nv_score``)
-----------------------------------------------------------------

The number-correlation method is used throughout the low-cost sensor
community, especially since it is quite often the only option.
Essentially, you simply correlate the output variable of the OPC
(typically something related to number concentration) to PM2.5 (or other
Mass loading) from a reference instrument.

In doing so, you are inherently assuming the underlying particle
distribution does not vary too much from whatever distribution you
calibrated against. To calculate the number-volume correlation, we use
the following relationship:

.. math:: \frac{N}{V}=\frac{Number\;of\;particles\;seen\;by\;the\;OPC}{total\;integrated\;volume\;in\;the\;underlying\;distribution}

To compute this value, we can use the ``opcsim.metrics.nv_score``
function which requires just the ``opcsim.OPC`` model and the
``opcsim.AerosolDistribution``. Optionally, you can limit the minimum
and maximum diameters to evaluate the volume CDF under; however, this
will not change the number calculation from the OPC.

For example, let’s compute the ``nv_score`` for the 1-bin OPC that is
default against the urban distribution.

.. code:: ipython3

    # Make imports
    import opcsim
    import numpy as np
    
    # turn off warnings temporarily
    import warnings
    warnings.simplefilter('ignore')

.. code:: ipython3

    # Build the OPC
    opc = opcsim.OPC()
    
    # load the distribution
    urban = opcsim.load_distribution("Urban")
    
    # compute the nv_score
    opcsim.metrics.nv_score(opc, urban)




.. parsed-literal::

    0.9543539696612989



What is the meaning of this number? Nothing, really. On it’s own, this
value is pretty useless. It is the slope between the number of particles
your OPC “sees” and the total volume under 2.5 microns in the urban
distribution. Now, if we change the distribution, what happens? Let’s
try computing the same score for the rural distribution:

.. code:: ipython3

    rural = opcsim.load_distribution("Rural")
    
    opcsim.metrics.nv_score(opc, rural)




.. parsed-literal::

    1.2993300334451958



Ahh! Nearly a 30% difference just by slightly perturbing the underlying
particle size distribution.

2. Volume-Volume Correlation Method (``opcsim.metrics.vv_score``)
-----------------------------------------------------------------

The volume-correlation method is used by some of the more robust
low-cost OPCs (like the Alphasense OPC-N2) which have many size bins.
The idea behind the volume method is to actually do a step-wise
integration across the bins to calculate the total volume/mass present.
Unfortunately, most of these sensors still have relatively high
:math:`dmin` values which means they are blind to a large enough
fraction of the mass that it makes a difference.

Unlike the number-correlation method described above, this method has a
chance at observing changes in the underlying size distribution. To
score this method, we take the ratio of the total volume calculated
across the bins of the OPC to the total volume present in the underlying
particle size distribution.

To compute this value, we can use the ``opcsim.metrics.vv_score``
function which requires just the ``opcsim.OPC`` model and the
``opcsim.AerosolDistribution``. Optionally, you can limit the minimum
and maximum diameters to evaluate the volume CDF under; however, this
will not change the volume calculation from the OPC.

For example, let’s compute the ``vv_score`` for the 1-bin OPC that is
default against the urban distribution.

.. code:: ipython3

    opcsim.metrics.vv_score(opc, urban)




.. parsed-literal::

    0.698349981739988



What does this value mean? Well, it is simply a ratio of OPC volume to
Actual Volume, so this is the fraction of volume seen by the OPC.

How much does it change when the distribution changes?

.. code:: ipython3

    opcsim.metrics.vv_score(opc, rural)




.. parsed-literal::

    0.9507867457739008



While these numbers look relatively great, it may just be a fluke! What
happens when we score this method for a variety of different OPCs?

.. code:: ipython3

    models = []
    
    dmin = 0.3
    dmax = 2.5
    
    for i in range(1, 10):
        models.append(("{}-Bin OPC".format(i), opcsim.OPC(n_bins=i, dmin=dmin, dmax=dmax)))
    
    for model in models:
        nv = opcsim.metrics.nv_score(model[1], urban)
        vv = opcsim.metrics.vv_score(model[1], urban)
        
        print ("\n{}".format(model[0]))
        print ("\tN/V = {:.3f}".format(nv))
        print ("\tV/V = {:.3f}".format(vv))


.. parsed-literal::

    
    1-Bin OPC
    	N/V = 12.829
    	V/V = 4.363
    
    2-Bin OPC
    	N/V = 12.829
    	V/V = 0.918
    
    3-Bin OPC
    	N/V = 12.829
    	V/V = 0.605
    
    4-Bin OPC
    	N/V = 12.829
    	V/V = 0.516
    
    5-Bin OPC
    	N/V = 12.829
    	V/V = 0.478
    
    6-Bin OPC
    	N/V = 12.829
    	V/V = 0.458
    
    7-Bin OPC
    	N/V = 12.829
    	V/V = 0.446
    
    8-Bin OPC
    	N/V = 12.829
    	V/V = 0.439
    
    9-Bin OPC
    	N/V = 12.829
    	V/V = 0.433


The number correlation method stayed the same for all OPCs! Why? Well,
right now we have simulated each of these with a counting efficiency of
1 which means they see 100% of the particles. While this doesn’t change
the number of particles we see in total, it does change the volume!

This ends the introduction to using the metrics to score your
OPC/Distribution model.
