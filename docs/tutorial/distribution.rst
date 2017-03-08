
.. _distribution_tutorial:


Build an Aerosol Distribution
=============================

An aerosol distribution...

.. code:: ipython3

    import opcsim
    import pandas as pd
    import numpy as np
    
    import seaborn as sns

Urban Aerosol Distribution
--------------------------

To begin, we can initiate an instance of the ``AerosolDistribution``
class using parameters from Seinfeld and Pandis:

+--------+-------------+------------------+-----------------------+
| Mode   | :math:`N`   | :math:`D_{pg}`   | :math:`log\sigma_i`   |
+========+=============+==================+=======================+
| I      | 7100        | 0.0117           | 0.232                 |
+--------+-------------+------------------+-----------------------+
| II     | 6320        | 0.0373           | 0.250                 |
+--------+-------------+------------------+-----------------------+
| III    | 960         | 0.151            | 0.204                 |
+--------+-------------+------------------+-----------------------+

.. code:: ipython3

    d = opcsim.AerosolDistribution("Urban")
    
    d.add_mode(7100, 0.0117, 10**0.232, "Mode I")
    d.add_mode(6320, 0.0373, 10**0.25, "Mode II")
    d.add_mode(960, 0.151, 10**0.204, "Mode III")
    
    d




.. parsed-literal::

    AerosolDistribution: Urban
    
    Mode	N		GM	GSD
    Mode I	7.10e+03	0.012	1.706
    Mode II	6.32e+03	0.037	1.778
    Mode III	9.60e+02	0.151	1.600



.. code:: ipython3

    d.pdf(1)




.. parsed-literal::

    0.57136840961758029



.. code:: ipython3

    d.pdf(np.linspace(0.01, 1., 5))




.. parsed-literal::

    array([  1.24309860e+04,   1.02054540e+03,   6.94029532e+01,
             5.44233259e+00,   5.71368410e-01])



.. code:: ipython3

    opcsim.load_distribution("Urban")




.. parsed-literal::

    AerosolDistribution: Urban
    
    Mode	N		GM	GSD
    Mode I	7.10e+03	0.012	1.706
    Mode II	6.32e+03	0.037	1.778
    Mode III	9.60e+02	0.151	1.600



