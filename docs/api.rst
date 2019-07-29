.. api_ref:

.. currentmodule:: opcsim

API Reference
=============

.. _distributions_api:

Aerosol Distributions
---------------------

.. rubric:: Aerosol Distribution Class

.. autosummary::
    :toctree: generated/

    opcsim.AerosolDistribution

.. rubric:: AerosolDistribution Methods

.. autosummary::
    :toctree: generated/

    opcsim.AerosolDistribution.add_mode
    opcsim.AerosolDistribution.pdf
    opcsim.AerosolDistribution.cdf

.. _models_api:

Models
------

OPC
^^^

.. rubric:: OPC Class

.. autosummary::
    :toctree: generated/

    opcsim.OPC

.. rubric:: OPC Methods

.. autosummary::
    :toctree: generated/

    opcsim.OPC.calibrate
    opcsim.OPC.evaluate
    opcsim.OPC.histogram
    opcsim.OPC.integrate


Nephelometer
^^^^^^^^^^^^

.. rubric:: Nephelometer Class

.. autosummary::
    :toctree: generated/

    opcsim.Nephelometer

.. rubric:: Nephelometer Methods

.. autosummary::
    :toctree: generated/

    opcsim.Nephelometer.calibrate
    opcsim.Nephelometer.evaluate


.. _plots_api:

Visualization
-------------

.. autosummary::
    :toctree: generated/

    opcsim.plots.histplot
    opcsim.plots.pdfplot
    opcsim.plots.cdfplot
    opcsim.plots.calplot


.. metrics_api:

Metrics Functions
-----------------

.. autosummary::
    :toctree: generated/

    opcsim.metrics.compute_bin_assessment

.. _utils_api:

Utility Functions
-----------------

.. autosummary::
    :toctree: generated/

    opcsim.load_distribution

    opcsim.utils.make_bins
    opcsim.utils.midpoints

    opcsim.utils.k_kohler
    opcsim.utils.rho_eff
    opcsim.utils.k_eff
    opcsim.utils.ri_eff

    
.. _mie_theory_api:

Mie Theory Calculations
-----------------------

.. autosummary::
    :toctree: generated/

    opcsim.mie.coef_pi_tau
    opcsim.mie.coef_ab
    opcsim.mie.s1s2
    opcsim.mie.cscat


.. _equations_api:

Equations
---------

.. autosummary::
    :toctree: generated/

    opcsim.equations.pdf.dn_ddp
    opcsim.equations.pdf.ds_ddp
    opcsim.equations.pdf.dv_ddp

    opcsim.equations.pdf.dn_dlndp
    opcsim.equations.pdf.ds_dlndp
    opcsim.equations.pdf.dv_dlndp

    opcsim.equations.pdf.dn_dlogdp
    opcsim.equations.pdf.ds_dlogdp
    opcsim.equations.pdf.dv_dlogdp

    opcsim.equations.cdf.nt
    opcsim.equations.cdf.st
    opcsim.equations.cdf.vt



