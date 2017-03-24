"""Contains the scoring algorithms used in the model.
"""

import numpy as np

def nv_score(model, distribution, dmin=0.0, dmax=2.5, **kwargs):
    """Calculate and return the number-to-volume ratio.

    The total number of particles is calculated by calculating the total number
    of particles in each individual bin, and then summing them. The total volume
    in the distribution is calculated by integrating the Volume-weighted CDF
    between 0 and `dmax` microns.

    Parameters
    ----------
    model : OPC
        A valid OPC model describing an OPC that can be evaluated.
    distribution : AerosolDistribution
        A valid AerosolDistribution instance that can be evaluated.
    dmin : float
        The minimum particle size to integrate the CDF under. Default is 0.0
        microns.
    dmax : float
        The maximum particle size to integrate the CDF under. Default is 2.5
        microns.

    Returns
    -------
    N/V : float
        Returns the number-to-volume ratio as a single float.

    Examples
    --------

    Compute the number-to-volume ratio for a 2-bin OPC on the Urban distribution

    >>> opc = opcsim.OPC(n_bins=2)
    >>> urban = opcsim.load_distribution("Urban")
    >>> n_v = opcsim.metrics.nv_score(opc, urban)

    """
    # evaluate the total number of particles in each bin (then sum)
    total_number = model.number(distribution, **kwargs).sum()

    # evaluate the total volume in the distribution < dmax
    total_volume = distribution.cdf(weight='volume', dmax=dmax)

    return total_number / total_volume

def vv_score(model, distribution, dmin=0.0, dmax=2.5, **kwargs):
    """Calculate and return the volume-to-volume ratio.

    The total volume of particles per the OPC is calculated by calculating the
    total number of particles in each individual bin, and then multiplying each
    bin by a 'volume-factor'. The sum of individual bin volumes is then used.
    The total volume in the distribution is calculated by integrating the
    Volume-weighted CDF between 0 and `dmax` microns.

    Parameters
    ----------
    model : OPC
        A valid OPC model describing an OPC that can be evaluated.
    distribution : AerosolDistribution
        A valid AerosolDistribution instance that can be evaluated.
    dmin : float
        The minimum particle size to integrate the CDF under. Default is 0.0
        microns.
    dmax : float
        The maximum particle size to integrate the CDF under. Default is 2.5
        microns.

    Returns
    -------
    V/V : float
        Returns the volume-to-volume ratio as a single float.

    Examples
    --------

    Compute the number-to-volume ratio for a 2-bin OPC on the Urban distribution

    >>> opc = opcsim.OPC(n_bins=2)
    >>> urban = opcsim.load_distribution("Urban")
    >>> v_v = opcsim.metrics.vv_score(opc, urban)

    """
    # evaluate the total number of particles in each bin (then sum)
    measured_volume = model.volume(distribution, **kwargs).sum()

    # evaluate the total volume in the distribution < dmax
    total_volume = distribution.cdf(weight='volume', dmin=dmin, dmax=dmax)

    return measured_volume / total_volume
