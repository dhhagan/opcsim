"""Utility functions
"""
import numpy as np
import math

def make_bins(dmin, dmax, n_bins, base='log'):
    """Returns a 3xn array of bin diameters.

    Build a 3xn array of bin diameters that can then be fed directly into
    the OPC class. The default behaviour is to space the bins equally on a log10
    basis.

    Parameters
    ----------

    dmin : float
        Minimum particle diameter in microns
    dmax : float
        Maximum particle diameter in microns
    n_bins : int
        Number of bins
    base : {'log' | 'none' | None}

    Returns
    -------
    bins : array
        Returns a 3xn_bins array in the format of [left edge, midpoint, right edge]


    Examples
    --------

    Build a set of bins for an OPC with dmin=0.5, dmax=2.5, and 3 bins:

    >>> bins = opcsim.utils.make_bins(0.5, 2.5, 3)

    """
    if base is None:
        base = 'none'

    bins = np.zeros((n_bins, 2)) * np.nan

    # Set the left edge of bin0 and right edge of the last bin
    bins[0, 0] = dmin
    bins[-1, -1] = dmax

    mult = 1. / (( np.log10(dmax) - np.log10(dmin) ) / n_bins )

    for i in range(n_bins):
        bins[i, -1]  = 10**(np.log10( bins[i, 0]) + ( 1. / mult ))

        if i < n_bins - 1:
            bins[i + 1, 0] = bins[i, -1]

    return midpoints(bins, base=base)

def midpoints(bins, base='log'):
    """Returns a 3xn array of bin diameters.

    Build a 3xn array of bin diameters from a 2xn array that can then be fed
    directly into the OPC class. The default behaviour is to space the bins
    equally on a log10 basis.

    Parameters
    ----------

    bins : 2xn array
        A 2xn array with bin boundaries.
    base : {'log' | 'none' | None}

    Returns
    -------
    bins : array
        Returns a 3xn_bins array in the format of [left edge, midpoint, right edge]


    Examples
    --------

    Build a set of bins for an OPC similar to a Dylos DC1100 Pro:

    >>> arr = np.array([0.5, 2.5], [2.5, 10])
    >>> bins = opcsim.utils.midpoints(arr)

    """
    base = 'base' if base is None else base

    if bins.shape[1] != (2 or 3):
        raise ValueError("Invalid bins array. Must be either 2xn or 3xn.")

    tmp = np.zeros((bins.shape[0], 3)) * np.nan
    tmp[:, 0] = bins[:, 0]
    tmp[:, 2] = bins[:, -1]

    for i in range(tmp.shape[0]):
        if base.lower() == 'log':
            tmp[i, 1] = 10**np.log10([tmp[i,0], tmp[i, 2]]).mean()
        else:
            tmp[i, 1] = np.mean([tmp[i, 0], tmp[i, 2]])

    return tmp
