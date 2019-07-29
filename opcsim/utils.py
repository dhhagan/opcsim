"""Utility functions
"""
import numpy as np
from scipy.optimize import curve_fit
import math

# assuming STP
RHO_H20 = 0.997


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


def k_kohler(diam_dry, kappa=0., rh=0.):
    """Calculate the wet diameter of a particle based on the hygroscopic growth 
    parameter, kappa (k-Kohler theory).

    .. math::

        D_w=D_d*\sqrt[3]{1 + \\frac{a_w}{1-a_w}\kappa_{eff}}

    Parameters
    ----------
    diameter_dry: float
        The dry diameter in any units (nm or um most likely)
    kappa: float, optional
        The effective kappa value
    rh: float: optional
        The relative humidity as a percentage (0.0-100.0)

    Returns
    -------
    diameter_wet
        The wet diameter in the same units supplied for the dry diameter

    """
    # calculate the water activity
    aw = rh / 100.

    return diam_dry * math.pow(1 + kappa * (aw / (1 - aw)), 1./3.)


def rho_eff(rho, weights=None, diams=None):
    """Calculate the effective density of a particle by calculating the 
    wet and dry percentages and taking the weighted sum. Alternatively,
    an array of diameters can be passed which will be used to 
    calculate the volumetric weights.

    Parameters
    ----------
    rho: ndarray of floats
        An array of particle densities.
    weights: ndarray of floats
        An array of volumetric weight percentages.
    diams: ndarray of floats
        An array of diameters for each species.
    
    Returns
    -------
    rho_eff: float
        The weighted density of the wet particle.
    """
    rho = np.asarray(rho)

    # if diams are present, compute their weights
    if diams is not None:
        diams = np.asarray(diams)
        weights = (diams**3) / (diams**3).sum()

    weights = np.asarray(weights)

    return (weights * rho).sum()


def k_eff(kappas, weights=None, diams=None):
    """Calculate the effective k-kohler coefficient from 
    an array of kappa values and their weights. Alternatively,
    an array of diameters can be passed which will be used to 
    calculate the volumetric weights.

    .. math::

        \kappa=\sum_{i=1} \epsilon_i \kappa_i

    Parameters
    ----------
    kappas: ndarray
        An array of k-kohler coefficients.
    weights: ndarray, optional
        An array of volumetric weights.
    diams: ndarray
        An array of diameters used to calculate the weights.
    
    Returns
    -------
    k_eff: float
        The effective k-kohler coefficient.
    """
    kappas = np.asarray(kappas)

    # if diams are present, compute their weights
    if diams is not None:
        diams = np.asarray(diams)
        weights = (diams**3) / (diams**3).sum()

    weights = np.asarray(weights)

    return (kappas * weights).sum()


def ri_eff(species, weights=None):
    """Calculate the effective refractive index for an 
    array of refractive indices and their respective weights. 
    Alternatively, an array of diameters can be passed which 
    will be used to calculate the volumetric weights.

    .. math::

        n_{eff}=\sum_{i=1}^{N} \\frac{V_i}{V_{total}}*n_i

    Parameters
    ----------
    species: ndarray
        An array of refractive indices.
    weights: ndarray, optional
        An array of volumetric weights.
    
    Returns
    -------
    ri_eff: float
        The effective refractive index.
    """
    species = np.asarray(species)
    weights = np.asarray(weights)

    # calculate the real and imag parts separately
    real = (species.real * weights).sum()
    imag = (species.imag * weights).sum()

    return complex(real, imag)


def power_law_fit(diams, cscat, fit_kws={}):
    """Generate a power-law fit (linear in log-log space)
    between bin midpoints and the scattering cross-section.

    Parameters
    ----------
    diams: ndarray
        An array of diameters corresponding to the boundaries of each bin.
    cscat: ndarray
        An array of Cscat values corresponding to the boundaries of each bin.
    fit_kws: dict
        A dictionary of keyword arguments that is passed directly to scipy.curve_fit
        for the optimization. Please see the scipy.optimize.curve_fit docs for details.

    Returns
    -------
    rv: ndarray
        An array of fitted Cscat values for each bin boundary diameter.

    Examples
    --------

    """
    def f(dp, a, b):
        return a*np.power(dp, b)
    
    # fit the function
    popt, _ = curve_fit(f, diams, cscat, **fit_kws)

    return f(diams, *popt)


def squash_dips(cscat_vals):
    """Remove any dips in an array by interpolating 
    around those points. If there are no dips and all 
    points are monotonically increasing, the same 
    array should be returned.

    Parameters
    ----------
    cscat: ndarray
        An array of scattering cross-section values.
    
    Returns
    -------
    rv: ndarray
        An array of smoothed scattering cross-section values.

    """
    cpy = cscat_vals.copy()

    if (np.diff(cpy) < 0).any():
        idx = np.where(np.diff(cpy) < 0)[0]

        if idx.shape[0] > 0:
            for i in np.nditer(idx):
                cpy[i] = np.mean([cpy[i-1], cpy[i+1]])
    
    return cpy
