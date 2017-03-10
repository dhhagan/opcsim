import math
import numpy as np
from scipy.special import erf, erfc

def nt(n, gm, gsd, dmin=None, dmax=10.):
    """Evaluate the total number of particles between two diameters.

    The CDF of a lognormal distribution as calculated using equations 8.2 and
    8.39 from Seinfeld and Pandis.

    Mathematically, it is represented as:

    .. math::

        N_t(D_p)=∫_{dmin}^{dmax}n_N(D_p^*)dD^*_p

    Parameters
    ----------
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.
    dmin : float
        The minimum particle diameter in microns. Default value is 0.
    dmax : float
        The maximum particle diameter in microns. Default value is 10.

    Returns
    -------
    N | float
        Returns the total number of particles between dmin and dmax in units of
        [particles cm-3]

    See Also
    --------
    opcsim.equations.pdf.dn_ddp
    opcsim.equations.pdf.dn_dlndp
    opcsim.equations.pdf.dn_dlogdp

    Examples
    --------

    Integrate a sample distribution between 0 and 2.5 microns:

    >>> d = opcsim.AerosolDistribution()
    >>> d.add_mode(1e3, 100, 1.5, "mode 1")
    >>> n = opcsim.equations.cdf.nt(1e3, 0.1, 1.5, dmax=2.5)

    """
    res = (n/2.) * (1 + erf((np.log(dmax/gm)) / (np.sqrt(2) * np.log(gsd))))

    if dmin is not None:
        res -= (n/2.) * (1 + erf((np.log(dmin/gm)) / (np.sqrt(2) * np.log(gsd))))

    return res

def st(n, gm, gsd, dmin=None, dmax=10.):
    """Evaluate the total surface area of the particles between dmin and dmax.

    The CDF of a lognormal distribution as calculated using equation 8.11
    from Seinfeld and Pandis.

    Mathematically, it is represented as:

    .. math::

        S_t=π∫_{-∞}^{∞}D_p^2n_N^e(ln D_p)d lnD_p

    Parameters
    ----------
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.
    dmin : float
        The minimum particle diameter in microns. Default value is 0.
    dmax : float
        The maximum particle diameter in microns. Default value is 10.

    Returns
    -------
    SA | float
        Returns the total surface area of particles between dmin and dmax in
        units of [um2 cm-3]

    See Also
    --------
    opcsim.equations.pdf.ds_ddp
    opcsim.equations.pdf.ds_dlndp
    opcsim.equations.pdf.ds_dlogdp

    Examples
    --------

    Integrate a sample distribution between 0 and 2.5 microns:

    >>> d = opcsim.AerosolDistribution()
    >>> d.add_mode(1e3, 100, 1.5, "mode 1")
    >>> n = opcsim.equations.cdf.st(1e3, 0.1, 1.5, dmax=2.5)

    """
    res = (np.pi/2.)*n*(gm**2) * np.exp(2*(np.log(gsd)** 2)) * \
                erfc((np.sqrt(2) * np.log(gsd)) - (np.log(dmax/gm) / (np.sqrt(2) * np.log(gsd))))

    if dmin is not None:
        res -= (np.pi/2.)*n*(gm**2) * np.exp(2*(np.log(gsd)** 2)) * \
                    erfc((np.sqrt(2) * np.log(gsd)) - (np.log(dmin/gm) / (np.sqrt(2) * np.log(gsd))))

    return res

def vt(n, gm, gsd, dmin=None, dmax=10.):
    """Evaluate the total volume of the particles between dmin and dmax.

    The CDF of a lognormal distribution as calculated using equation 8.12
    from Seinfeld and Pandis.

    Mathematically, it is represented as:

    .. math::

        V_t=\\frac{π}{6}∫_{-∞}^{∞}D_p^3n_N^e(ln D_p)d lnD_p

    Parameters
    ----------
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.
    dmin : float
        The minimum particle diameter in microns. Default value is 0.
    dmax : float
        The maximum particle diameter in microns. Default value is 10.

    Returns
    -------
    V | float
        Returns the total volume of particles between dmin and dmax in
        units of [um3 cm-3]

    See Also
    --------
    opcsim.equations.pdf.dv_ddp
    opcsim.equations.pdf.dv_dlndp
    opcsim.equations.pdf.dv_dlogdp

    Examples
    --------

    Integrate a sample distribution between 0 and 2.5 microns:

    >>> d = opcsim.AerosolDistribution()
    >>> d.add_mode(1e3, 100, 1.5, "mode 1")
    >>> n = opcsim.equations.cdf.vt(1e3, 0.1, 1.5, dmax=2.5)

    """

    res = (np.pi/12.)*n*(gm**3) * np.exp(9./2.*(np.log(gsd)**2)) * \
                erfc((1.5*np.sqrt(2) * np.log(gsd)) - (np.log(dmax/gm) / (np.sqrt(2) * np.log(gsd))))

    if dmin is not None:
        res -= (np.pi/12.)*n*(gm**3) * np.exp(9./2.*(np.log(gsd)**2)) * \
                    erfc((1.5*np.sqrt(2) * np.log(gsd)) - (np.log(dmin/gm) / (np.sqrt(2) * np.log(gsd))))

    return res
