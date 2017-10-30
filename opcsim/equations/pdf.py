#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file contains all of the probability distribution functions
"""
import math
import numpy as np

def dn_ddp(dp, n, gm, gsd):
    """Evaluate the number distribution as a lognormal PDF.

    The PDF of a lognormal distribution is calculated using equation 8.34
    from Seinfeld and Pandis.

    .. math::

        n_N(D_p)=\\frac{dN}{dD_p}=\\frac{N_t}{\sqrt{2π}D_p lnσ_g}exp\Big(-\\frac{(lnD_p - lnD̄_{pg})^2}{2ln^2σ_g}\Big)

    Parameters
    ----------
    dp : float or an array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    n(Dp) | float or an array of floats
        Returns the total number of particles at diameter dp in units of
        :math:`\mu m^{-1} cm^{-3}`.

    See Also
    --------
    opcsim.equations.pdf.dn_dlndp
    opcsim.equations.pdf.dn_dlogdp
    opcsim.equations.cdf.nt

    """
    res = (n / (np.sqrt(2*np.pi)*dp*np.log(gsd))) * \
            np.exp(-(np.log(dp) - np.log(gm))** 2 / (2*np.log(gsd)** 2))

    return res

def dn_dlndp(dp, n, gm, gsd):
    """The PDF of a lognormal distribution as calculated using equation 8.33 by way
    of 8.21 from Seinfeld and Pandis.

    .. math::

        n_N^e(lnD_p)=\\frac{dN}{dlnD_p}=\\frac{N_t}{\sqrt{2π} lnσ_g}exp(-\\frac{(lnD_p - lnD̄_{pg})^2}{2ln^2σ_g})

    Parameters
    ----------
    dp : float or an array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    dn/dlndp | float
        Returns the total number of particles at diameter dp in units of
        [cm-3]

    See Also
    --------
    opcsim.equations.pdf.dn_ddp
    opcsim.equations.pdf.dn_dlogdp
    opcsim.equations.cdf.nt

    """
    return dp * dn_ddp(dp, n, gm, gsd)

def dn_dlogdp(dp, n, gm, gsd):
    """The PDF of a lognormal distribution on a log10 basis as calculated using
    equation 8.18 from Seinfeld and Pandis.

    .. math::

        n_N^o(logD_p)=\\frac{dN}{dlogD_p}=ln(10)D_pn_N(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    dn/dlogdp | float
        Returns the total number of particles at diameter dp in units of
        [cm-3]

    See Also
    --------
    opcsim.equations.pdf.dn_ddp
    opcsim.equations.pdf.dn_dlndp
    opcsim.equations.cdf.nt

    """
    return np.log(10) * dp * dn_ddp(dp, n, gm, gsd)

def ds_ddp(dp, n, gm, gsd):
    """The surface-area weighted PDF of a lognormal distribution as calculated
    using equation 8.4 from Seinfeld and Pandis.

    .. math::

        n_s(D_p)=πD_p^2n_N(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    ds/ddp | float
        Returns the surface area of particles at diameter dp in units of
        [um cm-3]

    See Also
    --------
    opcsim.equations.pdf.ds_dlndp
    opcsim.equations.pdf.ds_dlogdp
    opcsim.equations.cdf.st

    """
    return np.pi * dp ** 2 * dn_ddp(dp, n, gm, gsd)

def ds_dlndp(dp, n, gm, gsd):
    """The surface-area weighted PDF of a lognormal distribution as calculated
    using equation 8.10 from Seinfeld and Pandis.

    .. math::

        n_s^e(ln D_p)=πD_p^2n_N^e(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    ds/dlndp | float
        Returns the surface area of particles at diameter dp on a log basis
        in units of [um2 cm-3]

    See Also
    --------
    opcsim.equations.pdf.ds_ddp
    opcsim.equations.pdf.ds_dlogdp
    opcsim.equations.cdf.st

    """
    return np.pi * dp ** 2 * dn_dlndp(dp, n, gm, gsd)

def ds_dlogdp(dp, n, gm, gsd):
    """The surface-area weighted PDF of a lognormal distribution as calculated
    using equation 8.19 from Seinfeld and Pandis.

    .. math::

        n_s^o(log D_p)=ln(10)D_p n_s(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    ds/dlogdp | float
        Returns the surface area of particles at diameter dp on a log10 basis
        in units of [um2 cm-3]

    See Also
    --------
    opcsim.equations.pdf.ds_ddp
    opcsim.equations.pdf.ds_dlndp
    opcsim.equations.cdf.st

    """
    return np.log(10) * dp * ds_ddp(dp, n, gm, gsd)

def dv_ddp(dp, n, gm, gsd):
    """The volume weighted PDF of a lognormal distribution as calculated
    using equation 8.6 from Seinfeld and Pandis.

    .. math::

        n_V(D_p)=\\frac{π}{6}D_p^3 n_N(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    dv/ddp | float
        Returns the volume of particles at diameter dp
        in units of [um2 cm-3]

    See Also
    --------
    opcsim.equations.pdf.dv_dlndp
    opcsim.equations.pdf.dv_dlogdp
    opcsim.equations.cdf.vt

    """
    return (np.pi / 6.) * dp ** 3 * dn_ddp(dp, n, gm, gsd)

def dv_dlndp(dp, n, gm, gsd):
    """The volume weighted PDF of a lognormal distribution as calculated
    using equation 8.10 from Seinfeld and Pandis.

    .. math::

        n_V^e(ln D_p)=\\frac{π}{6}D_p^3 n_N^e(ln D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    dV/dlndp | float
        Returns the volume of particles at diameter dp on a log basis
        in units of [um3 cm-3]

    See Also
    --------
    opcsim.equations.pdf.dv_ddp
    opcsim.equations.pdf.dv_dlogdp
    opcsim.equations.cdf.vt

    """
    return (np.pi / 6.) * dp ** 3 * dn_dlndp(dp, n, gm, gsd)

def dv_dlogdp(dp, n, gm, gsd):
    """The volume weighted PDF of a lognormal distribution as calculated
    using equation 8.20 from Seinfeld and Pandis.

    .. math::

        n_V^o(log D_p)=log(10)*D_p n_V(D_p)

    Parameters
    ----------
    dp : float or array of floats
        Particle diameter in microns.
    n : float
        Total aerosol number concentration in units of #/cc
    gm : float
        Median particle diameter (geometric mean) in units of microns.
    gsd : float
        Geometric Standard Deviation of the distribution.

    Returns
    -------
    dV/dlogdp | float
        Returns the volume of particles at diameter dp on a log10 basis
        in units of [um3 cm-3]

    See Also
    --------
    opcsim.equations.pdf.dv_ddp
    opcsim.equations.pdf.dv_dlndp
    opcsim.equations.cdf.vt

    """
    return np.log(10) * dp * dv_ddp(dp, n, gm, gsd)
