# -*- coding: utf-8 -*-
"""
"""
from .equations.pdf import *
from .equations.cdf import *
from .utils import k_kohler, rho_eff, RHO_H20

DISTRIBUTION_DATA = {
    'urban': [
        (7100, 0.0117, 0.232, "Mode I"),
        (6320, 0.0373, 0.250, "Mode II"),
        (960, 0.151, 0.204, "Mode III")
        ],
    'marine': [
        (133, 0.008, 0.657, "Mode I"),
        (66.6, 0.266, 0.210, "Mode II"),
        (3.1, 0.58, 0.396, "Mode III")
        ],
    'rural': [
        (6650, 0.015, 0.225, "Mode I"),
        (147, 0.054, 0.557, "Mode II"),
        (1990, 0.084, 0.266, "Mode III")
        ],
    'remote continental': [
        (3200, 0.02, 0.161, "Mode I"),
        (2900, 0.116, 0.217, "Mode II"),
        (0.3, 1.8, 0.380, "Mode III")
        ],
    'free troposphere': [
        (129, 0.007, 0.645, "Mode I"),
        (59.7, 0.250, 0.253, "Mode II"),
        (63.5, 0.52, 0.425, "Mode III")
        ],
    'polar': [
        (21.7, 0.138, 0.245, "Mode I"),
        (0.186, 0.75, 0.300, "Mode II"),
        (3e-4, 8.6, 0.291, "Mode III")
        ],
    'desert': [
        (726, 0.002, 0.247, "Mode I"),
        (114, 0.038, 0.770, "Mode II"),
        (0.178, 21.6, 0.438, "Mode III")
        ],
}

def _get_pdf_func(base, weight, dp, n, gm, gsd, rho=1.):
    """"""
    weight  = weight.lower()

    if weight not in ['number', 'surface', 'volume', 'mass']:
        raise Exception("Invalid argument for weight: ['number', 'surface', 'volume', 'mass']")

    if base not in [None, 'none', 'log', 'log10']:
        raise Exception("Invalid argument for base: ['none', 'log', 'log10']")

    if base == 'none' or base == None:
        if weight == 'number':
            return dn_ddp(dp, n, gm, gsd)
        elif weight == 'surface':
            return ds_ddp(dp, n, gm, gsd)
        elif weight == 'volume':
            return dv_ddp(dp, n, gm, gsd)
        elif weight == 'mass':
            return dv_ddp(dp, n, gm, gsd) * rho * 1e-6
    elif base == 'log':
        if weight == 'number':
            return dn_dlndp(dp, n, gm, gsd)
        elif weight == 'surface':
            return ds_dlndp(dp, n, gm, gsd)
        elif weight == 'volume':
            return dv_dlndp(dp, n, gm, gsd)
        elif weight == 'mass':
            return dv_dlndp(dp, n, gm, gsd) * rho *1e-6
    elif base == 'log10':
        if weight == 'number':
            return dn_dlogdp(dp, n, gm, gsd)
        elif weight == 'surface':
            return ds_dlogdp(dp, n, gm, gsd)
        elif weight == 'volume':
            return dv_dlogdp(dp, n, gm, gsd)
        elif weight == 'mass':
            return dv_dlogdp(dp, n, gm, gsd) * rho * 1e-6


def _get_cdf_func(n, gm, gsd, dmin=None, dmax=10., weight='number', rho=1.):
    """"""
    weight = weight.lower()
    if weight not in ['number', 'surface', 'volume', 'mass']:
        raise Exception("Invalid argument for weight: ['number', 'surface', 'volume', 'mass']")

    if weight == 'number':
        return nt(n, gm, gsd, dmin, dmax)
    elif weight == 'surface':
        return st(n, gm, gsd, dmin, dmax)
    elif weight == 'volume':
        return vt(n, gm, gsd, dmin, dmax)
    elif weight == 'mass':
        return vt(n, gm, gsd, dmin, dmax) * rho


def load_distribution(label):
    """Load sample distributions as described by Seinfeld+Pandis Table 8.3.

    There are currently 7 options including: Urban, Marine, Rural, Remote
    continental, Free troposphere, Polar, and Desert.

    Parameters
    ----------

    label : {'Urban' | 'Marine' | 'Rural' | 'Remote Continental' | 'Free Troposphere' | 'Polar' | 'Desert'}
        Choose which sample distribution to load.

    Returns
    -------
    An instance of the AerosolDistribution class

    Examples
    --------

    >>> d = opcsim.load_distribution("Urban")

    """
    label = label.lower()

    if label not in DISTRIBUTION_DATA.keys():
        raise ValueError("Invalid label.")

    _tmp = AerosolDistribution(label)

    for each in DISTRIBUTION_DATA[label]:
        _tmp.add_mode(each[0], each[1], 10**each[2], each[3])

    return _tmp


class AerosolDistribution(object):
    """Define an aerosol distribution.

    Assuming an aerosol distribution can be described as the sum of n lognormal
    modes, define as (Seinfeld and Pandis equation 8.54):

    .. math::

        n_N^o(logD_p)=Σ_{i=1}^{n}\\frac{N_i}{\sqrt{2π} logσ_i}exp(-\\frac{(logD_p - logD̄_{pi})^2}{2log^2σ_i})

    """
    def __init__(self, label=None):
        """Initialize an Aerosol Distribution.

        Parameters
        ----------
        label : string, optional
            Label the distribution

        Returns
        -------
        AerosolDistribution
            An instance of the AerosolDistribution class

        Examples
        --------

        >>> d = AerosolDistribution("Urban")

        """
        self.label = label
        self.modes = []

    def _get_mode(self, label):
        """Return the mode with label=`label`. """
        for mode in self.modes:
            if mode['label'].lower() == label.lower():
                return mode

        return None

    def add_mode(self, n, gm, gsd, label=None, kappa=0., rho=1., refr=complex(1.5, 0)):
        """Add a mode to the distribution as defined using N, GM, and GSD. Additionally,
        each mode has optical and chemical properties that can be set including: the k-kohler 
        coefficient, the density, and the refractive index.

        Parameters
        ----------
        n : float
            Total number of particles (#/cc)
        gm : float
            Median particle diameter (Geometric Mean) in units of microns.
        gsd : float
            Geometric Standard Deviation
        label : string, optional
            Label for the mode
        kappa: float, optional
            The k-kohler coefficient to describe hygroscopic growth. Default is 0.
        rho: float: optional
            The particle density in units of g/cm3. Default is 1.
        refr: complex: optional
            The complex refractive index. Default is 1.5+0i. This should be the dry 
            refractive index. If/when the distribution is evaluated at a non-zero RH,
            the refractive index will be adjusted based on it's water content at that 
            RH using the defined kappa value.

        Examples
        --------

        Create a single-mode aerosol distribution with parameters N=7100,
        GM=11.7 nm, GSD=1.706. These are the sample parameters for the first
        mode of a typical Urban aerosol distribution as described in Seinfeld
        and Pandis (Table 8.3).

        >>> d = AerosolDistribution("Urban")
        >>> d.add_mode(n=7100, gm=0.0117, gsd=1.706, label="Mode 1")

        We can also go ahead and add the second and third modes if we so choose:

        >>> d.add_mode(n=6320, gm=0.0373, gsd=1.778, label="Mode 2")
        >>> d.add_mode(n=960, gm=0.151, gsd=1.599, label="Mode 3")

        """
        self.modes.append(
            {
                'label':label if label else "Mode {}".format(len(self.modes)), 
                'N':n, 'GM':gm, 
                'GSD':gsd, 
                "kappa": kappa, 
                "rho": rho, 
                "refr": refr
                }
            )

    def pdf(self, dp, base='log10', weight='number', mode=None, rh=0., rho=None):
        """Evaluate and return the probability distribution function at
        particle diameter `dp`.

        Using equation 8.54 from Seinfeld and Pandis, we can evaluate the
        probability distribution function for a multi-modal aerosol
        distribution by summing the individual pdf's. By default, the calculations
        are based on the dry diameter of the particle (assuming RH=0); however,
        you can adjust the relative humidity using the **rh** parameter which 
        will calculate the size based on growth per k-kohler theory using the 
        k-kohler coefficient set for each individual mode.

        Parameters
        ----------
        dp : float or an array of floats
            Particle diameter(s) to evaluate the pdf (um)
        base : {None | 'none' | 'log' | 'log10'}
            Base algorithm to use. Default is 'log10'
        weight : {'number' | 'surface' | 'volume' | 'mass'}
            Choose how to weight the pdf. Default is `number`.
        mode : string or None
            Choose to only evaluate the pdf for a single mode
            of the entire distribution. If set to `None`, the entire
            distribution will be evaluated.
        rh: float: optional
            The relative humidity as a percentage (0-100). Default is 0.
        rho: float: optional
            The particle density. If set, this will override the density set 
            for individual modes.

        Returns
        -------
        float
            The evaluated pdf at particle diameter dp

        Examples
        --------

        Evaluate the PDF at 100 nm

        >>> d = opcsim.load_distribution("Urban")
        >>> d.pdf(0.1)

        Evaluate the PDF at a range of particle diameters

        >>> d.pdf(np.linespace(0.1, 1., 100))

        Evaluate the PDF for a volume-weighted distribution -> returns dVdlogDp

        >>> d.pdf(0.1, weight='volume')

        Evaluate the PDF for a volume-weighted distribution in ln-space

        >>> d.pdf(0.1, weight='volume', base='log')

        """
        value = 0.0

        if mode is not None:
            modes = [self._get_mode(mode)]
        else:
            modes = self.modes

        for mode in modes:
            # calculate the wet diameter of the particle based on k-kohler theory
            gm = k_kohler(diam_dry=mode['GM'], kappa=mode["kappa"], rh=rh)

            # override rho if set 
            rho = rho if rho else mode["rho"]

            # calculate the effective density, taking into account hygroscopic growth
            rho = rho_eff([rho, RHO_H20], diams=[mode['GM'], gm - mode['GM']])

            value += _get_pdf_func(base, weight, dp, mode['N'], gm, mode['GSD'], rho)

        return value

    def cdf(self, dmax, dmin=None, weight='number', mode=None, rh=0., rho=None):
        """Evaluate and return the cumulative probability distribution function
        between `dmin` and `dmax`.

        Using equation _ from Seinfeld and Pandis, we can evaluate the cdf of a
        a multi-modal particle size distribution by summing the individual
        distributions. For example, evaluating the cdf over the entire size
        range will return the total number of particles in that size range. If
        weighted by surface area or volume, it will return the integrated
        surface area or volume respectively. By default, the calculations
        are based on the dry diameter of the particle (assuming RH=0); however,
        you can adjust the relative humidity using the **rh** parameter which 
        will calculate the size based on growth per k-kohler theory using the 
        k-kohler coefficient set for each individual mode.

        Parameters
        ----------
        dmin : float
            The minimum particle diameter in the integration (um)
        dmax : float
            The maximum particle diameter in the integration (um)
        weight : {'number' | 'surface' | 'volume' | 'mass'}
            Choose how to weight the pdf. Default is `number`
        mode : string or None
            Choose to only evaluate the pdf for a single mode
            of the entire distribution. If set to `None`, the entire
            distribution will be evaluated.
        rh: float: optional
            The relative humidity as a percentage (0-100). Default is 0.
        rho: float: optional
            The particle density. If set, this will override the density set 
            for individual modes.

        Returns
        -------
        float
            The integrated distribution function representing the total number of
            {particles, surface area, volume} between dmin and dmax.

        Examples
        --------

        Evaluate the CDF between over the entire distribution:

        >>> d = opcsim.load_distribution("Urban")
        >>> d.cdf()

        Evaluate the CDF up to 100 nm

        >>> d.cdf(dmax=0.1)

        Evaluate the total volume/mass under 2.5 microns

        >>> d.cdf(dmax=2.5, weight='volume')

        """
        if dmin is not None:
            if not dmin < dmax:
                raise ValueError("dmin must be less than dmax")

        value = 0.0

        if mode is not None:
            modes = [self._get_mode(mode)]
        else:
            modes = self.modes

        for mode in modes:
            # calculate the wet diameter of the particle based on k-kohler theory
            gm = k_kohler(diam_dry=mode['GM'], kappa=mode["kappa"], rh=rh)

            # override rho if set
            rho = rho if rho else mode["rho"]

            # calculate the effective density, taking into account hygroscopic growth
            rho = rho_eff([rho, RHO_H20], diams=[mode['GM'], gm - mode['GM']])

            value += _get_cdf_func(mode['N'], gm, mode['GSD'], dmin, dmax, weight, rho)

        return value

    def __repr__(self):
        return "AerosolDistribution: {}".format(self.label)
