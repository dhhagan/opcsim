# -*- coding: utf-8 -*-
"""
"""
from .equations.pdf import *
from .equations.cdf import *

DISTRIBUTION_DATA = {
    'Urban': [
        (7100, 0.0117, 0.232, "Mode I"),
        (6320, 0.0373, 0.250, "Mode II"),
        (960, 0.151, 0.204, "Mode III")
        ],
    'Marine': [
        (133, 0.008, 0.657, "Mode I"),
        (66.6, 0.266, 0.210, "Mode II"),
        (3.1, 0.58, 0.396, "Mode III")
        ],
    'Rural': [
        (6650, 0.015, 0.225, "Mode I"),
        (147, 0.054, 0.557, "Mode II"),
        (1990, 0.084, 0.266, "Mode III")
        ],
    'Remote Continental': [
        (3200, 0.02, 0.161, "Mode I"),
        (2900, 0.116, 0.217, "Mode II"),
        (0.3, 1.8, 0.380, "Mode III")
        ],
    'Free Troposphere': [
        (129, 0.007, 0.645, "Mode I"),
        (59.7, 0.250, 0.253, "Mode II"),
        (63.5, 0.52, 0.425, "Mode III")
        ],
    'Polar': [
        (21.7, 0.138, 0.245, "Mode I"),
        (0.186, 0.75, 0.300, "Mode II"),
        (3e-4, 8.6, 0.291, "Mode III")
        ],
    'Desert': [
        (726, 0.002, 0.247, "Mode I"),
        (114, 0.038, 0.770, "Mode II"),
        (0.178, 21.6, 0.438, "Mode III")
        ],
}

def _get_pdf_func(base, weight, D, N, GM, GSD):
    """
    """
    base    = base.lower()
    weight  = weight.lower()

    if weight not in ['number', 'surface', 'volume']:
        raise Exception("Invalid argument for weight: ['number', 'surface', 'volume']")

    if base not in [None, 'none', 'log', 'log10']:
        raise Exception("Invalid argument for base: ['none', 'log', 'log10']")

    if base == 'none' or base == None:
        if weight == 'number':
            return dn_ddp(D, N, GM, GSD)
        elif weight == 'surface':
            return ds_ddp(D, N, GM, GSD)
        elif weight == 'volume':
            return dv_ddp(D, N, GM, GSD)
    elif base == 'log':
        if weight == 'number':
            return dn_dlndp(D, N, GM, GSD)
        elif weight == 'surface':
            return ds_dlndp(D, N, GM, GSD)
        elif weight == 'volume':
            return dv_dlndp(D, N, GM, GSD)
    elif base == 'log10':
        if weight == 'number':
            return dn_dlogdp(D, N, GM, GSD)
        elif weight == 'surface':
            return ds_dlogdp(D, N, GM, GSD)
        elif weight == 'volume':
            return dv_dlogdp(D, N, GM, GSD)

def _get_cdf_func(weight, D, N, GM, GSD):
    """
    """
    weight = weight.lower()
    if weight not in ['number', 'surface', 'volume']:
        raise Exception("Invalid argument for weight: ['number', 'surface', 'volume']")

    if weight == 'number':
        return Nt(D, N, GM, GSD)
    elif weight == 'surface':
        return St(D, N, GM, GSD)
    elif weight == 'volume':
        return Vt(D, N, GM, GSD)

def load_distribution(label):
    """Load sample distributions as described by S+P Table 8.3.

    Parameters
    ----------

    label : {'Urban' | 'Marine' | 'Rural' | 'Remote Continental' | 'Free Troposphere' | 'Polar' | 'Desert'}
        Choose which sample distribution to load.

    Returns
    -------
    Instance of the AerosolDistribution class

    Examples
    --------

    >>> d = opcsim.load_distribution("Urban")

    """
    if label not in DISTRIBUTION_DATA.keys():
        raise ValueError("Invalid label.")

    _tmp = AerosolDistribution(label)

    for each in DISTRIBUTION_DATA[label]:
        _tmp.add_mode(each[0], each[1], 10**each[2], each[3])

    return _tmp

class AerosolDistribution(object):
    """Define an aerosol distribution.

    We are operating on the assumption that any aerosol distribution can be
    described as the sum of n lognormal modes following Seinfeld and Pandis
    equation 8.54.

    .. math::

        n_N^o(logD_p)=Σ_{i=1}^{n}\\frac{N_i}{\sqrt{2π} logσ_i}exp(-\\frac{(logD_p - logD̄_p)^2}{2log^2σ_i})

    """
    def __init__(self, label=None):
        """Initialize the Aerosol Distribution.

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

    def add_mode(self, n, gm, gsd, label=None):
        """Add a mode to the distribution.

        Parameters
        ----------
        n : float
            Total number of particles (#/cc)
        gm : float
            Geometric Mean particle diameter (um)
        gsd : float
            Geometric Standard Deviation
        label : string, optional
            Label for the mode

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
        self.modes.append({'label':label, 'N':n, 'GM':gm, 'GSD':gsd})

        return

    def pdf(self, dp, base='log10', weight='number', mode=None):
        """Evaluate and return the probability distribution function at
        particle diameter `dp`.

        Using equation 8.54 from Seinfeld and Pandis, we can evaluate the
        probability distribution function for the multi-modal aerosol
        distribution by summing the individual pdfs.

        Parameters
        ----------
        dp : float or array of floats
            Particle diameter(s) to evaluate the pdf (um)
        base : {None | 'none' | 'log' | 'log10'}
            Base algorithm to use. Default is 'log10'
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`.
        mode : string or None
            Choose to only evaluate the pdf for a single mode
            of the entire distribution. If set to `None`, the entire
            distribution will be evaluated.

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
            value += _get_pdf_func(base, weight, dp, mode['N'], mode['GM'], mode['GSD'])

        return value

    def cdf(self, dmin=0.0, dmax=None, weight='number', mode=None):
        """Evaluate and return the cumulative probability distribution function
        between dmin and dmax.

        Using equation _ from Seinfeld and Pandis, we can evaluate the cdf of a
        a multi-modal particle size distribution by summing the individual
        distributions. For example, evaluating the cdf over the entire size
        range will return the total number of particles in that size range. If
        weighted by surface area or volume, it will return the integrated
        surface area or volume respectively.

        Mathematically, it is represented as:

        .. math::

            Σ

        Parameters
        ----------
        dmin : float
            The minimum particle diameter in the integration (um)
        dmax : float
            The maximum particle diameter in the integration (um)
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`
        mode : string or None
            Choose to only evaluate the pdf for a single mode
            of the entire distribution. If set to `None`, the entire
            distribution will be evaluated.

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
        if not dmin < dmax:
            raise ValueError("dmin must be less than dmax")

        value = 0.0

        if mode is not None:
            modes = [self._get_mode(mode)]
        else:
            modes = self.modes

        for mode in modes:
            value += _get_cdf_func(weight, dmax, mode['N'], mode['GM'], mode['GSD'])

            if dmin > 0.0:
                value -= _get_cdf_func(weight, dmin, mode['N'], mode['GM'], mode['GSD'])

        return value

    def mean(self, mode = None, weight = 'number', diameter = True):
        """
            Calculate the mean diameter, surface area, or volume weighted by
            number, surface area, or volume.

            weight = 'number',  diameter = True     Mean Diameter               S+P 8.44
            weight = 'number',  diameter = False    Mean Diameter               S+P 8.44
            weight = 'surface', diameter = True     Surface Area Mean Diameter  S+P T8.2
            weight = 'surface', diameter = False    Mean Surface Area           S+P T8.2
            weight = 'volume',  diameter = True     Volume Mean Diameter        S+P T8.2
            weight = 'volume',  diameter = False    Mean Volume                 S+P T8.2
        """
        # Get the mode
        if len(self.modes) > 1:
            _m_ = self._get_mode(mode)
        else:
            _m_ = self.modes[0]

        if weight == 'number':
            res = _m_['GM'] * np.exp( np.log( _m_['GSD'] ) ** 2 / 2. )

        elif weight == 'surface':
            res = (1. / _m_['N']) * self.cdf(1000., weight = 'surface', mode = _m_['label'])

            if diameter == True:
                res = np.sqrt( res / np.pi )

        elif weight == 'volume':
            res = (1. / _m_['N']) * self.cdf(1000., weight = 'volume', mode = _m_['label'])

            if diameter == True:
                res = ( (6. / np.pi) * res ) ** ( 1. / 3. )

        else:
            raise Exception("Invalid weight: ['number', 'surface', 'volume']")

        return res

    def median(self, mode = None, weight = 'number'):
        """
            Calculate the weighted median diameter

            weight = 'number'   Median Particle Diameter        argument
            weight = 'surface'  Surface Area Median Diameter    S+P 8.50
            weight = 'volume'   Volume Median Diameter          S+P 8.53

        """
        # Get the mode
        if len(self.modes) > 1:
            _m_ = self._get_mode(mode)
        else:
            _m_ = self.modes[0]

        if weight == 'number':
            res = _m_['GM']
        elif weight == 'surface':
            res = np.exp( np.log(_m_['GM']) + 2 * np.log( _m_['GSD'] ) ** 2 )
        elif weight == 'volume':
            res = np.exp( np.log(_m_['GM']) + 3 * np.log( _m_['GSD'] ) ** 2 )
        else:
            raise Exception("Invalid weight: ['number', 'surface', 'volume']")

        return res

    def __repr__(self):
        res = "AerosolDistribution: {}\n\n".format(self.label)
        res = res + "Mode\tN\t\tGM\tGSD\n"

        for mode in self.modes:
            res = res + "{}\t{:.2e}\t{:.3f}\t{:.3f}\n".format(mode['label'], mode['N'], mode['GM'], mode['GSD'])

        return res
