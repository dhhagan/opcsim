"""
"""
from .equations.pdf import *
from .equations.cdf import *

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

class AerosolDistribution(object):
    """Creates an aerosol distribution as represented by the sum of n lognormal
    modes.
    """
    def __init__(self, label = None, **kwargs):
        self.label = label
        self.modes = []

    def _get_mode(self, label):
        """Return the mode which contains label = label
        """
        for mode in self.modes:
            if mode['label'].lower() == label.lower():
                return mode

        return None

    def add_mode(self, params, label = None):
        """Add an additional mode to the AerosolDistribution by providing the
        parameters as a tuple: (N, GM, GSD)
        """

        N, GM, GSD = params

        self.modes.append({ 'label': label, 'N': N, 'GM': GM, 'GSD': GSD })

        return

    def pdf(self, Dp, base = 'none', weight = 'number'):
        """
            Return the probability distribution function (PDF) as a sum of the
            individual PDF's.

            base: none, log/ln, log10

            Eq. 8.54 S+P
        """
        total = 0.0

        for mode in self.modes:
            total = total + _get_pdf_func(base, weight, Dp, mode['N'], mode['GM'], mode['GSD'])

        return total

    def cdf(self, Dmax, weight = 'number', Dmin = 0.0, mode = None):
        """
            Return the cumulative distribution function (CDF) as a sum of all
            individual CDF's.
        """
        assert (Dmin < Dmax), "Dmin must be less than Dmax"

        total = 0.0

        if mode is not None:
            modes = [self._get_mode(mode)]
        else:
            modes = self.modes

        for mode in modes:
            total = total + _get_cdf_func(weight, Dmax, mode['N'], mode['GM'], mode['GSD'])

            if Dmin > 0.0:
                total = total - _get_cdf_func(weight, Dmin, mode['N'], mode['GM'], mode['GSD'])

        return total

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
