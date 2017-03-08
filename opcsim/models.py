import numpy as np
import pandas as pd
import math

from .distributions import AerosolDistribution

_avail_params = [
    'nm',
    'na',
    'vm',
    'va',
    'nm_total',
    'na_total',
    'vm_total',
    'va_total',
    'nm_na',
    'nm_va',
    'vm_va'
]

def _set_bins(dmin, dmax, num_bins):
    """
        Create and return a list of bins that are equally spaced (in log space)
        between dmin and dmax assuming num_bins bins.
    """
    bins    = np.zeros(( num_bins, 3 )) * np.nan

    bins[0, 0]      = dmin
    bins[-1, 2]     = dmax

    mult = 1. / (( np.log10(dmax) - np.log10(dmin) ) / num_bins )

    for i in range(num_bins):
        bins[i, 2]  = math.pow(10, np.log10( bins[i, 0]) + ( 1. / mult ))
        bins[i, 1]  = math.pow(10, np.mean([np.log10(bins[i, 0]), np.log10(bins[i, 2])]))

        if i < num_bins - 1:
            bins[i + 1, 0]  = bins[i, 2]

    return bins

def constant(x):
    """Define an efficiency curve for the OPC. """
    return x * 0.0 + 1.

class OPC(object):
    """Simulate an Optical Particle Counter (OPC) as defined by key instrument
    variables.
    """
    def __init__(self, num_bins=1, dmin=0.5, dmax=2.5, ce=constant, bins=None, **kwargs):
        """Initialize the simulated OPC.

        Parameters
        ----------
        num_bins : int, optional
            The number of discrete size bins for the OPC
        dmin : float, optional
            Minimum particle diameter the OPC can "see" in units of microns.
        dmax : float, optional
            Maximum particle diameter the OPC can "see" in units of microns.
        ce : callable, optional
            Function to be used to calculate the counting efficiency. Should map
            bin diameters (Dp) to a single value.
        bins : 3xn array, optional
            Array of bin diameters. If not set, the bin diameters will be
            automatically set by assuming equalivant bins in logspace between
            dmin and dmax.
        kwargs : key, value pairings
            Additional keyword arguments.

        Returns
        -------
        Instance of the OPC class.

        Examples
        --------

        Initialize an OPC with 1 bin (Sharp, Shinyei, etc.)

        >>> opc = opcsim.OPC(n_bins=1)

        Initialize an OPC with 2 bins (0.5-2.5 um, 2.5-10um)

        >>> import numpy as np
        >>> bins = np.array([])
        >>> opc = opcsim.OPC(n_bins=2, bins=bins)

        """
        self.num_bins = num_bins
        self.dmin = dmin
        self.dmax = dmax

        if bins is None:
            self.bins = _set_bins(dmin, dmax, num_bins)
        else:
            # Make sure that the bins are in the desired format
            assert (bins.shape[1] == 2 or bins.shape[1] == 3), "Bins must be a 3xn array"

            # If the shape is 2, then add in the middle bin using the log mean
            if bins.shape[1] == 2:
                self.bins = np.zeros(( bins.shape[0], 3 )) * np.nan

                self.bins[:, 0] = bins[:, 0]
                self.bins[:, 2] = bins[:, -1]

                for i in range(bins.shape[0]):
                    self.bins[i, 1]  = math.pow(10, np.mean([np.log10(self.bins[i, 0]), np.log10(self.bins[i, 2])]))
            else:
                self.bins       = bins

            self.num_bins   = self.bins.shape[0]
            self.dmin       = self.bins[0, 0]
            self.dmax       = self.bins[-1, -1]

        try:
            self.ce         = ce(self.bins[:, 1])
        except:
            self.ce         = 1.

        self.dlogdp     = np.log10(self.bins[:, 2]) - np.log10(self.bins[:, 0])

    def histogram(self, distribution, weight='number', base='log10'):
        """Return the histogram
            Accept a valid AerosolDistribution and return the histogram for this
            OPC based on its parameters and counting efficiency.

            Add in support for CE soon!

        """
        assert isinstance(distribution, AerosolDistribution), "Invalid AerosolDistribution"

        return distribution.pdf(self.bins[:, 1] , base = base, weight = weight)

    def boxes(self, distribution, weight = 'number', base = 'log10'):
        """
            Return the bins as a 3xn array for the given distribution. Should
            be in a format that can be easily plotted as a bar chart (left, height, width)

            `left` is the left hand side of bins
            `width` is dlogdp
            `height` is the histogram evaluated at the correct spots
        """
        left    = self.bins[:, 0]
        width   = self.bins[:, 2] - self.bins[:, 0]
        height  = self.histogram(distribution, weight = weight, base = base)

        return left, height, width

    def evaluate(self, distribution, param = 'nm', **kwargs):
        """
            Evaluate an individual parameter
        """
        dmax = kwargs.pop('dmax', 100.)

        data = self.histogram(distribution, weight = 'number', base = 'log10')

        if param == 'nm':
            res = data * self.dlogdp
        elif param == 'na':
            res = [distribution.cdf(dmin = self.bins[i, 0], dmax = self.bins[i, 2]) for i in range(self.num_bins)]
        elif param == 'vm':
            res = data * self.dlogdp * ( self.bins[:, 1] ** 3 * np.pi / 6. )
        elif param == 'va':
            res = [distribution.cdf(dmin = self.bins[i, 0], dmax = self.bins[i, 2], weight = 'volume') for i in range(self.num_bins)]
        elif param == 'nm_total':
            res = np.sum(self.evaluate(distribution, param = 'nm'))
        elif param == 'na_total':
            res = distribution.cdf(dmax = dmax, weight = 'number')
        elif param == 'vm_total':
            res = np.sum(self.evaluate(distribution, param = 'vm'))
        elif param == 'va_total':
            res = distribution.cdf(dmax = dmax, weight = 'volume')
        elif param == 'nm_na':
            res = round(self.evaluate(distribution, param = 'nm_total') / self.evaluate(distribution, param = 'na_total'), 3)
        elif param == 'nm_va':
            res = round(self.evaluate(distribution, param = 'nm_total') / self.evaluate(distribution, param = 'va_total'), 3)
        elif param == 'vm_va':
            res = round(self.evaluate(distribution, param = 'vm_total') / self.evaluate(distribution, param = 'va_total'), 3)
        else:
            raise Exception("Invalid param: [{}]".format(_avail_params))

        return res

    def _set_bins(self):
        """Calculate bin midpoints."""
        return

__all__ = [
    'OPC'
]
