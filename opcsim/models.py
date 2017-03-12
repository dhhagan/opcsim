import numpy as np
import pandas as pd
import math

from .distributions import AerosolDistribution
from .utils import make_bins

def constant(dp, x=1.):
    """Define an efficiency curve for the OPC. """
    return dp * 0.0 + x

class OPC(object):
    """Simulate an Optical Particle Counter (OPC) as defined by key instrument
    variables.
    """
    def __init__(self, n_bins=1, dmin=0.5, dmax=2.5, ce=constant,
                    bins=None, **kwargs):
        """Initialize the simulated OPC.

        Parameters
        ----------
        n_bins : int, optional
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
        OPC : class instance
            Instance of the OPC class.

        Examples
        --------

        Initialize an OPC with 1 bin (Sharp, Shinyei, etc.)

        >>> opc = opcsim.OPC(n_bins=1)

        Initialize an OPC with 2 bins (0.5-2.5 um, 2.5-10um)

        >>> import numpy as np
        >>> bins = np.array([[0.5, 2.5], [2.5, 10]])
        >>> opc = opcsim.OPC(bins=bins)

        """
        self.n_bins = n_bins
        self.dmin = dmin
        self.dmax = dmax
        self.bins = bins
        self.ce = ce

        # Make sure the bins are set
        if bins is None:
            self.bins = make_bins(self.dmin, self.dmax, self.n_bins)
        else:
            if self.bins.shape[1] != 3:
                raise Exception("Bins must be a 3xn array.")

            # Set the number of bins based on the array
            self.n_bins = self.bins.shape[0]
            self.dmin = self.bins[0, 0]
            self.dmax = self.bins[-1, -1]

    @property
    def dlogdp(self):
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def midpoints(self):
        return self.bins[:, 1]

    def evaluate(self, distribution, weight='number', base='log10'):
        """Evaluate the OPC given the PDF of an aerosol distribution.

        In the future, we should account for intra-bin counting efficiency
        errors.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`.
        base : {'none' | 'log' | 'log10'}
            Base algorithm to use. Default is 'log10'.

        Returns
        -------
        d[weight]/d[base]Dp : array
            An array containing the evaluated AerosolDistribution PDF at all
            midpoint particle diameters.

        See Also
        --------
        opcsim.AerosolDistribution

        Examples
        --------

        Evaluate the Urban Aerosol Distribution using a 2-bin OPC:

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> dNdlogDp = opc.evaluate(d)

        Evaluate the same OPC, but weight by volume/mass:

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> dVdlogDp = opc.evaluate(d, weight='volume')

        """
        if not isinstance(distribution, AerosolDistribution):
            raise Exception("Invalid AerosolDistribution")

        return self.ce(self.midpoints) * distribution.pdf(self.midpoints, base=base, weight=weight)

    def histogram(self, distribution, weight='number', base='log10'):
        """Return a histogram corresponding to how an OPC "sees" an aerosol
        distribution.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`.
        base : {'none' | 'log' | 'log10'}
            Base algorithm to use. Default is 'log10'.

        Returns
        -------
        left_bin, d[weight]/d[base]Dp, bin_width : arrays
            Returns three arrays containing the left bin boundary, the evaluated
            PDF, and the width of the bin (dDp). This data can be
            directly plotted as a histogram using matplotlib bar plots.


        Examples
        --------

        Evaluate the Urban Aerosol Distribution using a 2-bin OPC:

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> bb, h, bw = opc.histogram(d)

        Evaluate the same OPC, but weight by volume/mass:

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> bb, h, bw = opc.histogram(d, weight='volume')

        """
        bb = self.bins[:, 0]
        dDp = self.bins[:, 2] - self.bins[:, 0]
        pdf = self.evaluate(distribution, weight, base)

        return bb, pdf, dDp

    def number(self, distribution, measured=True):
        """Return the total number of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        measured: bool
            If true, the result returns the value as "seen" by the OPC. If false,
            the true integrated CDF for each bin is returned.

        Returns
        -------
        number : array of floats
            Returns the total number of particles seen in each bin. Units are #/cc.

        See Also
        --------

        opcsim.equations.cdf.nt


        Examples
        --------

        Determine the number of particles a basic OPC sees in the Urban
        Distribution

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> nt = opc.number(d)

        What if we want to know how many particles were actually in each bin by
        integrating the underlying distribution?

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> nt = opc.number(d, measured=False)

        """
        if measured == True:
            vals = self.evaluate(distribution) * self.dlogdp
        else:
            vals = [distribution.cdf(dmin=self.bins[i, 0],
                    dmax=self.bins[i, -1]) for i in range(self.n_bins)]

        return vals

    def surface_area(self, distribution, measured=True):
        """Return the total surface area of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        measured: bool
            If true, the result returns the value as "seen" by the OPC. If false,
            the true integrated CDF for each bin is returned.

        Returns
        -------
        surface_area : array of floats
            Returns the total surface area of particles seen in each bin.
            Units are um2 cm-3.

        See Also
        --------

        opcsim.equations.cdf.st


        Examples
        --------

        Determine the surface area of particles a basic OPC sees in the Urban
        Distribution

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> st = opc.surface_area(d)

        What if we want to know how many particles were actually in each bin by
        integrating the underlying distribution?

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> st = opc.surface_area(d, measured=False)

        """
        if measured == True:
            vals = self.evaluate(distribution, weight='surface') * self.dlogdp
        else:
            vals = [distribution.cdf(dmin=self.bins[i, 0],
                    dmax=self.bins[i, -1], weight='surface') for i in range(self.n_bins)]

        return vals

    def volume(self, distribution, measured=True):
        """Return the total volume of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        measured: bool
            If true, the result returns the value as "seen" by the OPC. If false,
            the true integrated CDF for each bin is returned.

        Returns
        -------
        surface_area : array of floats
            Returns the total volume of particles seen in each bin.
            Units are um3 cm-3.

        See Also
        --------

        opcsim.equations.cdf.vt


        Examples
        --------

        Determine the volume of particles a basic OPC sees in the Urban
        Distribution

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> vt = opc.volume(d)

        What if we want to know the total volume of particles that were
        actually in each bin by integrating the underlying distribution?

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> vt = opc.volume(d, measured=False)

        What if we want to know the total volume across all bins?

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> vt = opc.volume(d).sum()

        """
        if measured == True:
            vals = self.evaluate(distribution, weight='volume') * self.dlogdp
        else:
            vals = [distribution.cdf(dmin=self.bins[i, 0],
                    dmax=self.bins[i, -1], weight='volume') for i in range(self.n_bins)]

        return vals

__all__ = [
    'OPC'
]
