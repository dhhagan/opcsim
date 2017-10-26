import numpy as np
import pandas as pd
import math

from .distributions import AerosolDistribution
from .utils import make_bins

def constant(dp, x=1.):
    """Define an efficiency curve for the OPC. """
    return dp * 0.0 + x

class OPC(object):
    """Define an Optical Particle Counter (OPC) using key instrument parameters.

    Key instrument parameters are defined as: the number of discrete size bins,
     the minimum particle diameter, the maximum particle diameter, and optionaly,
     the counting efficiency as a function of particle diameter.
    """
    def __init__(self, n_bins=1, dmin=0.5, dmax=2.5, ce=constant,
                    bins=None, **kwargs):
        """Initialize the simulated OPC.

        Parameters
        ----------
        n_bins : int, optional
            The number of discrete size bins
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
            An instance of the OPC class.

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
    def dlndp(self):
        return np.log(self.bins[:, -1] - np.log(self.bins[:, 0]))

    @property
    def ddp(self):
        return self.bins[:, -1] - self.bins[:, 0]

    @property
    def midpoints(self):
        return self.bins[:, 1]

    def evaluate(self, distribution, weight='number', base='log10', **kwargs):
        """Evaluate the OPC given the PDF of an aerosol distribution.

        Two methods of evaluation: 'simple' and 'subint'. The 'simple' method
        evaluates the PDF at each given bin midpoint, setting this value to be
        the number of particles the OPC senses in that bin.

        The 'subint' method performs an integration across each bin, taking into
        account the change in counting efficiency within the bin. This becomes
        much more important with wider bins, where the PDF can span orders of
        magnitude.

        The default is to perform the integration, as we feel this is the most
        accurate and honest way to evaluate what an OPC truly 'sees' within a
        given size bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.
        weight : {'number' | 'surface' | 'volume' | 'mass'}
            Choose how to weight the pdf. Default is `number`.
        base : {'none' | 'log' | 'log10'}
            Base algorithm to use. Default is 'log10'.
        method : {'subint' | 'simple'}
            Details the mechanism used for integrating the PDF/CDF. Default is
            'subint'.

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

        iters = kwargs.pop('n_iter', 50)
        method = kwargs.pop('method', 'subint')
        rho = kwargs.pop('rho', 1.)

        if method == 'simple':
            res = self.ce(self.midpoints) * distribution.pdf(
                                                self.midpoints,
                                                base=base,
                                                weight=weight)
        elif method == 'subint':
            # Iterate over each mode in the distribution
            # return an array of values (one for each bin)
            # account for number, surface area, volume, mass
            res = []

            for b in self.bins:
                _bins = make_bins(b[0], b[-1], iters)

                # Calculate CE at each subbin midpoint
                _eff = self.ce(_bins[:, 1])

                # Calculate the integrated total [W] within each bin
                # Always calculate number and then convert to SA or V if needed
                _vals = np.array([
                    distribution.cdf(
                                    dmin=_bins[k,0],
                                    dmax=_bins[k,-1],
                                    weight='number') for k in range(iters)])

                res.append((_vals*_eff).sum())

            res = np.array(res)

            # We calculate N here because that's what the OPC sees
            if weight == 'surface':
                res = (np.pi * self.midpoints**2) * res
            elif weight == 'volume':
                res = (self.midpoints**3*np.pi/6.) * res # dp3 * pi/6*N
            elif weight == 'mass':
                pass

            # Divide by appropriate width of bins
            if base == 'log':
                res = res / self.dlndp
            elif base == 'log10':
                res = res / self.dlogdp
            else:
                res = res / self.ddp
        else:
            raise Exception("Invalid method: {'subint' | 'simple'}")

        return res

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

    def number(self, distribution, **kwargs):
        """Return the total number of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.

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

        """
        vals = self.evaluate(distribution, weight='number', **kwargs) * self.dlogdp

        return np.array(vals)

    def surface_area(self, distribution, **kwargs):
        """Return the total surface area of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.

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

        """
        vals = self.evaluate(distribution, weight='surface', **kwargs) * self.dlogdp

        return np.array(vals)

    def volume(self, distribution, measured=True, **kwargs):
        """Return the total volume of particles an OPC 'sees' in each bin.

        Parameters
        ----------
        distribution : AerosolDistribution
            A valid AerosolDistribution instance that can be evaluated.

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


        What if we want to know the total volume across all bins?

        >>> d = opcsim.load_distribution("Urban")
        >>> opc = opcsim.OPC()
        >>> vt = opc.volume(d).sum()

        """
        vals = self.evaluate(distribution, weight='volume', **kwargs) * self.dlogdp

        return np.array(vals)

__all__ = [
    'OPC'
]
