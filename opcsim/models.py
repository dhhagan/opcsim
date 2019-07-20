import numpy as np
import pandas as pd
import math
import scipy

from .distributions import AerosolDistribution
from .utils import make_bins, midpoints, squash_dips, power_law_fit
from .mie import cscat
import functools

RI_COMMON = {
    "psl": complex(1.59, 0),
    "ammonium_sulfate": complex(1.521, 0),
    "sodium_chloride": complex(1.5405, 0),
    "sodium_nitrate": complex(1.448, 0),
    "black_carbon": complex(1.95, 0.79),
    "sulfuric_acid": complex(1.427, 0),
    "soa": complex(1.4, 0.002),
}


class OPC(object):
    """Define an Optical Particle Counter (OPC).
    """
    def __init__(self, wl, bins=None, n_bins=None, dmin=0.5, 
            dmax=2.5, theta=(30., 90.), **kwargs):
        """

        Parameters
        ----------
        wl: float
            The laser wavelength in units of microns.
        bins: 3xn array
            An array of bin diameters containing the (left boundary, midpoint, right boundary).
        n_bins: int
            The number of desired bins. This should be used with dmin and dmax to generate 
            a 3xn array of bins.
        dmin: float
            The left-most bin boundary of the OPC.
        dmax: float
            The right-most bin boundary of the OPC.
        theta: tuple of floats
            The viewing range in units of degrees.
        
        Returns
        -------
        OPC2:
            An instance of the OPC2 class.
        
        Examples
        --------

        Initialize an OPC with 5 bins between 0.5 - 2.5 microns.

        >>> opc = opcsim.OPC2(wl=0.658, n_bins=5, dmin=0.5, dmax=2.5, theta=(30., 90.))

        Initialize an OPC with known bins.
        
        >>> opc = opcsim.OPC2(wl=0.658, bins=[0.38, 0.54, 0.78, 1.05, 1.5, 2.5], theta=(32., 88.))

        """
        # set some params
        self.n_bins = n_bins
        self.dmin = dmin
        self.dmax = dmax
        self.bins = bins
        self.wl = wl
        self.theta = theta
        self.label = kwargs.pop("label", None)
        self.calibration_function = None
        self._cscat_boundaries = None

        # generate the bins
        if self.bins is None:
            if self.n_bins is None:
                raise Exception("You must set either n_bins or bins.")

            self.bins = make_bins(self.dmin, self.dmax, self.n_bins)
        else:
            self.bins = np.asarray(self.bins)
            
            if self.bins.ndim == 1:
                self.bins = midpoints(np.array([self.bins[0:-1], self.bins[1:]]).T)
            
            # reset a few fields
            self.n_bins = self.bins.shape[0]
            self.dmin = self.bins[0, 0]
            self.dmax = self.bins[-1, -1]
        
        self.bin_boundaries = np.append(self.bins[:, 0], self.bins[-1, -1])

        return

    @property
    def dlogdp(self):
        return np.log10(self.bins[:, -1]) - np.log10(self.bins[:, 0])

    @property
    def ddp(self):
        return self.bins[:, -1] - self.bins[:, 0]

    @property
    def midpoints(self):
        return self.bins[:, 1]

    def calibrate(self, material, method='smooth', log_weight=True, mie_kws={}, fit_kws={}):
        """Calibrate the OPC assuming a specific material. 

        By calibration, we mean a method used to relate the peak height, 
        which is related to the scattering cross-section (Cscat) to the 
        particle size. At its simplest, we use this function to determine 
        which bin a particle belongs to based on its Cscat value. Once 
        calibrated, a calibration function is saved as a digitizer which 
        will take as an input a Cscat value and return the bin number of 
        the OPC it belongs to.

        Parameters
        ----------
        material: string or complex number
            Either a string containing the material available in the
            lookup table, or the complex refractive index. The option for 
            lookup values are: ['ammonium_sulfate', 'bc', 'psl']. Since 
            the refractive index is wavelength dependant, it is 
            recommended you define the refractive index at your OPC's 
            wavelength if you are looking for the best results.
        method: string or callable
            The method to use for creating a calibration curve. Options 
            include (1) 'smooth' which removes any non-monotonicly increasing 
            points from the Cscat to Dp curve; or (2) 'linear' fits a linear model 
            (in log-log space) between Cscat and Dp.
        log_weight: bool
            When generating a calibration curve using the 'linear_fit' method,
            setting log_weight=True will generate an error matrix that weights 
            the fit in log-units.
        mie_kws: dict
            Optional dictionary containing keyword arguments that is sent 
            directly to opcsim.mie.cscat when computing the scattering 
            cross section values used in the optimization.
        fit_kws: dict
            Optional dictionary containing keyword arguments that is 
            sent directly to scipy.optimize.curve_fit when generating 
            a fit using the 'fit_linear' approach. Please see the 
            scipy.optimize.curve_fit docs for more details.

        Examples
        --------

        """
        # determine the complex refractive index
        if type(material) == str:
            try:
                refr = RI_COMMON[material.lower()]
            except:
                raise ValueError("This material does not exist in our (limited) database.")
        else:
            if type(material) != complex:
                refr = complex(material, 0)
            else:
                refr = material
        
        # calculate Cscat at all bin boundaries
        yvals = np.array([cscat(dp=dp, wl=self.wl, refr=refr, theta1=self.theta[0],
                                theta2=self.theta[1], **mie_kws) for dp in self.bin_boundaries])
        
        # generate the fitted Cscat values based on the method chosen
        if method == "smooth":
            yvals = squash_dips(yvals)
        elif method == "fit_linear":
            # define the function we fit dp to Cscat
            def f(dp, a, b):
                return a * np.power(dp, b)
            
            # if log-weight, set sigma to weight in log-units
            if log_weight:
                fit_kws["sigma"] = np.power(10, np.log10(yvals) + 1)

            # fit the data with optional kwargs
            popt, _ = scipy.optimize.curve_fit(f, self.bin_boundaries, yvals, **fit_kws)

            # set the yvals to the fitted values
            yvals = f(self.bin_boundaries, *popt)
        else:
            raise ValueError("The `method` chosen is not currently supported.")

        # generate the digitizer and set as the calibration function
        # the calibration function must take just one argument (cscat(s)) and return 
        # a bin number(s)
        self.calibration_function = functools.partial(self._digitize_opc_bins, cscat_boundaries=yvals)

        # save the fitted boundaries for potential future use
        self._cscat_boundaries = yvals

    def evaluate(self, distribution, **kwargs):
        """Return the histogram for a given AerosolDistribution.

        We evaluate an OPC for a given distribution by calculating the Cscat value
        for every particle in the distribution and assigning it to a bin of the OPC.
        Since we are beginning with a PDF and not a distribution, we must first 
        discretize our PDF - we can choose the number of bins to convert of PDF 
        into which allows us to limit the computation needed by only performing 
        calculations for a small subset of the actual particles (i.e. we can do one
        calculation for a tiny bin and then replicate it without needing to re-do the 
        Mie calculations).


        Parameters
        ----------
        distribution: AerosolDistribution
            A valid instance of the AerosolDistribution class.

        Returns
        -------
        dN: array
            The number of particles in each OPC bin.

        Examples
        --------

        """
        if not self.calibration_function:
            raise Exception("The OPC must be calibrated before computing a histogram.")

        # create the bins for the distribution PDF
        bounds = kwargs.pop("bounds", (self.dmin / 2, 10.))

        # calculate the total number of particles in the distribution
        ntot = distribution.cdf(dmin=0, dmax=100., weight='number')

        bounds = np.logspace(start=np.log10(bounds[0]), stop=np.log10(
            bounds[1]), num=kwargs.pop("n_bins", int(min(ntot/10, 250))))

        # for each mode...
        rv = np.array([])
        for m in distribution.modes:
            # divide our PDF into bins to make the computations a bit easier
            n = np.array([round(distribution.cdf(dmin=a, dmax=b, mode=m["label"]), 0)
                                for a, b in zip(bounds[0:-1], bounds[1:])])
        
            # iterate over each bin and calculate the Cscat value and build an array
            for dp, dn in zip(self.midpoints, n):
                v = cscat(
                    dp, wl=self.wl, refr=m["refr"], theta1=self.theta[0], theta2=self.theta[1])
                
                rv = np.append(rv, np.repeat(np.array([v]), int(dn)))
        
        # convert the array of Cscat values into counts per bin
        rv = np.array([np.count_nonzero(self.calibration_function(values=rv) == x)
                       for x in np.arange(self.n_bins)])
        
        # force to be floats
        rv = rv.astype(float)

        return rv
    
    def histogram(self, distribution, weight="number", base="log10", **kwargs):
        """Return a histogram containing the [weight] of particles in each OPC bin.

        This represents what the OPC 'sees'. All calculations are made assuming the 
        center of the bin is reprentative of the entire bin.

        Parameters
        ----------
        distribution: AerosolDistribution
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`.
        base : {'none' | 'log10'}
            Base algorithm to use. Default is 'log10'.

        Returns
        -------
        left_bin, d[weight]/d[base]Dp, bin_width : arrays
            Returns three arrays containing the left bin boundary, the evaluated
            PDF, and the width of the bin (dDp). This data can be
            directly plotted as a histogram using matplotlib bar plots. By 
            default, dN/dlogDp is returned.

        See Also
        --------

        Examples
        --------

        """
        # get the density if needed [units of g/cc]
        rho = kwargs.pop("rho", 1.65)

        # get the binned values in units of dN/bin
        rv = self.evaluate(distribution, **kwargs)

        if weight not in ["number", "surface", "volume", "mass"]:
            raise ValueError("Invalid `weight` parameter")

        # convert the values into the proper weight
        if weight == "surface":
            rv *= (np.pi * self.midpoints**2)
        elif weight == "volume":
            rv *= (self.midpoints**3*np.pi/6.)  # dp3 * pi/6*N
        elif weight == "mass":
            rv *= rho * (self.midpoints**3*np.pi/6.)  # rho*dp3 * pi/6*N
        else:
            pass

        # divide the values by the proper base
        if base == "log10":
            rv /= self.dlogdp
        else:
            rv /= self.ddp

        return self.bins[:, 0], rv, self.ddp
    
    def integrate(self, distribution, dmin=0., dmax=1., weight="number", **kwargs):
        """Integrate the distribution according to the OPC for any [weight].

        By default, this method returns the total number of particles between 0-1 microns.
        It can be used to calculate PM values at any breakpoint. When calculating a PM value
        for a breakpoint that falls within a bin, a simple percentage of the bin is used. For
        example, if you calculate PM1 and the bin that contains 1-micron actually covers 0.5 
        to 1.5 microns, then 1/2 the mass from that bin will count towards the PM1 value.

        Parameters
        ----------
        distribution: AerosolDistribution
            The aerosol distribution to be evaluated.
        dmin: float
            The minimum particle diameter [microns] to integrate from.
        dmax: float
            The maximum particle diameter [microns] to integrate to.
        weight : {'number' | 'surface' | 'volume'}
            Choose how to weight the pdf. Default is `number`.

        Returns
        -------
        rv: float
            The total [weight] between dmin and dmax. By default, 
            the total number of particles (i.e. weight='number') 
            are returned.

        See Also
        --------

        Examples
        --------

        """
        rho = kwargs.pop("rho", 1.65)

        # calculate dN
        rv = self.evaluate(distribution, **kwargs)

        if weight not in ["number", "surface", "volume", "mass"]:
            raise ValueError("Invalid argument for `weight`")
        
        if weight == "surface":
            # convert dN to dS
            rv *= (np.pi * self.midpoints**2)
        elif weight == "volume":
            # convert dN to dV
            rv *= (self.midpoints**3*np.pi/6.)
        elif weight == "mass":
            # convert dN to dM
            rv *= rho*(self.midpoints**3*np.pi/6.)
        else:
            pass

        # cleave the array
        # for bins that are outside diameters of interest, set to 0
        # for bins that are partially within the range of diams, set to 
        # a fraction based on % within range
        factors = np.zeros(self.bins.shape[0])
        for i, b in enumerate(self.bins):
            if (b[0] >= dmin) and (b[-1] <= dmax):
                factors[i] = 1.
            elif (dmin >= b[0]) and (dmax <= b[-1]):
                factors[i] = (dmax - dmin) / (b[-1] - b[0])
            elif (dmin >= b[0]) and (dmax >= b[-1]):
                factors[i] = (b[-1] - dmin) / (b[-1] - b[0])
            elif (dmax >= b[0]) and (dmax <= b[-1]):
                factors[i] = (dmax - b[-1]) / (b[-1] / b[0])
            else:
                factors[i] = 0.

        # return the sum
        return (rv*factors).sum()

    def _digitize_opc_bins(self, cscat_boundaries, values):
        """Return the bin (or bins) corresponding to the :math:`C_{scat}` value(s).

        Parameters
        ----------
        cscat_boundaries: array-like
            An array of Cscat values at every bin boundary. The length of 
            the array should be one greater than the number of bins.
        values: array-like
            Either a single-value or array of computed Cscat values 
            that need to be assigned to OPC bins.

        Returns
        -------
        bins: np.ndarray
            An array of assigned bins associated with every value. If 
            a value was too small (i.e. value < smallest in cscat_boundaries)
            or too large, it will be dropped. Thus, the size of bins is less than 
            or equal to the size of values.

        Examples
        --------
        

        """
        _binb = np.asarray(cscat_boundaries)
        _vals = np.asarray(values)

        # create a digitized version of all values. 0 implies it is too small and
        # values that are too large will be assigned the value for the length of the array
        digitized = np.digitize(_vals, bins=_binb)

        # keep only valid bins
        digitized = digitized[np.where(
            (digitized != 0) & (digitized != len(_binb)))]

        # subtract 1 from every bin so they are 0-indexed
        digitized -= 1

        return digitized

    def __repr__(self): # pragma: no cover
        if self.label:
            return "{} (ƛ = {:.0f} nm)".format(self.label, self.wl*1e3)

        return self


__all__ = [
    'OPC'
]
