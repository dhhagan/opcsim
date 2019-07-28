"""Contains the scoring algorithms used in the model.
"""

import numpy as np
import pandas as pd
from .models import OPC
from .utils import k_kohler, ri_eff
from .mie import cscat


def compute_bin_assessment(opc, refr, kappa, rh_values=[0., 35., 95.]):
    """Assess the ability of an OPC to assign particles to their correct bin.

    Parameters
    ----------
    opc: opcsim.OPC
    refr: complex
        The complex refractive index of the material to assess
    kappa: float
        The kappa value to use for hygroscopic growth
    rh_values: list-like
        A list of relative humidities to assess the OPC at.

    Returns
    -------
    rv: pd.DataFrame
        A dataframe containing the results with self-explanatory columns.

    Examples
    --------

    """
    assert(isinstance(opc, OPC)), "opc must be an instance of the opcsim.OPC class"

    # init the dataframe to hold our results
    rv = pd.DataFrame()

    for rh in rh_values:
        for i, _bins in enumerate(opc.bins):
            # compute the wet diameter
            wet_diam_lo = k_kohler(diam_dry=_bins[0], kappa=kappa, rh=rh)
            wet_diam_hi = k_kohler(diam_dry=_bins[-1], kappa=kappa, rh=rh)

            # compute the pct_dry
            pct_dry = (_bins[0]**3) / (wet_diam_lo**3)

            # compute the effective RI
            ri = ri_eff(species=[refr, complex(1.333, 0)], weights=[pct_dry, 1-pct_dry])

            # compute the scattering cross-section
            cscat_lo_exp = cscat(
                dp=_bins[0], wl=opc.wl, refr=refr, theta1=opc.theta[0], theta2=opc.theta[1])
            cscat_hi_exp = cscat(
                dp=_bins[-1], wl=opc.wl, refr=refr, theta1=opc.theta[0], theta2=opc.theta[1])

            cscat_lo = cscat(
                dp=wet_diam_lo, wl=opc.wl, refr=ri, theta1=opc.theta[0], theta2=opc.theta[1])
            cscat_hi = cscat(
                dp=wet_diam_hi, wl=opc.wl, refr=ri, theta1=opc.theta[0], theta2=opc.theta[1])

            # assign bins
            bin_assign_lo = opc.calibration_function(values=[cscat_lo])
            bin_assign_hi = opc.calibration_function(values=[cscat_hi])

            # add results to the dataframe
            rv = rv.append({
                "bin_true": i,
                "bin_lo": bin_assign_lo[0] if len(bin_assign_lo) > 0 else -99,
                "bin_hi": bin_assign_hi[0] if len(bin_assign_hi) > 0 else -99,
                "refr_eff": ri,
                "rh": rh,
                "cscat_hi_ratio": cscat_hi / cscat_hi_exp,
                "cscat_lo_ratio": cscat_lo / cscat_lo_exp,
            }, ignore_index=True)
    
    # force datatypes to be correct
    rv["bin_true"] = rv["bin_true"].astype(int)
    rv["bin_lo"] = rv["bin_lo"].astype(int)
    rv["bin_hi"] = rv["bin_hi"].astype(int)
    rv["rh"] = rv["rh"].astype(float)
    rv["cscat_hi_ratio"] = rv["cscat_hi_ratio"].astype(float)
    rv["cscat_lo_ratio"] = rv["cscat_lo_ratio"].astype(float)

    return rv

# def nv_score(model, distribution, dmin=0.0, dmax=2.5, **kwargs):
#     """Calculate and return the number-to-volume ratio.

#     The total number of particles is calculated by calculating the total number
#     of particles in each individual bin, and then summing them. The total volume
#     in the distribution is calculated by integrating the Volume-weighted CDF
#     between 0 and `dmax` microns.

#     Parameters
#     ----------
#     model : OPC
#         A valid OPC model describing an OPC that can be evaluated.
#     distribution : AerosolDistribution
#         A valid AerosolDistribution instance that can be evaluated.
#     dmin : float
#         The minimum particle size to integrate the CDF under. Default is 0.0
#         microns.
#     dmax : float
#         The maximum particle size to integrate the CDF under. Default is 2.5
#         microns.

#     Returns
#     -------
#     N/V : float
#         Returns the number-to-volume ratio as a single float.

#     Examples
#     --------

#     Compute the number-to-volume ratio for a 2-bin OPC on the Urban distribution

#     >>> opc = opcsim.OPC(n_bins=2)
#     >>> urban = opcsim.load_distribution("Urban")
#     >>> n_v = opcsim.metrics.nv_score(opc, urban)

#     """
#     # evaluate the total number of particles in each bin (then sum)
#     total_number = model.number(distribution, **kwargs).sum()

#     # evaluate the total volume in the distribution < dmax
#     total_volume = distribution.cdf(weight='volume', dmax=dmax)

#     return total_number / total_volume

# def vv_score(model, distribution, dmin=0.0, dmax=2.5, **kwargs):
#     """Calculate and return the volume-to-volume ratio.

#     The total volume of particles per the OPC is calculated by calculating the
#     total number of particles in each individual bin, and then multiplying each
#     bin by a 'volume-factor'. The sum of individual bin volumes is then used.
#     The total volume in the distribution is calculated by integrating the
#     Volume-weighted CDF between 0 and `dmax` microns.

#     Parameters
#     ----------
#     model : OPC
#         A valid OPC model describing an OPC that can be evaluated.
#     distribution : AerosolDistribution
#         A valid AerosolDistribution instance that can be evaluated.
#     dmin : float
#         The minimum particle size to integrate the CDF under. Default is 0.0
#         microns.
#     dmax : float
#         The maximum particle size to integrate the CDF under. Default is 2.5
#         microns.

#     Returns
#     -------
#     V/V : float
#         Returns the volume-to-volume ratio as a single float.

#     Examples
#     --------

#     Compute the number-to-volume ratio for a 2-bin OPC on the Urban distribution

#     >>> opc = opcsim.OPC(n_bins=2)
#     >>> urban = opcsim.load_distribution("Urban")
#     >>> v_v = opcsim.metrics.vv_score(opc, urban)

#     """
#     # evaluate the total number of particles in each bin (then sum)
#     measured_volume = model.volume(distribution, **kwargs).sum()

#     # evaluate the total volume in the distribution < dmax
#     total_volume = distribution.cdf(weight='volume', dmin=dmin, dmax=dmax)

#     return measured_volume / total_volume
