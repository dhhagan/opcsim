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
