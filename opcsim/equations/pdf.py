"""
This file contains all of the probability distribution functions
"""

import math
import numpy as np

def dn_ddp(Dp, N, GM, GSD):
    """S+P 8.34 [um-1 cm-3]
    n_N(D_p)
    """

    res = ( N / ( np.sqrt( 2 * np.pi ) * Dp * np.log( GSD ))) * \
            np.exp( -( np.log( Dp ) - np.log( GM )) ** 2 / ( 2 * np.log( GSD ) ** 2 ))

    return res

def dn_dlndp(Dp, N, GM, GSD):
    """S+P 8.21 [# cm-3]
    """
    return Dp * dn_ddp(Dp, N, GM, GSD)

def dn_dlogdp(Dp, N, GM, GSD):
    """S+P 8.18 [# cm-3]
    """
    return np.log(10) * Dp * dn_ddp(Dp, N, GM, GSD)

def ds_ddp(Dp, N, GM, GSD):
    """S+P 8.4 [um cm-3]
    """
    return np.pi * Dp ** 2 * dn_ddp(Dp, N, GM, GSD)

def ds_dlndp(Dp, N, GM, GSD):
    """S+P 8.10 [um2 cm-3]
    """
    return np.pi * Dp ** 2 * dn_dlndp(Dp, N, GM, GSD)

def ds_dlogdp(Dp, N, GM, GSD):
    """S+P 8.19
    """
    return np.log(10) * Dp * ds_ddp(Dp, N, GM, GSD)

def dv_ddp(Dp, N, GM, GSD):
    """S+P 8.6 [um2 cm-3]
    """
    return (np.pi / 6.) * Dp ** 3 * dn_ddp(Dp, N, GM, GSD)

def dv_dlndp(Dp, N, GM, GSD):
    """S+P 8.10 [um3 cm-3]
    """
    return (np.pi / 6.) * Dp ** 3 * dn_dlndp(Dp, N, GM, GSD)

def dv_dlogdp(Dp, N, GM, GSD):
    """S+P 8.20 [um3 cm-3]
    """
    return np.log(10) * Dp * dv_ddp(Dp, N, GM, GSD)
