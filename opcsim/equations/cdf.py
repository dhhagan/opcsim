import math
import numpy as np
from scipy.special import erf, erfc

def Nt(Dp, N, GM, GSD):
    """S+P 8.39

    Calculate the CDF

    [# cm-3]
    """
    res = ( N / 2. ) * (1 + erf( (np.log( Dp / GM )) / ( np.sqrt( 2 ) * np.log( GSD ))))

    return res

def St(Dp, N, GM, GSD):
    """Calculate the CDF for Surface Area. S+P 8.11

    Units: [um2 cm-3]
    """

    res = (np.pi / 2.) * N * ( GM ** 2 ) * np.exp( 2 * (np.log( GSD ) ** 2)) * \
                erfc((np.sqrt( 2 ) * np.log( GSD )) - ( np.log( Dp / GM ) / ( np.sqrt( 2 ) * np.log( GSD ) ) ))

    return res

def Vt(Dp, N, GM, GSD):
    """Calculate the CDF for Volume. S+P 8.12

    Units: [um3 cm-3]
    """

    res = (np.pi / 12.) * N * ( GM ** 3 ) * np.exp( 9. / 2. * (np.log( GSD ) ** 2)) * \
                erfc((1.5 * np.sqrt( 2 ) * np.log( GSD )) - ( np.log( Dp / GM ) / ( np.sqrt( 2 ) * np.log( GSD ) ) ))

    return res
