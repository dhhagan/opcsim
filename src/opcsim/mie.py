import numpy as np
import math
import scipy.special as bessel


def coef_pi_tau(theta, x):
    """Compute the angle-dependant functions (pi and tau) using upward recurrence.

    Per Bohren and Huffman (1983) eq. 4.47, the angle-dependant functions :math:`\pi_n` and 
    :math:`\\tau_n` can be computed using upward recurrance from the relations:

    .. math::

        \pi_n=\\frac{2n-1}{n-1}\mu\pi_{n-1}-\\frac{n}{n-1}\pi_{n-2} 
        
    .. math::
        
        \\tau_n=n\mu\pi_n - (n+1)\pi_{n-1}

    

    Parameters
    ----------
    theta: float
        The scattering angle in degrees.
    x: float
        The dimensionless size parameter, used to determine the number of elements to compute.

    Returns
    -------
    `\pi_n`, `\\tau_n`: np.ndarray of floats

    """
    # compute mu = cos(theta)
    mu = np.cos(math.radians(theta))

    # compute the max number of iterations
    nc = int(np.round(2 + x + 4*np.power(x, 1./3.)))

    # init arrays to hold the values
    pi, tau = np.zeros(nc), np.zeros(nc)

    # set the initial params
    pi[0] = 1
    pi[1] = 3*mu
    tau[0] = mu
    tau[1] = 3*np.cos(2*np.arccos(mu))

    # iterate and solve
    for n in range(2, nc):
        pi[n] = (mu*pi[n-1]*(2*n+1) - (pi[n-2]*(n+1)))/n
        tau[n] = (n+1)*mu*pi[n] - (n+2)*pi[n-1]

    return pi, tau


def coef_ab(refr, x):
    """Compute the external field coefficients using the logarithmic derivative.

    Bohren and Huffman (1983) equations 4.88 and 4.89 show how to compute :math:`a_n` and :math:`b_n`
    using the logarithmic derivative (Aden, 1951).

    .. math::

        a_n=\\frac{[D_n(mx)/m + n/x]\psi_n(x) - \psi_{n-1}(x)}{[D_n(mx)/m + n/x]\\xi_n(x) - \\xi_{n-1}(x)},

    .. math::

        b_n=\\frac{[mD_n(mx) + n/x]\psi_n(x) - \psi_{n-1}(x)}{[mD_n(mx) + n/x]\\xi_n(x) - \\xi_{n-1}(x)}

    where the logarithmic derivative is computed as:

    .. math::

        D_{n-1}=\\frac{n}{\\rho}-\\frac{1}{D_n + n/\\rho}

    
    Parameters
    ----------
    refr: complex
        The complex refractive index of the material.
    x: float
        The dimensionless size parameter.

    Returns
    -------
    `a_n`, `b_n`: np.ndarray of floats
        The external field coefficients.

    """
    # compute the number of values to calculate
    nc = int(np.round(2 + x + 4*np.power(x, 1/3)))

    # calculate z, the product of the RI and dimensionless size parameter
    z = refr*x

    nmx = int(np.round(max(nc, np.abs(z)) + 16))

    #
    n = np.arange(1, nc + 1)
    nu = n + 0.5

    # use scipy's bessel functions to compute
    sqx = np.sqrt(0.5 * np.pi * x)

    px = sqx * bessel.jv(nu, x)
    p1x = np.append(np.sin(x), px[0:nc-1])

    chx = -sqx*bessel.yv(nu, x)
    ch1x = np.append(np.cos(x), chx[0:nc-1])

    gsx = px - (0 + 1j)*chx
    gs1x = p1x - (0 + 1j)*ch1x

    # Bohren & Huffman eq. 4.89
    dn = np.zeros(nmx, dtype=np.complex128)

    for i in range(nmx-1, 1, -1):
        dn[i-1] = (i/z) - (1 / (dn[i] + i/z))

    # drop terms beyond nc
    d = dn[1:nc+1]

    da = d/refr + n/x
    db = refr*d + n/x

    an = (da*px - p1x) / (da*gsx - gs1x)
    bn = (db*px - p1x) / (db*gsx - gs1x)

    return an, bn


def s1s2(refr, x, theta):
    """Compute the complex scattering amplitudes S1 and S2 at angle theta.

    Bohren and Huffman (1983) list the equations for computing the complex scattering 
    amplitudes as Eq. 4.74:

    .. math::

        S_1=\sum_{n=1}^{n_c}\\frac{2n+1}{n(n+1)}(a_n\pi_n + b_n\\tau_n),

    .. math::

        S_2=\sum_{n=1}^{n_c}\\frac{2n+1}{n(n+1)}(a_n\\tau_n + b_n\pi_n)

    Parameters
    ----------
    refr: complex
        The complex refractive index.
    x: float
        The dimensionless size parameter.
    theta: float
        The scattering angle in degrees.

    Returns
    -------
    `S_1`, `S_2`: complex
        The complex scattering amplitudes.

    """
    # compute the number of coefficients to calculate
    nc = int(np.round(2 + x + 4*np.power(x, 1/3)))

    # compute the external field coefficients
    an, bn = coef_ab(refr=refr, x=x)

    # compute the pi and tau coefficients
    pi, tau = coef_pi_tau(theta=theta, x=x)

    # init an array for holding the S1 and S2 values
    n = np.arange(1, nc+1)
    
    # compute the coef for the series
    cn = (2*n + 1) / (n*(n+1))

    # compute the series and sum
    S1 = (cn * (an*pi + bn*tau)).sum()
    S2 = (cn * (an*tau + bn*pi)).sum()

    return S1, S2


def cscat(dp, wl, refr, theta1, theta2, nsteps=100, **kwargs):
    """Compute the scattering cross section between two angles according to Jaenicke and Hanusch (1993).

    Following the lead of Jaenicke and Hanusch (1993), we can compute the scattering cross section for 
    a given viewing angle [:math:`\Theta_1` - :math:`\Theta_2`] as:

    .. math::

        C_{sca}=\\frac{\lambda^2}{4\pi} \int_{\Theta_1}^{\Theta_2}[i_1(\Theta) + i_2(\Theta)]sin\Theta d\Theta

    where :math:`\lambda` is the incident wavelength of light and :math:`i_1` and :math:`i_2` are the intensity 
    distribution functions, calculated as:

    .. math::

        i_1(\Theta)=\mid S_1(\Theta) \mid^2,

    .. math::

        i_2(\Theta)=\mid S_2(\Theta) \mid^2


    The integral is calculated step-wise using the numpy.trapz function.
    

    Parameters
    ----------
    dp: float
        The particle diameter in microns.
    wl: float
        The wavelength of incident light in microns.
    refr: complex
        The complex refractive index of the material.
    theta1: float
        The angle from which to begin the integration.
    theta1: float
        The angle from which to end the integration.
    nsteps: int
        The number of steps in theta to use in performing the step-wise integration.
    
    Returns
    -------
    `C_{scat}`: float
        The scattering cross-section.

    """
    # build an array of angles
    thetas = np.linspace(theta1, theta2, nsteps)

    # compute the dimensionless parameter x
    x = dp*np.pi / wl

    # init an array to hold the inside of the integral
    rv = np.zeros(nsteps)
    
    # iterate over each step and compute the inside part of the integral
    for i in range(nsteps):
        s1, s2 = s1s2(refr=refr, x=x, theta=thetas[i])

        # compute i1 and i2
        i1 = (s1 * np.conjugate(s1)).real
        i2 = (s2 * np.conjugate(s2)).real

        rv[i] = (i1 + i2) * np.sin(math.radians(thetas[i]))
    
    # convert the array of angles from degrees to radians
    thetas = np.array([math.radians(x) for x in thetas])

    # compute cscat (convert the wavelength to cm to make units match common literature values)
    rv = ((wl*1e-4)**2 / (4*np.pi)) * np.trapz(rv, thetas)

    return rv
