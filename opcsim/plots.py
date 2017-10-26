import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import itertools

from .distributions import AerosolDistribution
from .models import OPC

lrg_number_fmt = mtick.ScalarFormatter()
lrg_number_fmt.set_powerlimits((-3, 4))

rc_log = {
    'xtick.major.size': 10.0,
    'xtick.minor.size': 6.0,
    'ytick.major.size': 10.0,
    'ytick.minor.size': 8.0,
    'xtick.color': '0.1',
    'ytick.color': '0.1',
    'axes.linewidth': 1.75
}

YLABEL = {
    'none': {
        'number': "$dN/dD_p \; [\mu m^{-1} cm^{-3}]$",
        'surface': "$dS/dD_p \; [\mu m cm^{-3}]$",
        'volume': "$dV/dD_p \; [\mu m^2 cm^{-3}]$",
        'mass': "$dM/dD_p$"
    },
    'log': {
        'number': "$dN/dlnD_p \; [# cm^{-3}]$",
        'surface': "$dS/dlnD_p \; [\mu m^2 cm^{-3}]$",
        'volume': "$dV/dlnD_p \; [\mu m^3 cm^{-3}]$",
        'mass': "$dM/dlnD_p \; [\mu g m^{-3}]$"
    },
    'log10': {
        'number': "$dN/dlogD_p \; [cm^{-3}]$",
        'surface': "$dS/dlogD_p \; [\mu m^2 cm^{-3}]$",
        'volume': "$dV/dlogD_p \; [\mu m^3 cm^{-3}]$",
        'mass': "$dM/dlogD_p \; [\mu g m^{-3}]$"
    }
}

YLABEL_CDF = {
    'number': 'Total Number [$cm^{-3}$]',
    'surface': 'Total Surface Area [$\mu m^2cm^{-3}$]',
    'volume': 'Total Volume [$\mu m^3 cm^{-3}$]',
    'mass': "Total Mass [$\mu g m^{-3}$]"
}

def histplot(data, bins, ax=None, plot_kws={}, fig_kws={}, **kwargs):
    """Plot the particle size distribution as a histogram/bar chart.

    Parameters
    ----------
    data : array of floats
        An array containing the y variable you are plotting (ex. dNdlogDp)
    bins : 3xn array
        A 3xn array containing the bin values for an OPC
    ax : matplotlib axis
        If an axis is provided, the histogram will be plotted on this axis.
        Otherwise, a new axis object will be created.
    plot_kws : dict
        Optional keyword arguments to include. They are sent as an argument to
        the matplotlib bar plot.
    fig_kws : dict
        Optional keyword arguments to include for the figure.

    Returns
    -------
    ax : matplotlib axis object

    Examples
    --------

    Plot a 10-bin OPC's response to the Urban Distribution

    .. plot::
        :context: close-figs

        >>> import opcsim, seaborn as sns
        >>> opc = opcsim.OPC(n_bins=10, dmin=0.3)
        >>> d = opcsim.load_distribution("Urban")
        >>> ax = opcsim.plots.histplot(opc.evaluate(d), opc.bins)
        >>> ax.set_ylabel("$dN/dlogD_p$")

    We can also plot the same OPC in volume (mass) space

    .. plot::
        :context: close-figs

        >>> ax = opcsim.plots.histplot(
        ...         opc.evaluate(d, weight='volume'), opc.bins)
        >>> ax.set_ylabel("$dV/dlogD_p$")

    How about overlaying two OPC's

    .. plot::
        :context: close-figs

        >>> opcb = opcsim.OPC(n_bins=5, dmin=0.3)
        >>> ax = opcsim.plots.histplot(opc.evaluate(d),
        ...             opc.bins, label="10 bin OPC")
        >>> ax = opcsim.plots.histplot(opcb.evaluate(d), opcb.bins,
        ...             label="5 bin OPC", ax=ax)
        >>> ax.set_ylabel("$dN/dlogD_p$")
        >>> ax.legend(loc='best')


    What if we want to fill in the boxes?

    .. plot::
        :context: close-figs

        >>> plot_kws = dict(fill=True)
        >>> ax = opcsim.plots.histplot(opc.evaluate(d),
        ...             opc.bins, plot_kws=plot_kws)
        >>> ax.set_ylabel("$dN/dlogD_p$")

    """
    # Set the default figure kws
    default_fig_kws = dict()
    fig_kws = dict(default_fig_kws, **fig_kws)

    # Make a new axis object if one wasnt' set
    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    # Get the next color available in the palette
    nc = next(ax._get_lines.prop_cycler)['color']

    # Set the plot_kws as a mapping of default and kwargs
    default_plot_kws = dict(
                        alpha=1,
                        edgecolor=nc,
                        color=nc,
                        linewidth=5,
                        fill=False,
                        label=kwargs.pop('label', None)
                        )

    # Set the plot_kws
    plot_kws = dict(default_plot_kws, **plot_kws)

    # Plot the bar plot
    ax.bar(left=bins[:, 0], height=data, width=bins[:, -1] - bins[:, 0],
            align='edge', **plot_kws)

    # Set the xaxis to be log10
    ax.semilogx()

    # Set the xlabel
    ax.set_xlabel("$D_p \; [\mu m]$")

    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.3g"))

    return ax

def pdfplot(distribution, ax=None, weight='number', base='log10', with_modes=False,
            fill=False, plot_kws={}, fig_kws={}, fill_kws={}, **kwargs):
    """Plot the PDF of a particle size distribution.

    Parameters
    ----------
    distribution : valid `AerosolDistribution`
        An aerosol distribution with the method `pdf` that can be evaluated at
        an array of particle diameters.
    weight : {'number' | 'surface' | 'volume' | 'mass'}
        Choose how to weight the pdf. Default is `number`.
    base : {'none' | 'log' | 'log10'}
        Base algorithm to use. Default is 'log10'.
    ax : matplotlib axis
        If an axis is provided, the histogram will be plotted on this axis.
        Otherwise, a new axis object will be created.
    with_modes : bool
        If true, all modes of a given distribution will be plotted along with
        their sum. If false, only the sum will be plotted.
    fill : bool
        If true, the area under the PDF will be filled. Cannot be used while
        plotting individual modes (with_modes=True).
    plot_kws : dict
        Optional keyword arguments to include. They are sent as an argument to
        the matplotlib plot call.
    fig_kws : dict
        Optional keyword arguments to include for the figure.
    fill_kws : dict
        Optional keyword arguments to include for the fill. Sent to
        matplotlib.axes.fill_between.

    Returns
    -------
    ax : matplotlib axis object

    Examples
    --------

    Plot the number-weighted Urban aerosol distribution

    .. plot::
        :context: close-figs

        >>> import opcsim, seaborn as sns
        >>> d = opcsim.load_distribution("Urban")
        >>> ax = opcsim.plots.pdfplot(d)
        >>> ax.set_title("Urban Aerosol Distribution", fontsize=16)
        >>> sns.despine()

    Let's take a look at the individual modes as well

    .. plot::
        :context: close-figs

        >>> ax = opcsim.plots.pdfplot(d, with_modes=True)
        >>> ax.set_title("Urban Aerosol Distribution", fontsize=16)
        >>> ax.legend(loc='best')
        >>> sns.despine()

    Let's plot the volume weighted version

    .. plot::
        :context: close-figs

        >>> ax = opcsim.plots.pdfplot(d, weight='volume', with_modes=True)
        >>> ax.set_title("Volume Weighted Urban Aerosol Distribution", fontsize=16)
        >>> ax.legend(loc='best')
        >>> sns.despine()

    Let's plot a few different distributions together

    .. plot::
        :context: close-figs

        >>> d2 = opcsim.load_distribution("Marine")
        >>> d3 = opcsim.load_distribution("Rural")
        >>> ax = opcsim.plots.pdfplot(d)
        >>> ax = opcsim.plots.pdfplot(d2, ax=ax)
        >>> ax = opcsim.plots.pdfplot(d3, ax=ax)
        >>> ax.set_title("Various Aerosol Distributions", fontsize=16)
        >>> ax.legend(loc='best')
        >>> sns.despine()

    """
    if not hasattr(distribution, 'pdf'):
        raise Exception("Invalid AerosolDistribution.")

    if weight not in ['number', 'surface', 'volume', 'mass']:
        raise ValueError("Invalid weight: ['number', 'surface', 'volume', 'mass']")

    # Set the default dp values to plot against
    dp = kwargs.get('dp', np.logspace(-3, 1, 1000))
    rho = kwargs.get('rho', 1.)

    # Set the default figure kws
    default_fig_kws = dict()
    fig_kws = dict(default_fig_kws, **fig_kws)

    # Make a new axis object if one wasnt' set
    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    # Plot the compete distribution
    nc = next(ax._get_lines.prop_cycler)['color']

    # Set the plot_kws as a mapping of default and kwargs
    default_plot_kws = dict(
                        alpha=1,
                        linewidth=4,
                        color=nc
                        )

    # Set the plot_kws
    plot_kws = dict(default_plot_kws, **plot_kws)

    # Set the default fill_kws
    default_fill_kws = dict(color=plot_kws['color'])
    fill_kws = dict(default_fill_kws, **fill_kws)

    # If label kwarg is present, use -> otherwise use default
    label = kwargs.pop('label', distribution.label)

    # Get the data to plot
    data = distribution.pdf(dp, base=base, weight=weight, rho=rho)

    # If fill is selected, fill the gap, otherwise just plot a line
    ax.plot(dp, data, c=nc, label=label, **plot_kws)

    if fill:
        ax.fill_between(dp, 0, data, label=label, **fill_kws)

    # Get data
    if with_modes is True:
        for m in distribution.modes:
            nc = next(ax._get_lines.prop_cycler)['color']

            data = distribution.pdf(dp, base=base, weight=weight, mode=m['label'])

            ax.plot(dp, data, label=m['label'], c=nc, ls='--', **plot_kws)

    if base is not ('none' or None):
        ax.semilogx()
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.4g"))

    ax.set_xlabel("$D_p \; [\mu m]$")
    ax.set_ylabel(YLABEL[base][weight])

    # Set the yaxis to show scientific notation when numbers > 1e4
    ax.yaxis.set_major_formatter(lrg_number_fmt)

    return ax

def cdfplot(distribution, ax=None, weight='number', plot_kws={},
            fig_kws={}, **kwargs):
    """Plot the CDF of a particle size distribution.

    Parameters
    ----------
    distribution : valid `AerosolDistribution`
        An aerosol distribution with the method `pdf` that can be evaluated at
        an array of particle diameters.
    weight : {'number' | 'surface' | 'volume' | 'mass'}
        Choose how to weight the pdf. Default is `number`.
    ax : matplotlib axis
        If an axis is provided, the histogram will be plotted on this axis.
        Otherwise, a new axis object will be created.
    plot_kws : dict
        Optional keyword arguments to include. They are sent as an argument to
        the matplotlib bar plot.
    fig_kws : dict
        Optional keyword arguments to include for the figure.

    Returns
    -------
    ax : matplotlib axis object

    Examples
    --------

    Plot the number-weighted Urban aerosol distribution

    .. plot::
        :context: close-figs

        >>> import opcsim, seaborn as sns
        >>> d = opcsim.load_distribution("Urban")
        >>> ax = opcsim.plots.cdfplot(d)
        >>> ax.set_title("Urban Aerosol Distribution", fontsize=16)
        >>> sns.despine()

    Let's plot the volume weighted version

    .. plot::
        :context: close-figs

        >>> ax = opcsim.plots.cdfplot(d, weight='volume')
        >>> ax.set_title("Volume Weighted Urban Aerosol Distribution", fontsize=16)
        >>> ax.legend(loc='best')
        >>> sns.despine()

    """
    if not hasattr(distribution, 'cdf'):
        raise Exception("Invalid AerosolDistribution.")

    if weight not in ['number', 'surface', 'volume', 'mass']:
        raise ValueError("Invalid weight: ['number', 'surface', 'volume', 'mass']")

    # Set the default dp values to plot against
    dp = kwargs.pop('dp', np.logspace(-3, 1, 1000))

    # Set the default figure kws
    default_fig_kws = dict()
    fig_kws = dict(default_fig_kws, **fig_kws)

    # Make a new axis object if one wasnt' set
    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    # Set the plot_kws as a mapping of default and kwargs
    default_plot_kws = dict(
                        alpha=1,
                        linewidth=5
                        )

    # Set the plot_kws
    plot_kws = dict(default_plot_kws, **plot_kws)

    # Plot the compete distribution
    nc = next(ax._get_lines.prop_cycler)['color']

    # If label kwarg is present, use -> otherwise use default
    label = kwargs.pop('label', distribution.label)

    data = distribution.cdf(dp, weight=weight, **kwargs)

    ax.plot(dp, data, c=nc, label=label, **plot_kws)

    ax.semilogx()
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.3g"))

    ax.set_xlabel("$D_p \; [\mu m]$")
    ax.set_ylabel(YLABEL_CDF[weight])

    # Set the yaxis to show scientific notation when numbers > 1e4
    ax.yaxis.set_major_formatter(lrg_number_fmt)

    return ax

__all__ = [
    'histplot',
    'pdfplot',
    'cdfplot'
]
