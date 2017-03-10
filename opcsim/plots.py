import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import itertools

from .distributions import AerosolDistribution
from .models import OPC

rc_log = {
    'xtick.major.size': 8.0,
    'xtick.minor.size': 5.0,
    'ytick.major.size': 8.0,
    'ytick.minor.size': 5.0
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

    with sns.axes_style('ticks', rc_log):

        ax.bar(bins[:, 0], data, bins[:, -1] - bins[:, 0], **plot_kws)

        ax.semilogx()

        # Set the xlabel
        ax.set_xlabel("$D_p \; [\mu m]$")

        ax.xaxis.set_major_formatter(ScalarFormatter())

    return ax

"""
def distplot(distribution, x=np.linspace(0.01, 10., 1000), weight='number', base='log10', with_modes=False, **kwargs):
    #distplot

    #with_modes only works for an individual weight
    assert hasattr(distribution, 'pdf'), "Invalid AerosolDistribution"

    if type(weight) is not list:
        weight = [weight]

    for w in weight:
        if w not in ['number', 'surface', 'volume', 'mass']:
            raise Exception("Invalid weight: ['number', 'surface', 'volume', 'mass']")

    figsize     = kwargs.pop('figsize', (14, 7))
    fig         = kwargs.pop('fig', None)
    ax          = kwargs.pop('ax', None)
    plot_kws    = kwargs.pop('plot_kws', {'lw': 3.})

    with sns.axes_style('white', rc_log):
        if fig is None or ax is None:
            fig, ax = plt.subplots(len(weight), figsize = figsize, sharex = True)

        if len(weight) == 1:
            ax.plot(x, distribution.pdf(x, base, weight[0]), **plot_kws)

            if with_modes:
                for mode in distribution.modes:
                    ax.plot(x, distribution.pdf(x, base, weight[0], mode = mode['label']), label = mode['label'], linestyle = '--')

            ax.semilogx()
            ax.xaxis.set_major_formatter(ScalarFormatter())
        else:
            i = 0
            for w in weight:
                ax[i].plot(x, distribution.pdf(x, base, w))

                ax[i].set_ylabel(w)
                ax[i].semilogx()

                i = i + 1

            ax[0].xaxis.set_major_formatter(ScalarFormatter())

        sns.despine()

    return fig, ax
"""

__all__ = [
    'histplot'
]
