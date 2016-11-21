import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from .distributions import AerosolDistribution
from .models import OPC

rc_log = {
    'xtick.major.size':     8.0,
    'xtick.minor.size':    5.0,
    'ytick.major.size':     8.0
}

def distplot(distribution, x = np.linspace(0.01, 10., 1000), weight = 'number', base = 'log10', with_modes = False, **kwargs):
    """
        distplot

        with_modes only works for an individual weight
    """
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


__all__ = [
    'distplot'
]
