"""Make changes to the default styling which seaborn doesn't support
"""
import matplotlib as mpl

def set(**kwargs):
    mpl.rcParams['figure.autolayout'] = kwargs.pop('figure.autolayout', False)
    mpl.rcParams['mathtext.default'] = kwargs.pop('mathtext.default', 'regular')

    return
