"""
10-Bin OPC Response with Distribution
=====================================
_thumb: .6, .5
"""
import seaborn as sns
import opcsim

sns.set(style='ticks', font_scale=1.5)

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Build a 10-bin OPC with a dmin of 300 nm
opc = opcsim.OPC(n_bins=10, dmin=0.3)

# Plot the histogram response
ax = opcsim.plots.histplot(opc.evaluate(d), bins=opc.bins)

# Overlay the distribution
ax = opcsim.plots.pdfplot(d, ax=ax, fill=True, fill_kws=dict(alpha=.2),
                            plot_kws=dict(linewidth=1))

# Remove the top and right spines
sns.despine()
