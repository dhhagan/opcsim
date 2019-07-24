"""
Urban Response of a 10-bin OPC
==============================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks', font_scale=1.5)

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Build a 10-bin OPC
opc = opcsim.OPC(wl=0.658, n_bins=10, dmin=0.3)

# Plot the histogram response
# ax = opcsim.plots.histplot(opc.evaluate(d), bins=opc.bins)

# Remove the top and right spines
sns.despine()
