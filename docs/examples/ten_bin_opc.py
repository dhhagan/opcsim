"""
Urban Response of a 10-bin OPC
==============================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks')

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Build a 10-bin OPC
opc = opcsim.OPC(n_bins=10, dmin=0.3)

# Plot the histogram response
opcsim.plots.histplot(opc.evaluate(d), bins=opc.bins)

# Remove the top and right spines
sns.despine()
