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

# calibrate the OPC using PSl's
opc.calibrate("psl", method='spline')

# compute the values
vals = opc.histogram(d, weight="number", rh=0.0)

# Plot the histogram response
ax = opcsim.plots.histplot(vals, bins=opc.bins)

# Remove the top and right spines
sns.despine()
