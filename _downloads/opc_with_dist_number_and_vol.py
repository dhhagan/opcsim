"""
10-Bin OPC Response in Number and Volume Space
==============================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks')

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Build a 10-bin OPC
opc = opcsim.OPC(n_bins=10, dmin=0.3)

# Set up the subplots
fig, (ax1, ax2) = plt.subplots(2)

# Plot the histogram response
ax1 = opcsim.plots.histplot(
                    opc.evaluate(d),
                    bins=opc.bins, ax=ax1)

# Overlay the distribution
ax1 = opcsim.plots.pdfplot(d, ax=ax1)

# Repeat the above step but weight by volume
ax2 = opcsim.plots.histplot(
                    opc.evaluate(d, weight='volume'),
                    bins=opc.bins, ax=ax2)

# Overlay the distribution
ax2 = opcsim.plots.pdfplot(d, weight='volume', ax=ax2)

# Remove the top and right spines
sns.despine()
