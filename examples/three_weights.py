"""
Plot the Urban Distribution in N, SA, and V
===========================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks')

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Set up a subplot with three subplots
fig, ax = plt.subplots(3, sharex=True)

# Plot the number-weighted pdf
opcsim.plots.pdfplot(d, ax=ax[0])

# Plot the SA-weighted pdf
opcsim.plots.pdfplot(d, weight='surface', ax=ax[1])

# Plot the volume-weighted pdf
opcsim.plots.pdfplot(d, weight='volume', ax=ax[2])

# Remove the top and right spines
sns.despine()
