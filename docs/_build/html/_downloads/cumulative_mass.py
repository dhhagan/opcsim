"""
Compare Mass Loading Distributions for Several Samples
======================================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim

sns.set(style='ticks', font_scale=2)

# Load the example urban distribution
urban = opcsim.load_distribution("Urban")
rural = opcsim.load_distribution("Rural")

# Plot the mass-weighted cdf [urban]
ax = opcsim.plots.cdfplot(urban, weight='mass')

# Plot the mass-weighted cdf [rural]
ax = opcsim.plots.cdfplot(rural, weight='mass', ax=ax)

# Add a legend
ax.legend(loc='best')

# Remove the top and right spines
sns.despine()
