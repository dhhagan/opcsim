"""
Visualize the Urban Aerosol Distribution
========================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks', font_scale=1.25)

# Load the example urban distribution
d = opcsim.load_distribution("Urban")

# Plot the number-weighted pdf with modes
ax = opcsim.plots.pdfplot(d, with_modes=True)

# Set the title and axes labels
ax.set_title("Urban Distribution", fontsize=18)

# Add a legend
ax.legend(loc='best')

# Set the ylim
ax.set_ylim(0, None)

# Remove the top and right spines
sns.despine()
