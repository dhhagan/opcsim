"""
Compare the Number and Mass-weighted Urban Distribution
=======================================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks', font_scale=1.25)

# Load the example urban distribution
urban = opcsim.load_distribution("Urban")

# Set up a subplot with three subplots
fig, ax = plt.subplots(2, sharex=True)

# Plot the number-weighted pdf
opcsim.plots.pdfplot(urban, ax=ax[0])

# Plot the mass-weighted pdf
opcsim.plots.pdfplot(urban, weight='mass', ax=ax[1])

# Remove axis label from top subplot
ax[0].set_xlabel("")

# Remove the top and right spines
sns.despine()
