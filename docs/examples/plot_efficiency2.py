"""
Visualize OPC Counting Efficiency
=================================
_thumb: .55, .6
"""
import seaborn as sns
import numpy as np
import opcsim

sns.set(style='ticks', font_scale=1.5)

# Build an array of particle diameters
diams = np.logspace(-2, 1, 250)

# define an efficiency function
eff = lambda dp: 1 - np.exp(-5*dp)

# Load the example urban distribution
d = opcsim.AerosolDistribution()

# Add a mode with 1000 particles/cc, GM=100nm, and GSD=1.5
d.add_mode(n=1000, gm=0.1, gsd=1.5, label="Mode I")

# Overlay the distribution
ax = opcsim.plots.pdfplot(d, plot_kws={'linewidth': 6})

# Add the distribution that the OPC "sees" accounting for CE
ax.plot(diams, eff(diams)*d.pdf(diams), label="$\eta=1-exp(-5*D_p)$")

ax.set_ylim(0, None)
ax.set_xlim(.01, 1)

# Add the legend to the plot
ax.legend(loc='best')

# Remove the top and right spines
sns.despine(offset=5)
