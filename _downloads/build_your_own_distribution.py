"""
Build Your Own Aerosol Distribution
===================================
_thumb: .4, .4
"""
import seaborn as sns
import opcsim
sns.set(style='ticks')

# Load the example urban distribution
d = opcsim.AerosolDistribution()

# Add a mode with 1000 particles/cc, GM=100nm, and GSD=1.5
d.add_mode(n=1000, gm=0.1, gsd=1.5, label="Mode I")

# Overlay the distribution
ax = opcsim.plots.pdfplot(d)

# Remove the top and right spines
sns.despine()
