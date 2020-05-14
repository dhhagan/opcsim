"""
Build Your Own Aerosol Distribution
===================================

_thumb: .6, .6
"""
import seaborn as sns
import opcsim

sns.set(style='ticks', font_scale=1.75)

# Load the example urban distribution
d = opcsim.AerosolDistribution()

# Add a mode with 1000 particles/cc, GM=100nm, and GSD=1.5
d.add_mode(n=1000, gm=0.1, gsd=1.5, label="Mode I", 
    kappa=0., rho=1.5, refr=complex(1.5, 0))

# Overlay the distribution
ax = opcsim.plots.pdfplot(d)

# Fix the y-axis to begin at 0
ax.set_ylim(0, None)

# Set the xlim to the "visible" region
ax.set_xlim(0.01, 1)

# Remove the top and right spines
sns.despine()