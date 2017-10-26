import opcsim, seaborn as sns
d = opcsim.load_distribution("Urban")
ax = opcsim.plots.pdfplot(d)
ax.set_title("Urban Aerosol Distribution", fontsize=16)
sns.despine()
