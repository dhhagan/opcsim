d2 = opcsim.load_distribution("Marine")
d3 = opcsim.load_distribution("Rural")
ax = opcsim.plots.pdfplot(d)
ax = opcsim.plots.pdfplot(d2, ax=ax)
ax = opcsim.plots.pdfplot(d3, ax=ax)
ax.set_title("Various Aerosol Distributions", fontsize=16)
ax.legend(loc='best')
sns.despine()
