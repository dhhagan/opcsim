ax = opcsim.plots.cdfplot(d, weight='volume')
ax.set_title("Volume Weighted Urban Aerosol Distribution", fontsize=16)
ax.legend(loc='best')
sns.despine()
