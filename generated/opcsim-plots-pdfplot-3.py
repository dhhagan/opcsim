ax = opcsim.plots.pdfplot(d, weight='volume', with_modes=True)
ax.set_title("Volume Weighted Urban Aerosol Distribution", fontsize=16)
ax.legend(loc='best')
sns.despine()
