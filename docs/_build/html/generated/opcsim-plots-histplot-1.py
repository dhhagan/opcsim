import opcsim, seaborn as sns
opc = opcsim.OPC(wl=0.658, n_bins=10, dmin=0.3)
opc.calibrate("psl")
d = opcsim.load_distribution("Urban")
ax = opcsim.plots.histplot(opc.evaluate(d), opc.bins)
ax.set_ylabel("$dN/dlogD_p$")
sns.despine()
