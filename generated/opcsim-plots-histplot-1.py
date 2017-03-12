import opcsim, seaborn as sns
opc = opcsim.OPC(n_bins=10, dmin=0.3)
d = opcsim.load_distribution("Urban")
ax = opcsim.plots.histplot(opc.evaluate(d), opc.bins)
ax.set_ylabel("$dN/dlogD_p$")
