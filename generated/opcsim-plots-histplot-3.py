opcb = opcsim.OPC(n_bins=5, dmin=0.3)
ax = opcsim.plots.histplot(opc.evaluate(d),
            opc.bins, label="10 bin OPC")
ax = opcsim.plots.histplot(opcb.evaluate(d), opcb.bins,
            label="5 bin OPC", ax=ax)
ax.set_ylabel("$dN/dlogD_p$")
ax.legend(loc='best')
sns.despine()
