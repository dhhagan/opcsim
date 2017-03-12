plot_kws = dict(fill=True)
ax = opcsim.plots.histplot(opc.evaluate(d),
            opc.bins, plot_kws=plot_kws)
ax.set_ylabel("$dN/dlogD_p$")
