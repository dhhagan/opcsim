ax = opcsim.plots.histplot(
        opc.evaluate(d, weight='volume'), opc.bins)
ax.set_ylabel("$dV/dlogD_p$")
