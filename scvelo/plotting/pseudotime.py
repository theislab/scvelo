import numpy as np

import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator


def principal_curve(adata):
    X_curve = adata.uns["principal_curve"]["projections"]
    ixsort = adata.uns["principal_curve"]["ixsort"]
    pl.plot(X_curve[ixsort, 0], X_curve[ixsort, 1], c="k", lw=3, zorder=2000000)


def pseudotime(adata, gene_list, ckey="velocity", reverse=False):
    ixsort = adata.uns["principal_curve"]["ixsort"]
    arclength = adata.uns["principal_curve"]["arclength"]
    if reverse:
        arclength /= np.max(arclength)
    else:
        arclength = (np.max(arclength) - arclength) / np.max(arclength)
    cell_subset = adata.uns["principal_curve"]["cell_subset"]

    adata_subset = adata[cell_subset].copy()

    gs = pl.GridSpec(1, len(gene_list))
    for n, gene in enumerate(gene_list):
        i = adata_subset.var_names.get_loc(gene)
        ax = pl.subplot(gs[n])

        lb, ub = np.percentile(adata_subset.obsm[ckey][:, i], [0.5, 99.5])
        c = np.clip(adata_subset.obsm[ckey][ixsort, i], lb, ub)
        # pl.scatter(arclength[ixsort], adata2.obsm['Mu'][ixsort, i], label="unspliced")
        pl.scatter(
            arclength[ixsort],
            adata_subset.obsm["Ms"][ixsort, i]
            * adata_subset.uns["velocity_pars"]["gamma"][i],
            c=c,
            cmap="coolwarm",
            alpha=1,
            s=1,
            label="spliced",
        )

        c /= abs(c).max()
        c *= (
            adata_subset.obsm["Ms"][ixsort, i]
            * adata_subset.uns["velocity_pars"]["gamma"][i]
        ).max()

        z = np.ma.polyfit(arclength[ixsort], c, 8)
        fun = np.poly1d(z)
        pl.plot(arclength[ixsort], fun(arclength[ixsort]), label=ckey)

        # Hide the right and top spines
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))

        pl.ylabel(gene)
        pl.title(f"Colored by {ckey}")
