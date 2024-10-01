import numpy as np

import matplotlib.pyplot as pl

from scvelo.core import sum
from .utils import savefig_or_show


def proportions(
    adata,
    groupby="clusters",
    layers=None,
    highlight="unspliced",
    add_labels_pie=True,
    add_labels_bar=True,
    fontsize=8,
    figsize=(10, 2),
    dpi=100,
    use_raw=True,
    show=True,
    save=None,
):
    """Plot pie chart of spliced/unspliced proprtions.

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    groupby: `str` (default: 'clusters')
        Key of observations grouping to consider.
    layers: list of `str`(default: `['spliced#, 'unspliced', 'ambiguous']`)
        Specify the layers of count matrices for computing proportions.
    highlight: `str` (default: 'unspliced')
        Which proportions to highlight in pie chart.
    add_labels_pie: `bool` (default: True)
        Whether to add percentage labels in pie chart.
    add_labels_bar: `bool` (default: True)
        Whether to add percentage labels in bar chart.
    fontsize: `float` (default: 8)
        Label font size.
    figsize: tuple (default: `(10,2)`)
        Figure size.
    dpi: `int` (default: 80)
        Figure dpi.
    use_raw : `bool` (default: `True`)
        Use initial cell sizes before normalization and filtering.
    show: `bool` (default: True)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the default
        filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.

    Returns
    -------
    Plots the proportions of abundances as pie chart.
    """
    # get counts per cell for each layer
    if layers is None:
        layers = ["spliced", "unspliced", "ambiguous"]
    layers_keys = [key for key in layers if key in adata.layers.keys()]
    counts_layers = [sum(adata.layers[key], axis=1) for key in layers_keys]

    if use_raw:
        ikey, obs = "initial_size_", adata.obs
        counts_layers = [
            obs[ikey + layer_key] if ikey + layer_key in obs.keys() else c
            for layer_key, c in zip(layers_keys, counts_layers)
        ]
    counts_total = np.sum(counts_layers, 0)
    counts_total += counts_total == 0
    counts_layers = np.array([counts / counts_total for counts in counts_layers])

    gspec = pl.GridSpec(1, 2, pl.figure(None, figsize, dpi=dpi))
    colors = pl.get_cmap("tab20b")(np.linspace(0.10, 0.65, len(layers_keys)))

    # pie chart of total abundances
    ax = pl.subplot(gspec[0])
    if highlight is None:
        highlight = "none"
    explode = [
        0.1 if (layer_key == highlight or layer_key in highlight) else 0
        for layer_key in layers_keys
    ]

    autopct = "%1.0f%%" if add_labels_pie else None
    pie = ax.pie(
        np.mean(counts_layers, axis=1),
        colors=colors,
        explode=explode,
        autopct=autopct,
        shadow=True,
        startangle=45,
    )
    if autopct is not None:
        for pct, color in zip(pie[-1], colors):
            r, g, b, _ = color
            pct.set_color("white" if r * g * b < 0.5 else "darkgrey")
            pct.set_fontweight("bold")
            pct.set_fontsize(fontsize)
    ax.legend(
        layers_keys,
        ncol=len(layers_keys),
        bbox_to_anchor=(0, 1),
        loc="lower left",
        fontsize=fontsize,
    )

    # bar chart of abundances per category
    if groupby is not None and groupby in adata.obs.keys():
        counts_groups = {}
        for cluster in adata.obs[groupby].cat.categories:
            counts_groups[cluster] = np.mean(
                counts_layers[:, adata.obs[groupby] == cluster], axis=1
            )

        labels = list(counts_groups.keys())
        data = np.array(list(counts_groups.values()))
        data_cum = data.cumsum(axis=1)

        ax2 = pl.subplot(gspec[1])
        for i, (colname, color) in enumerate(zip(layers_keys, colors)):
            starts, widths = data_cum[:, i] - data[:, i], data[:, i]
            xpos = starts + widths / 2
            curr_xpos = xpos[0]
            for i, (x, w) in enumerate(zip(xpos, widths)):
                if not (x - w / 2 + 0.05 < curr_xpos < x + w / 2 - 0.05):
                    curr_xpos = x
                xpos[i] = curr_xpos

            ax2.barh(
                labels, widths, left=starts, height=0.9, label=colname, color=color
            )

            if add_labels_bar:
                r, g, b, _ = color
                text_color = "white" if r * g * b < 0.5 else "darkgrey"
                for y, (x, c) in enumerate(zip(xpos, widths)):
                    ax2.text(
                        x,
                        y,
                        f"{(c * 100):.0f}%",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=fontsize,
                        fontweight="bold",
                    )

        ax2.legend(
            ncol=len(layers_keys),
            bbox_to_anchor=(0, 1),
            loc="lower left",
            fontsize=fontsize,
        )
        ax2.invert_yaxis()
        ax2.set_xlim(0, np.nansum(data, axis=1).max())
        ax2.margins(0)

        ax2.set_xlabel("proportions", fontweight="bold", fontsize=fontsize * 1.2)
        ax2.set_ylabel(groupby, fontweight="bold", fontsize=fontsize * 1.2)
        ax2.tick_params(axis="both", which="major", labelsize=fontsize)
        ax = [ax, ax2]

    savefig_or_show("proportions", dpi=dpi, save=save, show=show)
    if show is False:
        return ax
