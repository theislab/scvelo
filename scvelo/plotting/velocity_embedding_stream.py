import numpy as np

import matplotlib.pyplot as pl
from matplotlib import rcParams

from scvelo.tools.utils import groups_to_bool
from scvelo.tools.velocity_embedding import velocity_embedding
from .docs import doc_params, doc_scatter
from .scatter import scatter
from .utils import (
    default_basis,
    default_color,
    default_size,
    get_ax,
    get_basis,
    get_components,
    get_figure_params,
    make_unique_list,
    savefig_or_show,
    velocity_embedding_changed,
)
from .velocity_embedding_grid import compute_velocity_on_grid


@doc_params(scatter=doc_scatter)
def velocity_embedding_stream(
    adata,
    basis=None,
    vkey="velocity",
    density=2,
    smooth=None,
    min_mass=None,
    cutoff_perc=None,
    arrow_color=None,
    arrow_size=1,
    arrow_style="-|>",
    max_length=4,
    integration_direction="both",
    linewidth=None,
    n_neighbors=None,
    recompute=None,
    color=None,
    use_raw=None,
    layer=None,
    color_map=None,
    colorbar=True,
    palette=None,
    size=None,
    alpha=0.3,
    perc=None,
    X=None,
    V=None,
    X_grid=None,
    V_grid=None,
    sort_order=True,
    groups=None,
    components=None,
    legend_loc="on data",
    legend_fontsize=None,
    legend_fontweight=None,
    xlabel=None,
    ylabel=None,
    title=None,
    fontsize=None,
    figsize=None,
    dpi=None,
    frameon=None,
    show=None,
    save=None,
    ax=None,
    ncols=None,
    **kwargs,
):
    """\
    Stream plot of velocities on the embedding.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    density: `float` (default: 2)
        Controls the closeness of streamlines. When density = 2 (default), the domain
        is divided into a 60x60 grid, whereas density linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use a tuple (density_x, density_y).
    smooth: `float` (default: 0.5)
        Multiplication factor for scale in Gaussian kernel around grid point.
    min_mass: `float` (default: 1)
        Minimum threshold for mass to be shown.
        It can range between 0 (all velocities) and 5 (large velocities only).
    cutoff_perc: `float` (default: `None`)
        If set, mask small velocities below a percentile threshold (between 0 and 100).
    linewidth: `float` (default: 1)
        Line width for streamplot.
    arrow_color: `str` or 2D array (default: 'k')
        The streamline color. If given an array, it must have the same shape as u and v.
    arrow_size: `float` (default: 1)
        Scaling factor for the arrow size.
    arrow_style: `str` (default: '-|>')
        Arrow style specification, '-|>' or '->'.
    max_length: `float` (default: 4)
        Maximum length of streamline in axes coordinates.
    integration_direction: `str` (default: 'both')
        Integrate the streamline in 'forward', 'backward' or 'both' directions.
    n_neighbors: `int` (default: None)
        Number of neighbors to consider around grid point.
    X: `np.ndarray` (default: None)
        Embedding coordinates. Using `adata.obsm['X_umap']` per default.
    V: `np.ndarray` (default: None)
        Embedding velocity coordinates. Using `adata.obsm['velocity_umap']` per default.

    {scatter}

    Returns
    -------
    `matplotlib.Axis` if `show==False`
    """

    basis = default_basis(adata, **kwargs) if basis is None else get_basis(adata, basis)
    if vkey == "all":
        lkeys = list(adata.layers.keys())
        vkey = [key for key in lkeys if "velocity" in key and "_u" not in key]
    color, color_map = kwargs.pop("c", color), kwargs.pop("cmap", color_map)
    colors = make_unique_list(color, allow_array=True)
    layers, vkeys = make_unique_list(layer), make_unique_list(vkey)

    if V is None:
        for key in vkeys:
            if recompute or velocity_embedding_changed(adata, basis=basis, vkey=key):
                velocity_embedding(adata, basis=basis, vkey=key)

    color, layer, vkey = colors[0], layers[0], vkeys[0]
    color = default_color(adata) if color is None else color

    if X_grid is None or V_grid is None:
        _adata = (
            adata[groups_to_bool(adata, groups, groupby=color)]
            if groups is not None and color in adata.obs.keys()
            else adata
        )
        comps, obsm = get_components(components, basis), _adata.obsm
        X_emb = np.array(obsm[f"X_{basis}"][:, comps]) if X is None else X[:, :2]
        V_emb = np.array(obsm[f"{vkey}_{basis}"][:, comps]) if V is None else V[:, :2]
        X_grid, V_grid = compute_velocity_on_grid(
            X_emb=X_emb,
            V_emb=V_emb,
            density=1,
            smooth=smooth,
            min_mass=min_mass,
            n_neighbors=n_neighbors,
            autoscale=False,
            adjust_for_stream=True,
            cutoff_perc=cutoff_perc,
        )
        lengths = np.sqrt((V_grid ** 2).sum(0))
        linewidth = 1 if linewidth is None else linewidth
        linewidth *= 2 * lengths / lengths[~np.isnan(lengths)].max()

    scatter_kwargs = {
        "basis": basis,
        "perc": perc,
        "use_raw": use_raw,
        "sort_order": sort_order,
        "alpha": alpha,
        "components": components,
        "legend_loc": legend_loc,
        "groups": groups,
        "legend_fontsize": legend_fontsize,
        "legend_fontweight": legend_fontweight,
        "palette": palette,
        "color_map": color_map,
        "frameon": frameon,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "colorbar": colorbar,
        "dpi": dpi,
        "fontsize": fontsize,
        "show": False,
        "save": False,
    }

    stream_kwargs = {
        "linewidth": linewidth,
        "density": density or 2,
        "zorder": 3,
        "arrow_color": arrow_color or "k",
        "arrowsize": arrow_size or 1,
        "arrowstyle": arrow_style or "-|>",
        "maxlength": max_length or 4,
        "integration_direction": integration_direction or "both",
    }

    multikey = (
        colors
        if len(colors) > 1
        else layers
        if len(layers) > 1
        else vkeys
        if len(vkeys) > 1
        else None
    )
    if multikey is not None:
        if title is None:
            title = list(multikey)
        elif isinstance(title, (list, tuple)):
            title *= int(np.ceil(len(multikey) / len(title)))
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams["figure.figsize"] if figsize is None else figsize
        figsize, dpi = get_figure_params(figsize, dpi, ncols)
        gs_figsize = (figsize[0] * ncols, figsize[1] * nrows)
        ax = []
        for i, gs in enumerate(
            pl.GridSpec(nrows, ncols, pl.figure(None, gs_figsize, dpi=dpi))
        ):
            if i < len(multikey):
                ax.append(
                    velocity_embedding_stream(
                        adata,
                        size=size,
                        smooth=smooth,
                        n_neighbors=n_neighbors,
                        ax=pl.subplot(gs),
                        color=colors[i] if len(colors) > 1 else color,
                        layer=layers[i] if len(layers) > 1 else layer,
                        vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                        title=title[i] if isinstance(title, (list, tuple)) else title,
                        X_grid=None if len(vkeys) > 1 else X_grid,
                        V_grid=None if len(vkeys) > 1 else V_grid,
                        **scatter_kwargs,
                        **stream_kwargs,
                        **kwargs,
                    )
                )
        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax

    else:
        ax, show = get_ax(ax, show, figsize, dpi)

        for arg in list(kwargs):
            if arg in stream_kwargs:
                stream_kwargs.update({arg: kwargs[arg]})
            else:
                scatter_kwargs.update({arg: kwargs[arg]})

        stream_kwargs["color"] = stream_kwargs.pop("arrow_color", "k")
        ax.streamplot(X_grid[0], X_grid[1], V_grid[0], V_grid[1], **stream_kwargs)

        size = 8 * default_size(adata) if size is None else size
        ax = scatter(
            adata,
            layer=layer,
            color=color,
            size=size,
            title=title,
            ax=ax,
            zorder=0,
            **scatter_kwargs,
        )

        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax
