import numpy as np
from scipy.stats import norm as normal
from sklearn.neighbors import NearestNeighbors

import matplotlib.pyplot as pl
from matplotlib import rcParams

from scvelo.tools.utils import groups_to_bool
from scvelo.tools.velocity_embedding import quiver_autoscale, velocity_embedding
from .docs import doc_params, doc_scatter
from .scatter import scatter
from .utils import (
    default_arrow,
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


# TODO: Add docstrings
def compute_velocity_on_grid(
    X_emb,
    V_emb,
    density=None,
    smooth=None,
    n_neighbors=None,
    min_mass=None,
    autoscale=True,
    adjust_for_stream=False,
    cutoff_perc=None,
):
    """TODO."""
    # remove invalid cells
    idx_valid = np.isfinite(X_emb.sum(1) + V_emb.sum(1))
    X_emb = X_emb[idx_valid]
    V_emb = V_emb[idx_valid]

    # prepare grid
    n_obs, n_dim = X_emb.shape
    density = 1 if density is None else density
    smooth = 0.5 if smooth is None else smooth

    grs = []
    for dim_i in range(n_dim):
        m, M = np.min(X_emb[:, dim_i]), np.max(X_emb[:, dim_i])
        m = m - 0.01 * np.abs(M - m)
        M = M + 0.01 * np.abs(M - m)
        gr = np.linspace(m, M, int(50 * density))
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    X_grid = np.vstack([i.flat for i in meshes_tuple]).T

    # estimate grid velocities
    if n_neighbors is None:
        n_neighbors = int(n_obs / 50)
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=-1)
    nn.fit(X_emb)
    dists, neighs = nn.kneighbors(X_grid)

    scale = np.mean([(g[1] - g[0]) for g in grs]) * smooth
    weight = normal.pdf(x=dists, scale=scale)
    p_mass = weight.sum(1)

    V_grid = (V_emb[neighs] * weight[:, :, None]).sum(1)
    V_grid /= np.maximum(1, p_mass)[:, None]
    if min_mass is None:
        min_mass = 1

    if adjust_for_stream:
        X_grid = np.stack([np.unique(X_grid[:, 0]), np.unique(X_grid[:, 1])])
        ns = int(np.sqrt(len(V_grid[:, 0])))
        V_grid = V_grid.T.reshape(2, ns, ns)

        mass = np.sqrt((V_grid**2).sum(0))
        min_mass = 10 ** (min_mass - 6)  # default min_mass = 1e-5
        min_mass = np.clip(min_mass, None, np.max(mass) * 0.9)
        cutoff = mass.reshape(V_grid[0].shape) < min_mass

        if cutoff_perc is None:
            cutoff_perc = 5
        length = np.sum(np.mean(np.abs(V_emb[neighs]), axis=1), axis=1).T
        length = length.reshape(ns, ns)
        cutoff |= length < np.percentile(length, cutoff_perc)

        V_grid[0][cutoff] = np.nan
    else:
        min_mass *= np.percentile(p_mass, 99) / 100
        X_grid, V_grid = X_grid[p_mass > min_mass], V_grid[p_mass > min_mass]

        if autoscale:
            V_grid /= 3 * quiver_autoscale(X_grid, V_grid)

    return X_grid, V_grid


@doc_params(scatter=doc_scatter)
def velocity_embedding_grid(
    adata,
    basis=None,
    vkey="velocity",
    density=None,
    smooth=None,
    min_mass=None,
    arrow_size=None,
    arrow_length=None,
    arrow_color=None,
    scale=None,
    autoscale=True,
    n_neighbors=None,
    recompute=None,
    X=None,
    V=None,
    X_grid=None,
    V_grid=None,
    principal_curve=False,
    color=None,
    use_raw=None,
    layer=None,
    color_map=None,
    colorbar=True,
    palette=None,
    size=None,
    alpha=0.2,
    perc=None,
    sort_order=True,
    groups=None,
    components=None,
    projection="2d",
    legend_loc="none",
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
    """Scatter plot of velocities on a grid.

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    density: `float` (default: 1)
        Amount of velocities to show - 0 none to 1 all
    arrow_size: `float` or triple `headlength, headwidth, headaxislength` (default: 1)
        Size of arrows.
    arrow_length: `float` (default: 1)
        Length of arrows.
    scale: `float` (default: 1)
        Length of velocities in the embedding.
    min_mass: `float` or `None` (default: `None`)
        Minimum threshold for mass to be shown.
        It can range between 0 (all velocities) and 100 (large velocities).
    smooth: `float` (default: 0.5)
        Multiplication factor for scale in Gaussian kernel around grid point.
    n_neighbors: `int` (default: None)
        Number of neighbors to consider around grid point.
    X: `np.ndarray` (default: None)
        embedding grid point coordinates
    V: `np.ndarray` (default: None)
        embedding grid velocity coordinates
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
            density=density,
            autoscale=autoscale,
            smooth=smooth,
            n_neighbors=n_neighbors,
            min_mass=min_mass,
        )

    scatter_kwargs = {
        "basis": basis,
        "perc": perc,
        "use_raw": use_raw,
        "sort_order": sort_order,
        "alpha": alpha,
        "components": components,
        "projection": projection,
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
                    velocity_embedding_grid(
                        adata,
                        density=density,
                        scale=scale,
                        size=size,
                        min_mass=min_mass,
                        smooth=smooth,
                        n_neighbors=n_neighbors,
                        principal_curve=principal_curve,
                        ax=pl.subplot(gs),
                        arrow_size=arrow_size,
                        arrow_length=arrow_length,
                        color=colors[i] if len(colors) > 1 else color,
                        layer=layers[i] if len(layers) > 1 else layer,
                        vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                        title=title[i] if isinstance(title, (list, tuple)) else title,
                        X_grid=None if len(vkeys) > 1 else X_grid,
                        V_grid=None if len(vkeys) > 1 else V_grid,
                        autoscale=False if len(vkeys) > 1 else autoscale,
                        **scatter_kwargs,
                        **kwargs,
                    )
                )
        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax

    else:
        ax, show = get_ax(ax, show, figsize, dpi)
        hl, hw, hal = default_arrow(arrow_size)
        if arrow_length is not None:
            scale = 1 / arrow_length
        if scale is None:
            scale = 1
        if arrow_color is None:
            arrow_color = "grey"
        quiver_kwargs = {"angles": "xy", "scale_units": "xy", "edgecolors": "k"}
        quiver_kwargs.update({"scale": scale, "width": 0.001, "headlength": hl / 2})
        quiver_kwargs.update({"headwidth": hw / 2, "headaxislength": hal / 2})
        quiver_kwargs.update({"color": arrow_color, "linewidth": 0.2, "zorder": 3})

        for arg in list(kwargs):
            if arg in quiver_kwargs:
                quiver_kwargs.update({arg: kwargs[arg]})
            else:
                scatter_kwargs.update({arg: kwargs[arg]})

        ax.quiver(
            X_grid[:, 0], X_grid[:, 1], V_grid[:, 0], V_grid[:, 1], **quiver_kwargs
        )

        if principal_curve:
            curve = adata.uns["principal_curve"]["projections"]
            pl.plot(curve[:, 0], curve[:, 1], c="w", lw=6, zorder=4)
            pl.plot(curve[:, 0], curve[:, 1], c="k", lw=3, zorder=5)

        size = 4 * default_size(adata) if size is None else size

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
