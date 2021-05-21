import numpy as np

import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import is_color_like

from scvelo.tools.utils import groups_to_bool
from scvelo.tools.velocity_embedding import (
    velocity_embedding as compute_velocity_embedding,
)
from .docs import doc_params, doc_scatter
from .scatter import scatter
from .utils import (
    default_arrow,
    default_basis,
    default_color,
    default_color_map,
    default_size,
    get_ax,
    get_components,
    get_figure_params,
    interpret_colorkey,
    make_unique_list,
    make_unique_valid_list,
    savefig_or_show,
    velocity_embedding_changed,
)


@doc_params(scatter=doc_scatter)
def velocity_embedding(
    adata,
    basis=None,
    vkey="velocity",
    density=None,
    arrow_size=None,
    arrow_length=None,
    scale=None,
    X=None,
    V=None,
    recompute=None,
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
    """\
    Scatter plot of velocities on the embedding.

    Arguments
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
    {scatter}

    Returns
    -------
    `matplotlib.Axis` if `show==False`
    """

    if vkey == "all":
        lkeys = list(adata.layers.keys())
        vkey = [key for key in lkeys if "velocity" in key and "_u" not in key]
    color, color_map = kwargs.pop("c", color), kwargs.pop("cmap", color_map)
    layers, vkeys = make_unique_list(layer), make_unique_list(vkey)
    colors = make_unique_list(color, allow_array=True)
    bases = make_unique_valid_list(adata, basis)
    bases = [default_basis(adata, **kwargs) if b is None else b for b in bases]

    if V is None:
        for key in vkeys:
            for bas in bases:
                if recompute or velocity_embedding_changed(adata, basis=bas, vkey=key):
                    compute_velocity_embedding(adata, basis=bas, vkey=key)

    scatter_kwargs = {
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
        else bases
        if len(bases) > 1
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
                    velocity_embedding(
                        adata,
                        density=density,
                        scale=scale,
                        size=size,
                        ax=pl.subplot(gs),
                        arrow_size=arrow_size,
                        arrow_length=arrow_length,
                        basis=bases[i] if len(bases) > 1 else basis,
                        color=colors[i] if len(colors) > 1 else color,
                        layer=layers[i] if len(layers) > 1 else layer,
                        vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                        title=title[i] if isinstance(title, (list, tuple)) else title,
                        **scatter_kwargs,
                        **kwargs,
                    )
                )
        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax

    else:
        ax, show = get_ax(ax, show, figsize, dpi, projection)

        color, layer, vkey, basis = colors[0], layers[0], vkeys[0], bases[0]
        color = default_color(adata) if color is None else color
        color_map = default_color_map(adata, color) if color_map is None else color_map
        size = default_size(adata) / 2 if size is None else size
        if use_raw is None and "Ms" not in adata.layers.keys():
            use_raw = True
        _adata = (
            adata[groups_to_bool(adata, groups, groupby=color)]
            if groups is not None and color in adata.obs.keys()
            else adata
        )

        quiver_kwargs = {
            "scale": scale,
            "cmap": color_map,
            "angles": "xy",
            "scale_units": "xy",
            "edgecolors": "k",
            "linewidth": 0.1,
            "width": None,
        }
        if basis in adata.var_names:
            if use_raw:
                x = adata[:, basis].layers["spliced"]
                y = adata[:, basis].layers["unspliced"]
            else:
                x = adata[:, basis].layers["Ms"]
                y = adata[:, basis].layers["Mu"]
            dx = adata[:, basis].layers[vkey]
            dy = np.zeros(adata.n_obs)
            if f"{vkey}_u" in adata.layers.keys():
                dy = adata[:, basis].layers[f"{vkey}_u"]
            X = np.stack([np.ravel(x), np.ravel(y)]).T
            V = np.stack([np.ravel(dx), np.ravel(dy)]).T
        else:
            x = None if X is None else X[:, 0]
            y = None if X is None else X[:, 1]
            comps = get_components(components, basis, projection)
            X = _adata.obsm[f"X_{basis}"][:, comps] if X is None else X[:, :2]
            V = _adata.obsm[f"{vkey}_{basis}"][:, comps] if V is None else V[:, :2]

            hl, hw, hal = default_arrow(arrow_size)
            if arrow_length is not None:
                scale = 1 / arrow_length
            if scale is None:
                scale = 1
            quiver_kwargs.update({"scale": scale, "width": 0.0005, "headlength": hl})
            quiver_kwargs.update({"headwidth": hw, "headaxislength": hal})

        for arg in list(kwargs):
            if arg in quiver_kwargs:
                quiver_kwargs.update({arg: kwargs[arg]})
            else:
                scatter_kwargs.update({arg: kwargs[arg]})

        if (
            basis in adata.var_names
            and isinstance(color, str)
            and color in adata.layers.keys()
        ):
            c = interpret_colorkey(_adata, basis, color, perc)
        else:
            c = interpret_colorkey(_adata, color, layer, perc)

        if density is not None and 0 < density < 1:
            s = int(density * _adata.n_obs)
            ix_choice = np.random.choice(_adata.n_obs, size=s, replace=False)
            c = c[ix_choice] if len(c) == _adata.n_obs else c
            X = X[ix_choice]
            V = V[ix_choice]

        if projection == "3d" and X.shape[1] > 2 and V.shape[1] > 2:
            V, size = V / scale / 5, size / 10
            x0, x1, x2 = X[:, 0], X[:, 1], X[:, 2]
            v0, v1, v2 = V[:, 0], V[:, 1], V[:, 2]
            quiver3d_kwargs = {"zorder": 3, "linewidth": 0.5, "arrow_length_ratio": 0.3}
            c = list(c) + [element for element in list(c) for _ in range(2)]
            if is_color_like(c[0]):
                ax.quiver(x0, x1, x2, v0, v1, v2, color=c, **quiver3d_kwargs)
            else:
                ax.quiver(x0, x1, x2, v0, v1, v2, c, **quiver3d_kwargs)
        else:
            quiver_kwargs.update({"zorder": 3})
            if is_color_like(c[0]):
                ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=c, **quiver_kwargs)
            else:
                ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], c, **quiver_kwargs)

        scatter_kwargs.update({"basis": basis, "x": x, "y": y, "color": color})
        scatter_kwargs.update({"vkey": vkey, "layer": layer})
        ax = scatter(adata, size=size, title=title, ax=ax, zorder=0, **scatter_kwargs)

        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax
