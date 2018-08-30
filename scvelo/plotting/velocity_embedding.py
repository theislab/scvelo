from .utils import *
from scanpy.api.pl import scatter
from matplotlib.colors import is_color_like


def velocity_embedding(adata, basis='umap', vbasis='velocity', layer=None, density=1, color=None,
                       use_raw=True, sort_order=True, alpha=.2, groups=None, components=None, projection='2d',
                       legend_loc='right margin', legend_fontsize=None, legend_fontweight=None,
                       color_map=None, palette=None, frameon=False, right_margin=None, left_margin=None,
                       size=None, title=None, show=True, save=None, ax=None, **kwargs):
    """Scatter plot with velocities along `.obs` or `.var` axes.
    Color the plot using annotations of observations (`.obs`), variables (`.var`) or expression of genes (`.var_names`).

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'umap'`)
        Key for embedding coordinates.
    vbasis: `str` (default: `'velocity'`)
        Key for velocity embedding coordinates.
    color : `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    if ax is None: ax = pl.figure(None, (14, 10), dpi=160).gca()

    scatter(adata, color=color, use_raw=use_raw, sort_order=sort_order, alpha=alpha, basis=basis,
            groups=groups, components=components, projection=projection, legend_loc=legend_loc,
            legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, color_map=color_map,
            palette=palette, frameon=frameon, right_margin=right_margin, left_margin=left_margin,
            size=size, title=title, show=False, save=save, ax=ax)

    vbasis += '_' + basis
    if vbasis not in adata.obsm_keys():
        raise ValueError(
            'You need to run `tl.velocity_embedding` first to compute embedded velocity vectors.')

    X_emb = adata.obsm['X_' + basis][:, :2]
    V_emb = adata.obsm[vbasis]

    _kwargs = {"scale": 1, "width": .0005, "edgecolors": 'k', "headwidth": 9, "headlength": 10, "headaxislength": 6, "linewidth": .25}
    _kwargs.update(kwargs)

    ix_choice = np.random.choice(adata.n_obs, size=int(density * adata.n_obs), replace=False)

    X, V, C = X_emb[ix_choice], V_emb[ix_choice], interpret_colorkey(adata, color, layer, perc=[2, 98])[ix_choice]

    if is_color_like(C[0]):
        pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=C, angles='xy', scale_units='xy', cmap=color_map, **_kwargs)
    else:
        pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], C, angles='xy', scale_units='xy', cmap=color_map, **_kwargs)

    if show: pl.show()
    else: return ax