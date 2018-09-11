from .utils import interpret_colorkey, get_components, plot_colorbar, savefig
import matplotlib.pyplot as pl
import scanpy.api.pl as scpl
from matplotlib.ticker import MaxNLocator


def scatter(adata, x=None, y=None, basis='umap', layer=None, color=None, xlabel=None, ylabel=None, color_map=None,
            perc=None, size=5, alpha=1, fontsize=None, colorbar=False, groups=None, use_raw=None, sort_order=True,
            legend_loc='none', legend_fontsize=None, legend_fontweight=None, projection='2d', palette=None,
            right_margin=None, left_margin=None, components=None, frameon=False, title=None, show=True, figsize=(7,5),
            dpi=80, save=None, ax=None, zorder=0, **kwargs):
    """Scatter plot along observations or variables axes.
    Color the plot using annotations of observations (`.obs`), variables (`.var`) or expression of genes (`.var_names`).

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    basis: `str` (default='tsne')
        plots embedding obsm['X_' + basis]

    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate

    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate

    color : `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """
    colors = color if isinstance(color, (list, tuple)) else [color]
    layers = layer if isinstance(layer, (list, tuple)) else [layer]

    if len(colors) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0]*len(colors), figsize[1]), dpi=dpi))):
            scatter(adata, basis=basis, layer=layer, color=color[i], xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                    perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                    colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=pl.subplot(gs),
                    use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                    legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                    palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(layers) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            scatter(adata, basis=basis, layer=layers[i], color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                    perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                    colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=pl.subplot(gs),
                    use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                    legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                    palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    else:
        if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca()
        if color_map is None: color_map = 'viridis_r' if (color == 'root' or color == 'end') else 'RdBu_r'
        is_embedding = (x is None) | (y is None)

        if color in adata.obs and adata.obs[color].dtype.name == 'category' and is_embedding:
            ax = scpl.scatter(adata, color=color, use_raw=use_raw, sort_order=sort_order, alpha=alpha, basis=basis,
                              groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                              legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, color_map=color_map,
                              palette=palette, right_margin=right_margin, left_margin=left_margin,
                              size=size, title=title, frameon=frameon, show=False, save=None, ax=ax, **kwargs)
        else:
            if is_embedding:
                X_emb = adata.obsm['X_' + basis][:, get_components(components)]
                x, y = X_emb[:, 0], X_emb[:, 1]
                ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
            else:
                ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
                labelsize = int(fontsize * .75) if fontsize is not None else None
                ax.tick_params(axis='both', which='major', labelsize=labelsize)
            c = interpret_colorkey(adata, color, layer, perc)
            pl.scatter(x, y, c=c, cmap=color_map, s=size, alpha=alpha, zorder=zorder, **kwargs)

            if isinstance(xlabel, str) and isinstance(ylabel, str):
                pl.xlabel(xlabel, fontsize=fontsize)
                pl.ylabel(ylabel, fontsize=fontsize)
            elif basis is not None:
                component_name = ('DC' if basis == 'diffmap' else 'tSNE' if basis == 'tsne' else 'UMAP' if basis == 'umap'
                else 'PC' if basis == 'pca' else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis else basis)
                pl.xlabel(component_name + '1')
                pl.ylabel(component_name + '2')

            if isinstance(title, str): pl.title(title, fontsize=fontsize)
            elif isinstance(layer, str) and isinstance(color, str): pl.title(color + ' ' + layer, fontsize=fontsize)
            elif isinstance(color, str): pl.title(color, fontsize=fontsize)

            if not frameon: pl.axis('off')
            if colorbar: plot_colorbar(ax)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)

        if show: pl.show()
        else: return ax