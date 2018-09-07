from .utils import interpret_colorkey
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def scatter(adata, x=None, y=None, basis='umap', layer=None, color=None, xlabel=None, ylabel=None, color_map=None,
            perc=None, size=5, alpha=1, fontsize=None, frameon=False, title=None, show=True, colorbar=False,
            figsize=(7,5), dpi=80, save=None, ax=None, **kwargs):
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
    if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca()

    if (x is None) | (y is None):
        X_emb = adata.obsm['X_' + basis][:, :2]
        x, y = X_emb[:, 0], X_emb[:, 1]
        ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        labelsize = int(fontsize * .75) if fontsize is not None else None
        ax.tick_params(axis='both', which='major', labelsize=labelsize)

    pl.scatter(x, y, c=interpret_colorkey(adata, color, layer, perc), cmap=color_map, s=size, alpha=alpha, **kwargs)

    if isinstance(title, str):
        pl.title(title, fontsize=fontsize)

    if isinstance(xlabel, str) and isinstance(ylabel, str):
        pl.xlabel(xlabel, fontsize=fontsize)
        pl.ylabel(ylabel, fontsize=fontsize)
    else:
        pl.xlabel(basis + '1')
        pl.ylabel(basis + '2')

    if not frameon: pl.axis('off')
    if isinstance(save, str): pl.savefig(save)

    if colorbar:
        cb = pl.colorbar(orientation='vertical', cax=inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0))
        cb.locator = (MaxNLocator(nbins=3))
        cb.update_ticks()

    if show: pl.show()
    else: return ax