from ..tools.utils import groups_to_bool
from .utils import default_basis, default_size, savefig_or_show, \
    default_color, make_unique_list, make_unique_valid_list, get_components
from .scatter import scatter
from .docs import doc_scatter, doc_params

from matplotlib import rcParams
import matplotlib.pyplot as pl
import numpy as np

from scanpy.plotting._tools.paga import paga as scanpy_paga


@doc_params(scatter=doc_scatter)
def paga(adata, basis=None, vkey='velocity', color=None, layer=None, title=None, threshold=0.2, layout=None,
         layout_kwds={}, init_pos=None, root=0, labels=None, single_component=False, solid_edges='connectivities',
         dashed_edges='connectivities', transitions='transitions_confidence', node_size_scale=1, node_size_power=0.5,
         edge_width_scale=.4, min_edge_width=None, max_edge_width=1, arrowsize=15, random_state=0, pos=None,
         normalize_to_color=False, cmap=None, cax=None, cb_kwds={}, add_pos=True, export_to_gexf=False, plot=True,
         use_raw=None, size=None, groups=None, components=None, figsize=None, dpi=None, show=True, save=None, ax=None,
         ncols=None, scatter_flag=None, **kwargs):
    """\
    PAGA plot on the embedding.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate
    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    {scatter}
    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """

    if scatter_flag is None:
        scatter_flag = ax is None
    vkey = [key for key in adata.layers.keys() if 'velocity' in key and '_u' not in key] if vkey is 'all' else vkey
    layers, vkeys, colors = make_unique_list(layer), make_unique_list(vkey), make_unique_list(color, allow_array=True)
    bases = [default_basis(adata) if basis is None else basis for basis in make_unique_valid_list(adata, basis)]

    paga_kwargs = {'threshold': threshold, 'layout': layout, 'layout_kwds': layout_kwds, 'init_pos': init_pos,
                   'root': root, 'labels': labels, 'single_component': single_component,
                   'solid_edges': solid_edges, 'dashed_edges': dashed_edges, 'transitions': transitions,
                   'node_size_scale': node_size_scale, 'node_size_power': node_size_power,
                   'edge_width_scale': edge_width_scale, 'min_edge_width': min_edge_width,
                   'max_edge_width': max_edge_width, 'arrowsize': arrowsize, 'random_state': random_state,
                   'pos': pos, 'normalize_to_color': normalize_to_color, 'cmap': cmap, 'cax': cax, 'cb_kwds': cb_kwds,
                   'add_pos': add_pos, 'export_to_gexf': export_to_gexf, 'colors': colors, 'plot': plot}

    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 \
        else vkeys if len(vkeys) > 1 else bases if len(bases) > 1 else None
    if multikey is not None:
        if title is None:
            title = list(multikey)
        elif isinstance(title, (list, tuple)):
            title *= int(np.ceil(len(multikey) / len(title)))
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        ax = []
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))):
            if i < len(multikey):
                ax.append(paga(adata, size=size, ax=pl.subplot(gs), scatter_flag=scatter_flag,
                               basis=bases[i] if len(bases) > 1 else basis,
                               color=colors[i] if len(colors) > 1 else color,
                               layer=layers[i] if len(layers) > 1 else layer,
                               vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                               title=title[i] if isinstance(title, (list, tuple)) else title,
                               **kwargs, **paga_kwargs))
        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax

    else:

        color, layer, vkey, basis = colors[0], layers[0], vkeys[0], basis
        color = default_color(adata) if color is None else color
        size = default_size(adata) / 2 if size is None else size
        _adata = adata[
            groups_to_bool(adata, groups, groupby=color)] if groups is not None and color in adata.obs.keys() else adata

        if basis in adata.var_names and basis is not None:
            x = adata[:, basis].layers['spliced'] if use_raw else adata[:, basis].layers['Ms']
            y = adata[:, basis].layers['unspliced'] if use_raw else adata[:, basis].layers['Mu']
        elif basis is not None:
            X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)]
            x, y = X_emb[:, 0], X_emb[:, 1]

        if basis is None and pos is None:
            pos = None  # default to paga embedding
        elif pos is None:
            if 'paga' in adata.uns:
                # Recompute the centroid positions
                categories = list(adata.obs[color].cat.categories)
                pos = np.zeros((len(categories), 2))
                for ilabel, label in enumerate(categories):
                    X_emb = adata.obsm['X_' + basis][adata.obs[color] == label, :2]
                    x_pos, y_pos = np.median(X_emb, axis=0)
                    pos[ilabel] = [x_pos, y_pos]
            else:
                raise ValueError(
                    'You need to run `scv.tl.paga` first.')
        paga_kwargs['pos'] = pos

        ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax
        if scatter_flag and basis is not None:
            if 'alpha' not in kwargs: kwargs['alpha'] = .5
            ax = scatter(adata, basis=basis, x=x, y=y, vkey=vkey, layer=layer, color=color, size=size, title=title,
                         ax=ax, save=None, zorder=0, show=False, **kwargs)
        scanpy_paga(adata, ax=ax, show=False, text_kwds={'alpha': 0}, **paga_kwargs)

        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax
