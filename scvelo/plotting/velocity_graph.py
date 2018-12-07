from .. import settings
from ..tools.transition_matrix import transition_matrix
from .utils import savefig, default_basis
from .scatter import scatter
from .docs import doc_scatter, doc_params

import warnings
import numpy as np
import matplotlib.pyplot as pl


@doc_params(scatter=doc_scatter)
def velocity_graph(adata, basis=None, vkey='velocity', which_graph='velocity', n_neighbors=10,
                   alpha=.8, perc=90, edge_width=.2, edge_color='grey', color=None, use_raw=None, layer=None,
                   color_map=None, colorbar=False, palette=None, size=None,  sort_order=True, groups=None,
                   components=None, projection='2d', legend_loc='on data', legend_fontsize=None, legend_fontweight=None,
                   right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None, fontsize=None,
                   figsize=None, dpi=None, frameon=None, show=True, save=None, ax=None):
    """\
    Plot of the velocity graph.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    which_graph: `'velocity'` or `'neighbors'` (default: `'velocity'`)
        Whether to show transitions from velocity graph or connectivities from neighbors graph.

    {scatter}

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    basis = default_basis(adata) if basis is None else basis
    title = which_graph + ' graph' if title is None else title
    scatter_kwargs = {"basis": basis, "perc": perc, "use_raw": use_raw, "sort_order": sort_order, "alpha": alpha,
                      "components": components, "projection": projection, "legend_loc": legend_loc, "groups": groups,
                      "legend_fontsize": legend_fontsize, "legend_fontweight": legend_fontweight, "palette": palette,
                      "color_map": color_map, "frameon": frameon, "title": title, "xlabel": xlabel, "ylabel": ylabel,
                      "right_margin": right_margin, "left_margin": left_margin, "colorbar": colorbar, "dpi": dpi,
                      "fontsize": fontsize, "show": False, "save": None, "figsize": figsize, }
    ax = scatter(adata, layer=layer, color=color, size=size, ax=ax, zorder=0, **scatter_kwargs)

    from networkx import Graph, draw_networkx_edges
    if which_graph == 'neighbors':
        T = adata.uns['neighbors']['connectivities']
        if perc is not None:
            threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    else:
        T = transition_matrix(adata, vkey=vkey, weight_indirect_neighbors=0, n_neighbors=n_neighbors, perc=perc)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        edges = draw_networkx_edges(Graph(T), adata.obsm['X_' + basis], width=edge_width, edge_color=edge_color, ax=ax)
        edges.set_zorder(-2)
        edges.set_rasterized(settings._vector_friendly)

    if isinstance(save, str):
        savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
    if show: pl.show()
    else: return ax
