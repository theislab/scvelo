from .. import settings
from ..tools.transition_matrix import transition_matrix
from .utils import savefig_or_show, default_basis, get_components, get_basis, groups_to_bool
from .scatter import scatter
from .docs import doc_scatter, doc_params

import warnings
import numpy as np
from scipy.sparse import issparse, csr_matrix


@doc_params(scatter=doc_scatter)
def velocity_graph(adata, basis=None, vkey='velocity', which_graph='velocity', n_neighbors=10, arrows=None,
                   alpha=.8, perc=90, edge_width=.2, edge_color='grey', edges_on_top=None, color=None, layer=None,
                   size=None, groups=None, components=None, title=None, dpi=None, show=True, save=None, ax=None, **kwargs):
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
    n_neighbors: `int` (default: 10)
        Number of neighbors to be included for generating connectivity / velocity graph.
    arrows: `bool` (default: `None`)
        Whether to display arrows instead of edges. Recommended to be used only on a cluster by setting groups parameter.

    {scatter}

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    basis = default_basis(adata) if basis is None else get_basis(adata, basis)
    kwargs.update({"basis": basis, "title": which_graph + ' graph' if title is None else title,
                   "alpha": alpha, "components": components, "groups": groups, "dpi": dpi, "show": False, "save": None})
    ax = scatter(adata, layer=layer, color=color, size=size, ax=ax, zorder=0, **kwargs)

    from networkx import Graph, DiGraph, draw_networkx_edges
    if which_graph in {'neighbors', 'connectivities'}:
        T = adata.uns['neighbors']['connectivities'].copy()
        if perc is not None:
            threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    elif which_graph in adata.uns.keys():
        T = adata.uns[which_graph].copy()
        if perc is not None:
            threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    else:
        T = transition_matrix(adata, vkey=vkey, weight_indirect_neighbors=0, n_neighbors=n_neighbors, perc=perc)

    if groups is not None:
        if issparse(T): T = T.A
        T[~groups_to_bool(adata, groups, color)] = 0
        T = csr_matrix(T)
        T.eliminate_zeros()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)]
        if arrows:
            edges = draw_networkx_edges(DiGraph(T), X_emb, width=edge_width, edge_color=edge_color, ax=ax)
        else:
            edges = draw_networkx_edges(Graph(T), X_emb, width=edge_width, edge_color=edge_color, ax=ax)
            if not edges_on_top:
                edges.set_zorder(-2)
                edges.set_rasterized(settings._vector_friendly)

    savefig_or_show(dpi=dpi, save=save, show=show)
    if not show: return ax
