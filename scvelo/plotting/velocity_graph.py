import warnings

import numpy as np
from scipy.sparse import csr_matrix, issparse

from scvelo import settings
from scvelo.preprocessing.neighbors import get_neighs
from scvelo.tools.transition_matrix import transition_matrix
from .docs import doc_params, doc_scatter
from .scatter import scatter
from .utils import (
    default_basis,
    default_size,
    get_basis,
    get_components,
    groups_to_bool,
    savefig_or_show,
)


@doc_params(scatter=doc_scatter)
def velocity_graph(
    adata,
    basis=None,
    vkey="velocity",
    which_graph=None,
    n_neighbors=10,
    arrows=None,
    arrowsize=3,
    alpha=0.8,
    perc=None,
    threshold=None,
    edge_width=0.2,
    edge_color="grey",
    edges_on_top=None,
    color=None,
    layer=None,
    size=None,
    groups=None,
    components=None,
    title=None,
    dpi=None,
    show=None,
    save=None,
    ax=None,
    **kwargs,
):
    """Plot of the velocity graph.

    Velocity graph with connectivities (dashed) and transitions (solid/arrows).

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    which_graph: `'velocity_graph'` or `'connectivities'`  (default: `None`)
        Whether to show transitions from velocity graph or neighbor connectivities.
    n_neighbors: `int` (default: 10)
        Number of neighbors to be included for generating connectivity / velocity graph.
    arrows: `bool` (default: `None`)
        Whether to display arrows instead of edges. Recommended to be used only on a
        cluster by setting groups parameter.
    arrowsize: `int` (default: 3)
        Size of the arrow heads.
    threshold: `float` (default: None)
        Threshold below which values are set to zero.
    edge_width: `float` (default: 0.2)
        Line width of edges.
    edge_color: `str` (default: "grey")
        Edge color. Can be a single color or a sequence of colors with the same length
        as edgelist. Color can be string or rgb (or rgba) tuple of floats from 0-1. If
        numeric values are specified they will be mapped to colors using the edge_cmap
        and edge_vmin,edge_vmax parameters.
    edges_on_top: `bool` (default: None)
        Whether or not to plot edges on top.

    {scatter}

    Returns
    -------
    If `show==False`, the `matplotlib.Axis` object.
    This is a test.
    """
    basis = default_basis(adata, **kwargs) if basis is None else get_basis(adata, basis)
    kwargs.update(
        {
            "basis": basis,
            "title": which_graph if title is None else title,
            "alpha": alpha,
            "components": components,
            "groups": groups,
            "dpi": dpi,
            "show": False,
            "save": None,
        }
    )
    ax = scatter(adata, layer=layer, color=color, size=size, ax=ax, zorder=0, **kwargs)

    from networkx import DiGraph, Graph

    if which_graph in {"neighbors", "connectivities"}:
        T = get_neighs(adata, "connectivities").copy()
        if perc is not None or threshold is not None:
            if threshold is None:
                threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    elif which_graph in adata.uns.keys():
        T = adata.uns[which_graph].copy()
        if perc is not None or threshold is not None:
            if threshold is None:
                threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    elif hasattr(adata, "obsp") and which_graph in adata.obsp.keys():
        T = adata.obsp[which_graph].copy()
        if perc is not None or threshold is not None:
            if threshold is None:
                threshold = np.percentile(T.data, perc)
            T.data[T.data < threshold] = 0
            T.eliminate_zeros()
    else:
        if threshold is None:
            threshold = 0.05
        T = transition_matrix(
            adata,
            vkey=vkey,
            weight_indirect_neighbors=0,
            n_neighbors=n_neighbors,
            perc=perc,
            threshold=threshold,
        )

    if groups is not None:
        if issparse(T):
            T = T.A
        T[~groups_to_bool(adata, groups, color)] = 0
        T = csr_matrix(T)
        T.eliminate_zeros()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X_emb = adata.obsm[f"X_{basis}"][:, get_components(components, basis)]
        node_size = (kwargs["size"] if "size" in kwargs else default_size(adata)) / 4
        edges = _draw_networkx_edges(
            DiGraph(T) if arrows else Graph(T),
            X_emb,
            node_size=node_size,
            width=edge_width,
            edge_color=edge_color,
            arrowsize=arrowsize,
            ax=ax,
        )
        if not arrows and not edges_on_top:
            edges.set_zorder(-2)
            edges.set_rasterized(settings._vector_friendly)

    savefig_or_show(dpi=dpi, save=save, show=show)
    if show is False:
        return ax


def _draw_networkx_edges(
    G,
    pos,
    edgelist=None,
    width=1.0,
    edge_color="k",
    style="solid",
    alpha=None,
    arrowstyle="-|>",
    arrowsize=3,
    edge_cmap=None,
    edge_vmin=None,
    edge_vmax=None,
    ax=None,
    arrows=True,
    label=None,
    node_size=300,
    nodelist=None,
    node_shape="o",
    connectionstyle=None,
    min_source_margin=0,
    min_target_margin=0,
):
    """Draw the edges of the graph G. Adjusted from networkx."""
    try:
        from numbers import Number

        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.colors import colorConverter, Colormap, Normalize
        from matplotlib.patches import FancyArrowPatch
    except ImportError:
        raise ImportError("Matplotlib required for draw()")
    except RuntimeError:
        print("Matplotlib unable to open display")
        raise

    if ax is None:
        ax = plt.gca()

    if edgelist is None:
        edgelist = list(G.edges())

    if not edgelist or len(edgelist) == 0:  # no edges!
        return None

    if nodelist is None:
        nodelist = list(G.nodes())

    # FancyArrowPatch handles color=None different from LineCollection
    if edge_color is None:
        edge_color = "k"

    # set edge positions
    edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])

    # Check if edge_color is an array of floats and map to edge_cmap.
    # This is the only case handled differently from matplotlib
    if (
        np.iterable(edge_color)
        and (len(edge_color) == len(edge_pos))
        and np.alltrue([isinstance(c, Number) for c in edge_color])
    ):
        if edge_cmap is not None:
            assert isinstance(edge_cmap, Colormap)
        else:
            edge_cmap = plt.get_cmap()
        if edge_vmin is None:
            edge_vmin = min(edge_color)
        if edge_vmax is None:
            edge_vmax = max(edge_color)
        color_normal = Normalize(vmin=edge_vmin, vmax=edge_vmax)
        edge_color = [edge_cmap(color_normal(e)) for e in edge_color]

    if not G.is_directed() or not arrows:
        edge_collection = LineCollection(
            edge_pos,
            colors=edge_color,
            linewidths=width,
            antialiaseds=(1,),
            linestyle=style,
            alpha=alpha,
        )

        edge_collection.set_cmap(edge_cmap)
        edge_collection.set_clim(edge_vmin, edge_vmax)

        edge_collection.set_zorder(1)  # edges go behind nodes
        edge_collection.set_label(label)
        ax.add_collection(edge_collection)

        return edge_collection

    arrow_collection = None

    if G.is_directed() and arrows:
        # Note: Waiting for someone to implement arrow to intersection with
        # marker.  Meanwhile, this works well for polygons with more than 4
        # sides and circle.

        def to_marker_edge(marker_size, marker):
            if marker in "s^>v<d":  # `large` markers need extra space
                return np.sqrt(2 * marker_size) / 2
            else:
                return np.sqrt(marker_size) / 2

        # Draw arrows with `matplotlib.patches.FancyarrowPatch`
        arrow_collection = []
        mutation_scale = arrowsize  # scale factor of arrow head

        # FancyArrowPatch doesn't handle color strings
        arrow_colors = colorConverter.to_rgba_array(edge_color, alpha)
        for i, (src, dst) in enumerate(edge_pos):
            x1, y1 = src
            x2, y2 = dst
            shrink_source = 0  # space from source to tail
            shrink_target = 0  # space from  head to target
            if np.iterable(node_size):  # many node sizes
                source, target = edgelist[i][:2]
                source_node_size = node_size[nodelist.index(source)]
                target_node_size = node_size[nodelist.index(target)]
                shrink_source = to_marker_edge(source_node_size, node_shape)
                shrink_target = to_marker_edge(target_node_size, node_shape)
            else:
                shrink_source = shrink_target = to_marker_edge(node_size, node_shape)

            if shrink_source < min_source_margin:
                shrink_source = min_source_margin

            if shrink_target < min_target_margin:
                shrink_target = min_target_margin

            if len(arrow_colors) == len(edge_pos):
                arrow_color = arrow_colors[i]
            elif len(arrow_colors) == 1:
                arrow_color = arrow_colors[0]
            else:  # Cycle through colors
                arrow_color = arrow_colors[i % len(arrow_colors)]

            if np.iterable(width):
                if len(width) == len(edge_pos):
                    line_width = width[i]
                else:
                    line_width = width[i % len(width)]
            else:
                line_width = width

            arrow = FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle=arrowstyle,
                shrinkA=shrink_source,
                shrinkB=shrink_target,
                mutation_scale=mutation_scale,
                color=arrow_color,
                linewidth=line_width,
                connectionstyle=connectionstyle,
                linestyle=style,
                zorder=1,
            )  # arrows go behind nodes

            # There seems to be a bug in matplotlib to make collections of
            # FancyArrowPatch instances. Until fixed, the patches are added
            # individually to the axes instance.
            arrow_collection.append(arrow)
            ax.add_patch(arrow)

    # update view
    minx = np.amin(np.ravel(edge_pos[:, :, 0]))
    maxx = np.amax(np.ravel(edge_pos[:, :, 0]))
    miny = np.amin(np.ravel(edge_pos[:, :, 1]))
    maxy = np.amax(np.ravel(edge_pos[:, :, 1]))

    w = maxx - minx
    h = maxy - miny
    padx, pady = 0.05 * w, 0.05 * h
    corners = (minx - padx, miny - pady), (maxx + padx, maxy + pady)
    ax.update_datalim(corners)
    ax.autoscale_view()

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return arrow_collection
