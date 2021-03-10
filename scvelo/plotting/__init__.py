from scanpy.plotting import paga_compare, rank_genes_groups

from .gridspec import gridspec
from .heatmap import heatmap
from .paga import paga
from .proportions import proportions
from .scatter import diffmap, draw_graph, pca, phate, scatter, tsne, umap
from .simulation import simulation
from .summary import summary
from .utils import hist, plot
from .velocity import velocity
from .velocity_embedding import velocity_embedding
from .velocity_embedding_grid import velocity_embedding_grid
from .velocity_embedding_stream import velocity_embedding_stream
from .velocity_graph import velocity_graph

__all__ = [
    "diffmap",
    "draw_graph",
    "gridspec",
    "heatmap",
    "hist",
    "paga",
    "paga_compare",
    "pca",
    "phate",
    "plot",
    "proportions",
    "rank_genes_groups",
    "scatter",
    "simulation",
    "summary",
    "tsne",
    "umap",
    "velocity",
    "velocity_embedding",
    "velocity_embedding_grid",
    "velocity_embedding_stream",
    "velocity_graph",
]
