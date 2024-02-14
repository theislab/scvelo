import contextlib

from scanpy.tools import diffmap, dpt, louvain, tsne, umap

from ._em_model import ExpectationMaximizationModel
from ._em_model_core import (
    align_dynamics,
    differential_kinetic_test,
    DynamicsRecovery,
    latent_time,
    rank_dynamical_genes,
    recover_dynamics,
    recover_latent_time,
)
from ._steady_state_model import SecondOrderSteadyStateModel, SteadyStateModel
from .paga import paga
from .rank_velocity_genes import rank_velocity_genes, velocity_clusters
from .score_genes_cell_cycle import score_genes_cell_cycle
from .terminal_states import eigs, terminal_states
from .transition_matrix import transition_matrix
from .velocity import velocity, velocity_genes
from .velocity_confidence import velocity_confidence, velocity_confidence_transition
from .velocity_embedding import velocity_embedding
from .velocity_graph import velocity_graph
from .velocity_pseudotime import velocity_map, velocity_pseudotime

with contextlib.suppress(ImportError):
    from ._vi_model import VELOVI


__all__ = [
    "align_dynamics",
    "differential_kinetic_test",
    "diffmap",
    "dpt",
    "DynamicsRecovery",
    "eigs",
    "latent_time",
    "louvain",
    "paga",
    "rank_dynamical_genes",
    "rank_velocity_genes",
    "recover_dynamics",
    "recover_latent_time",
    "score_genes_cell_cycle",
    "terminal_states",
    "transition_matrix",
    "tsne",
    "umap",
    "velocity",
    "velocity_clusters",
    "velocity_confidence",
    "velocity_confidence_transition",
    "velocity_embedding",
    "velocity_genes",
    "velocity_graph",
    "velocity_map",
    "velocity_pseudotime",
    "SteadyStateModel",
    "SecondOrderSteadyStateModel",
    "ExpectationMaximizationModel",
]
if "VELOVI" in locals():
    __all__ += ["VELOVI"]

del contextlib
