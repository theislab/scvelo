from ..preprocessing import filter_and_normalize, moments
from . import velocity, velocity_graph, velocity_embedding


def run(data, basis=None, mode='deterministic', min_counts=10, n_pcs=30, n_neighbors=30, copy=False):
    from time import time
    start = time()

    adata = data.copy() if copy else data
    filter_and_normalize(adata, min_counts=min_counts)
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    velocity(adata, mode=mode)
    velocity_graph(adata)
    if basis is not None: velocity_embedding(adata, basis=basis)

    print('/n Total time (seconds): ' + str(round(time() - start, 2)))

    return adata if copy else None