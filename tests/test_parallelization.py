import scvelo as scv
from time import time

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)


def test_pancreas():
    adata = scv.datasets.pancreas()

    start = time()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=200)
    # actually 2000 genes
    print(f"Time filter and normalization: {time() - start}")

    start = time()
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    print(f"Time moments: {time() - start}")

    start = time()
    scv.tl.recover_dynamics(adata, n_procs=2)
    print(f"Time recover_dynamics: {time() - start}")
