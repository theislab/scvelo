import sys
import cloudpickle as pickle
from time import time
import scvelo as scv


def test_todo():
    adata = scv.datasets.pancreas()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=200)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print(f"Size adata: {sys.getsizeof(adata)}")
    print(f"Size adata dump: {sys.getsizeof(pickle.dumps(adata))}")

    var_names = scv.tools.dynamical_model.parse_var_names(
        "velocity_genes", adata=adata, use_raw=False, n_top_genes=None
    )

    start = time()
    dms = []
    for i, gene in enumerate(var_names):
        dm = scv.tools.DynamicsRecovery(
            adata,
            gene,
            use_raw=False,
            load_pars=None,
            max_iter=10,
            fit_time=True,
            fit_steady_states=True,
            fit_connected_states=None,
            fit_scaling=True,
            fit_basal_transcription=None,
            steady_state_prior=None,
        )
        dms.append(dm)
    print(f"Time create: {time() - start}")
    print(f"Size dms: {sys.getsizeof(dms)}")
    print(f"Size dms dump: {sys.getsizeof([pickle.dumps(dm) for dm in dms])}")

    start = time()
    for dm in dms:
        if dm.recoverable:
            dm.fit(assignment_mode="projection")
    print(f"Time fit: {time() - start}")
