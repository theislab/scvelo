from time import time
import scvelo as scv


def test_profile_dynamics_recovery_components():
    adata = scv.datasets.pancreas()
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=200)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

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

    start = time()
    for dm in dms:
        if dm.recoverable:
            dm.fit(assignment_mode="projection")
    print(f"Time fit: {time() - start}")


def test_profile_dynamics_recovery_details():
    adata = scv.datasets.simulation(random_seed=0, n_obs=300, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    var_names = scv.tools.dynamical_model.parse_var_names(
        "velocity_genes", adata=adata, use_raw=False, n_top_genes=None
    )

    # select a subset of genes
    genes = var_names[:100]

    t_create = 0
    t_t_and_alpha = 0
    t_scaling = 0
    t_rates = 0
    t_t = 0
    t_t_and_rates = 0
    t_update = 0
    t_t_and_rates_2 = 0
    t_div = 0
    t_lh = 0
    t_var = 0

    for gene in genes:
        print(gene)

        start = time()
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
        t_create += time() - start
        if not dm.recoverable:
            continue

        # fit() components one at a time

        start = time()
        dm.fit_t_and_alpha()
        t_t_and_alpha += time() - start

        start = time()
        dm.fit_scaling_()
        t_scaling += time() - start

        start = time()
        dm.fit_rates()
        t_rates += time() - start

        start = time()
        dm.fit_t_()
        t_t += time() - start

        start = time()
        dm.fit_t_and_rates()
        t_t_and_rates += time() - start

        start = time()
        dm.update(adjust_t_=False)
        t_update += time() - start

        start = time()
        dm.fit_t_and_rates(refit_time=False)
        t_t_and_rates_2 += time() - start

        start = time()
        dm.update()
        t_update += time() - start

        start = time()
        dm.tau, dm.tau_ = dm.get_divergence(mode="tau")
        t_div += time() - start

        start = time()
        dm.likelihood = dm.get_likelihood(refit_time=False)
        t_lh += time() - start

        start = time()
        dm.varx = dm.get_variance()
        t_var += time() - start

        print(dm.pars[:, -1])

    print("create", t_create)
    print("t and alpha", t_t_and_alpha)
    print("scaling", t_scaling)
    print("rates", t_rates)
    print("t", t_t)
    print("t and alpha", t_t_and_rates)
    print("update", t_update)
    print("t and rates 2", t_t_and_rates_2)
    print("div", t_div)
    print("lh", t_lh)
    print("var", t_var)
