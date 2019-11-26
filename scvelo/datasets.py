"""Builtin Datasets.
"""
from .read_load import read, load
from .preprocessing.utils import cleanup
from anndata import AnnData
import numpy as np
import pandas as pd


def toy_data(n_obs=None):
    """
    Randomly sampled from the Dentate Gyrus dataset.

    Arguments
    ---------
    n_obs: `int` (default: `None`)
        Size of the sampled dataset

    Returns
    -------
    Returns `adata` object
    """

    """Random samples from Dentate Gyrus.
    """
    adata_dg = dentategyrus()

    if n_obs is not None:
        indices = np.random.choice(adata_dg.n_obs, n_obs)
        adata = adata_dg[indices]
    else:
        adata = adata_dg
    adata.obs_names_make_unique()
    return adata.copy()


def dentategyrus(adjusted=True):
    """Dentate Gyrus dataset from `Hochgerner et al. (2018) <https://doi.org/10.1038/s41593-017-0056-2>`_.

    Dentate gyrus is part of the hippocampus involved in learning, episodic memory formation and spatial coding.
    It is measured using 10X Genomics Chromium and described in Hochgerner et al. (2018).
    The data consists of 25,919 genes across 3,396 cells and provides several interesting characteristics.

    Returns
    -------
    Returns `adata` object
    """

    if adjusted:
        filename = 'data/DentateGyrus/10X43_1.h5ad'
        url = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/10X43_1.h5ad'
        adata = read(filename, backup_url=url, sparse=True, cache=True)

    else:
        filename = 'data/DentateGyrus/10X43_1.loom'
        url = 'http://pklab.med.harvard.edu/velocyto/DG1/10X43_1.loom'
        adata = read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
        cleanup(adata, clean='all', keep={'spliced', 'unspliced', 'ambiguous'})

        url_louvain = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/DG_clusters.npy'
        url_umap = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/DG_umap.npy'

        adata.obs['clusters'] = load('./data/DentateGyrus/DG_clusters.npy', url_louvain, allow_pickle=True)
        adata.obsm['X_umap'] = load('./data/DentateGyrus/DG_umap.npy', url_umap, allow_pickle=True)

        adata.obs['clusters'] = pd.Categorical(adata.obs['clusters'])

    return adata


def forebrain():
    """Developing human forebrain.

    Forebrain tissue of a week 10 embryo, focusing on the glutamatergic neuronal lineage.

    Returns
    -------
    Returns `adata` object
    """
    filename = 'data/ForebrainGlut/hgForebrainGlut.loom'
    url = 'http://pklab.med.harvard.edu/velocyto/hgForebrainGlut/hgForebrainGlut.loom'
    adata = read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def pancreatic_endocrinogenesis():
    """Pancreatic endocrinogenesis from `Bastidas-Ponce et al. (2019) <https://doi.org/10.1242/dev.173849>`_.

    Pancreatic epithelial and Ngn3-Venus fusion (NVF) cells during secondary transition / embryonic day 15.5.

    Returns
    -------
    Returns `adata` object
    """
    filename = 'data/Pancreas/endocrinogenesis_day15.h5ad'
    url = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad'
    adata = read(filename, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def simulation(n_obs=300, n_vars=None, alpha=None, beta=None, gamma=None, alpha_=None, t_max=None,
               noise_model='normal', noise_level=1, switches=None, random_seed=0):
    """Simulation of mRNA metabolism with transcription, splicing and degradation

    Returns
    -------
    Returns `adata` object
    """
    from .tools.dynamical_model_utils import unspliced, spliced, vectorize, mRNA
    np.random.seed(random_seed)

    def draw_poisson(n):
        from random import uniform, seed  # draw from poisson
        seed(random_seed)
        t = np.cumsum([- .1 * np.log(uniform(0, 1)) for _ in range(n - 1)])
        return np.insert(t, 0, 0)  # prepend t0=0

    def simulate_dynamics(tau, alpha, beta, gamma, u0, s0, noise_model, noise_level):
        ut, st = mRNA(tau, u0, s0, alpha, beta, gamma)
        if noise_model is 'normal':  # add noise
            ut += np.random.normal(scale=noise_level * np.percentile(ut, 99) / 10, size=len(ut))
            st += np.random.normal(scale=noise_level * np.percentile(st, 99) / 10, size=len(st))
        ut, st = np.clip(ut, 0, None), np.clip(st, 0, None)
        return ut, st

    def simulate_gillespie(alpha, beta, gamma):
        # update rules: transcription (u+1,s), splicing (u-1,s+1), degradation (u,s-1), nothing (u,s)
        update_rule = np.array([[1, 0], [-1, 1], [0, -1], [0, 0]])

        def update(props):
            if np.sum(props) > 0:
                props /= np.sum(props)
            p_cumsum = props.cumsum()
            p = np.random.rand()
            i = 0
            while p > p_cumsum[i]:
                i += 1
            return update_rule[i]

        u, s = np.zeros(len(alpha)), np.zeros(len(alpha))
        for i, alpha_i in enumerate(alpha):
            u_, s_ = (u[i - 1], s[i - 1]) if i > 0 else (0, 0)
            du, ds = update(props=np.array([alpha_i, beta * u_, gamma * s_]))
            u[i], s[i] = (u_ + du, s_ + ds)
        return u, s

    alpha = 5 if alpha is None else alpha
    beta = .5 if beta is None else beta
    gamma = .3 if gamma is None else gamma
    alpha_ = 0 if alpha_ is None else alpha_

    t = draw_poisson(n_obs)
    if t_max is not None:
        t *= t_max / np.max(t)
    t_max = np.max(t)

    def cycle(array, n_vars=None):
        if isinstance(array, (np.ndarray, list, tuple)):
            return array if n_vars is None else array * int(np.ceil(n_vars / len(array)))
        else:
            return [array] if n_vars is None else [array] * n_vars

    # switching time point obtained as fraction of t_max rounded down
    switches = cycle([.4, .7, 1, .1], n_vars) if switches is None else cycle(switches, n_vars)
    t_ = np.array([np.max(t[t < t_i * t_max]) for t_i in switches])

    noise_level = cycle(noise_level, len(switches) if n_vars is None else n_vars)

    n_vars = max(len(switches), len(noise_level)) if n_vars is None else n_vars
    U = np.zeros(shape=(len(t), n_vars))
    S = np.zeros(shape=(len(t), n_vars))

    for i in range(n_vars):
        alpha_i = alpha[i] if isinstance(alpha, (tuple, list, np.ndarray)) else alpha
        beta_i = beta[i] if isinstance(beta, (tuple, list, np.ndarray)) else beta
        gamma_i = gamma[i] if isinstance(gamma, (tuple, list, np.ndarray)) else gamma
        tau, alpha_vec, u0_vec, s0_vec = vectorize(t, t_[i], alpha_i, beta_i, gamma_i, alpha_=alpha_, u0=0, s0=0)

        if noise_model is 'gillespie':
            U[:, i], S[:, i] = simulate_gillespie(alpha_vec, beta, gamma)
        else:
            U[:, i], S[:, i] = simulate_dynamics(tau, alpha_vec, beta_i, gamma_i,
                                                 u0_vec, s0_vec, noise_model, noise_level[i])

    obs = {'true_t': t.round(2)}
    var = {'true_t_': t_[:n_vars],
           'true_alpha': np.ones(n_vars) * alpha,
           'true_beta': np.ones(n_vars) * beta,
           'true_gamma': np.ones(n_vars) * gamma,
           'true_scaling': np.ones(n_vars)}
    layers = {'unspliced': U, 'spliced': S}

    return AnnData(S, obs, var, layers=layers)
