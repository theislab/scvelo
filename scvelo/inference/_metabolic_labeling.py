from typing import Dict, List, Optional

from joblib import delayed, Parallel
from tqdm import tqdm

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, issparse
from sklearn.metrics import pairwise_distances

from anndata import AnnData


def _get_bounds(experiment: np.ndarray):
    """Returns bounds for constraint optimization problem.

    Arguments:
    ---------
    experiment
        Type of experiment (``"chase"`` or ``"pulse"``) that observation was generated from.

    Returns
    -------
    Bounds for constraint optimization problem.
    """
    if {"chase", "both"}.intersection(experiment):
        return ([1e-7, None], [1e-7, None], [1e-7, None])
    else:
        return ([1e-7, None], [1e-7, None])


def _get_chase_estimate(alpha, gamma, r0, labeling_time):
    return r0 - alpha / gamma * (1 - np.exp(-gamma * labeling_time))


def _get_pulse_estimate(alpha, gamma, labeling_time):
    return alpha / gamma * (1 - np.exp(-gamma * labeling_time))


def get_mse(
    x, measured_new: np.ndarray, labeling_time: np.ndarray, experiment: np.ndarray
):
    """Calculates mean squared error (MSE) between predicted and measurmed, newly synthesized mRNA.

    Arguments:
    ---------
    x
        Parameter estimate.
    measured_new
        Measured mRNA.
    labeling_time
        Labeling time of each observation.
    experiment
        Type of experiment (``"chase"`` or ``"pulse"``) that observation was generated from.

    Returns
    -------
        Mean squared error between measurement and prediction.
    """
    alpha = x[0]
    gamma = x[1]
    if {"chase", "both"}.intersection(experiment):
        r0 = x[2]

    estimated_new = np.full(labeling_time.shape, np.nan)

    pulse_mask = experiment == "pulse"
    estimated_new[pulse_mask] = _get_pulse_estimate(
        alpha=alpha, gamma=gamma, labeling_time=labeling_time[pulse_mask]
    )
    if {"chase", "both"}.intersection(experiment):
        estimated_new[~pulse_mask] = _get_chase_estimate(
            alpha=alpha, gamma=gamma, r0=r0, labeling_time=labeling_time[~pulse_mask]
        )

    return np.mean((measured_new - estimated_new) ** 2)


def _get_n_neighbors(
    X,
    labeling_time_mask,
    dist_argsort: Dict[float, np.ndarray],
    n_obs: Dict[float, int],
    sparse_op: bool,
    n_nontrivial_counts: int,
) -> pd.DataFrame:
    """Get number of neighbors required to include ``n_nontrivial_counts`` counts per labeling time.

    Arguments:
    ---------
    X
        Gene expression.
    labeling_times
        Array of labeling times.
    dist_argsort
        Argsorted distances for each labeling time point.
    n_obs
        Number of observations per labeling_time
    sparse_op
        Boolean flag to run operations on sparse or dense matrix.
    n_nontrivial_counts
        Number of non-trivial counts to consider for each labeling time point.

    Returns
    -------
    Number of neighbors to use for each cell. Returned as a Pandas DataFrame of size ``(n_obs x n_labeling_times)``.
    """
    n_neighbors = {}

    for labeling_time, mask in labeling_time_mask.items():
        if issparse(X):
            rep_X = csr_matrix(np.ones([X.shape[0], 1])) * X[mask, :].T
            if not sparse_op:
                rep_X = rep_X.toarray()
        else:
            rep_X = np.ones([X.size, 1]) * X[mask].T

        rows = np.tile(np.arange(X.shape[0]).reshape(-1, 1), n_obs[labeling_time])
        cols = dist_argsort[labeling_time]
        rep_X = rep_X[rows, cols]

        if sparse_op:
            n_neighbors_to_use = np.cumsum(rep_X.A > 0, axis=1)
        else:
            n_neighbors_to_use = np.cumsum(rep_X > 0, axis=1)

        # Check that number of non-trivial neighbors is large enough and reduce number of required neighbors
        # otherwise
        mask = n_neighbors_to_use.max(axis=1) < n_nontrivial_counts
        if mask.any():
            n_nontrivial_counts[mask] = n_neighbors_to_use[mask].max(axis=1)

        # First time `n_nontrivial_counts` neighbors with non-trivial counts
        n_neighbors_to_use = (
            n_neighbors_to_use == n_nontrivial_counts.reshape(-1, 1)
        ).argmax(axis=1)
        n_neighbors[labeling_time] = n_neighbors_to_use

    return pd.DataFrame(n_neighbors)


def get_labeling_times(adata, time_key) -> List:
    """Get labeling times in dataset.

    See :cite:p:`Weiler2023`.

    Arguments:
    ---------
    adata
        AnnData object.
    time_key
        Column of ``adata.obs`` containing labeling time information.

    Returns
    -------
    List of labeling times.
    """
    return adata.obs[time_key].unique()


def get_labeling_time_mask(
    adata: AnnData, time_key: str, labeling_times: List[float]
) -> Dict[float, np.ndarray]:
    """Get number of neighbors required to include ``n_nontrivial_counts`` counts per labeling time.

    See :cite:p:`Weiler2023`.

    Arguments:
    ---------
    adata
        AnnData object.
    time_key
        Column in `adata.obs` containing labeling time of each observation.
    labeling_times
        List of labeling times.

    Returns
    -------
    Dictionary with labeling times as keys and masks to subset to relevant observations as values.
    """
    # Mask for observation of each labeling time
    return {
        labeling_time: (adata.obs[time_key] == labeling_time).values
        for labeling_time in labeling_times
    }


def get_obs_dist_argsort(
    adata: AnnData, labeling_time_mask: Dict[float, np.ndarray]
) -> Dict[float, np.ndarray]:
    """Calculate argsorted pairwise distances per labeling_time_point.

    See :cite:p:`Weiler2023`.

    Arguments:
    ---------
    adata
        AnnData object.
    labeling_time_mask
        Dictionary assigning each labeling time a mask to subset observations.

    Returns
    -------
    Dictionary argsorted pairwise distances per labeling_time_point.
    """
    n_pcs = adata.uns["neighbors"]["params"]["n_pcs"]

    # Distance of each observation to all other observations
    obs_dist = pairwise_distances(adata.obsm["X_pca"][:, :n_pcs])

    # Distances of each observation to observations of one labeling time
    obs_dist = {
        labeling_time: obs_dist[:, obs_mask]
        for labeling_time, obs_mask in labeling_time_mask.items()
    }

    return {
        labeling_time: obs_dist[labeling_time].argsort(axis=1)
        for labeling_time in labeling_time_mask.keys()
    }


def get_n_neighbors(
    adata,
    labeling_time_mask: Dict[float, np.ndarray],
    obs_dist_argsort: Dict[float, np.ndarray],
    n_nontrivial_counts: int,
    use_rep="X",
    sparse_op: bool = False,
    n_jobs: Optional[int] = None,
) -> Dict[str, pd.DataFrame]:
    """Get number of neighbors required to include ``n_nontrivial_counts`` counts per labeling time.

    See :cite:p:`Weiler2023`.

    Arguments:
    ---------
    adata
        AnnData object.
    labeling_time_mask
        Dictionary with labeling times as keys and masks to subset to relevant observations as values.
    obs_dist_argsort
        Dictionary with argsorted pairwise distances per labeling_time_point.
    n_nontrivial_counts
        Number of non-trivial counts to consider for each labeling time point.
    use_rep
        Representation to use for identifying number of neighbors.
    sparse_op
        Boolean flag to run operations on sparse or dense matrix.

    Returns
    -------
    Number of neighbors to use for each gene and cell. Returned as a dictionary with variable names as keys
    and Pandas DataFrame of size ``(n_obs x n_labeling_times)`` as values.
    """
    if use_rep == "X":
        X = adata.X.copy()
    else:
        X = adata.layers[use_rep].copy()

    n_nontrivial_counts = np.array([n_nontrivial_counts] * adata.n_obs)

    # Number of observations per labeling time point
    n_obs_per_labeling_time = {
        labeling_time: argsorted_dists.shape[1]
        for labeling_time, argsorted_dists in obs_dist_argsort.items()
    }

    res = Parallel(n_jobs)(
        delayed(_get_n_neighbors)(
            X[:, var_id].copy(),
            labeling_time_mask=labeling_time_mask,
            dist_argsort=obs_dist_argsort,
            n_obs=n_obs_per_labeling_time,
            sparse_op=sparse_op,
            n_nontrivial_counts=n_nontrivial_counts.copy(),
        )
        for var_id in tqdm(range(adata.n_vars))
    )

    return dict(zip(adata.var_names, res))


def get_counts(gex, labeling_time_mask, obs_dist_argsort, obs_id, neighbors):
    """Return gex counts used for fitting parameters.

    Arguments:
    ---------
    gex
        GEX vector of a given gene.
    time_masks
        For each labeling time point, the masks to subset to relevant observations.
    obs_dist_argsort
        Argsorted distances for one observation.
    obs_id
        ID of observations for which to collect data.
    neighbors
        Number of neighbors to consider per labeling time point

    Returns
    -------
    Counts used for fitting parameters of cell ``obs_id``.
    """
    # IDs of observations in each labeling time to consider
    obs_ids = {
        idx: obs_dist_argsort[idx][obs_id, : (val + 1)]
        for idx, val in neighbors.iteritems()
    }
    # Stacked counts to consider
    counts = np.hstack(
        [
            gex[mask][obs_ids[labeling_time]]
            for labeling_time, mask in labeling_time_mask.items()
        ]
    )

    return counts


def _get_parameters(
    measured_labeled,
    labeling_times,
    experiment,
    constrained: bool = True,
    method: str = "SLSQP",
    x0: Optional[List] = None,
    **kwargs,
):
    """Estimates parameters of splicing kinetics from metabolic labeling data.

    Arguments:
    ---------
    measured_labeled
        Measured labeled RNA.
    labeling_times
        Labeling time of each observation.
    experiment
        Type of experiment (``"chase"``, ``"pulse"``, or ``both``) that observation was generated from.
    constrained
        Boolean flag to use constrained optimization or not.
    method
        Optimization method used to estimate parameters. Any valid option for ``scipy.optimize.minimize`` can be passed.
    x0
        The initial parameter estimate to use if provided.

    Returns
    -------
    Estimated transcription rate (alpha), degradation rate (gamma), initial RNA abundance (r0; if observations from chase experiment are
    included) and boolean flag if optimization terminated successfully.
    """
    if {"chase", "both"}.intersection(experiment):
        r0 = []

    if constrained:
        bounds = _get_bounds(experiment=experiment)

    if {"chase", "both"}.intersection(experiment):
        x0 = [10, 5, 17] if x0 is None else x0
    else:
        x0 = [10, 5] if x0 is None else x0

    res = minimize(
        get_mse,
        x0=x0,
        args=(measured_labeled, labeling_times, experiment),
        method=method,
        bounds=bounds,
        **kwargs,
    )

    alpha = res.x[0]
    gamma = res.x[1]
    if {"chase", "both"}.intersection(experiment):
        r0 = res.x[2]
    success = res.success

    return alpha, gamma, r0, success, res


def get_parameters(
    adata: AnnData,
    use_rep: Optional[str],
    time_key: str,
    experiment_key: str,
    n_neighbors,
    x0: List[float],
    n_jobs: Optional[int] = None,
):
    """Estimates parameters of splicing kinetics from metabolic labeling data.

    See :cite:p:`Weiler2023`.

    Arguments:
    ---------
    adata
        AnnData object containing data.
    use_rep
        Layer name containing labeled mRNA data.
    time_key
        Column in ``adata.obs`` with labeling time of observations.
    experiment_key
        Column in ``adata.obs`` with experimt type of observations.
    n_neighbors
        TODO
    x0
        The initial parameter estimate to use if provided.
    n_job
        Number of optimization problems to solve in parallel.

    Returns
    -------
        Estimated parameters alpha (transcription rate), gamma (degradation rate), r0 (initial GEX if
        data from chase experiment is included), and success (flag if optimization ran successfully).
    """
    if use_rep is None:
        X = adata.X.copy()
    else:
        X = adata.layers[use_rep].copy()

    time_to_experiment = dict(
        adata.obs[[time_key, experiment_key]].drop_duplicates().values
    )

    labeling_times = get_labeling_times(adata=adata, time_key=time_key)
    labeling_time_mask = get_labeling_time_mask(
        adata=adata, time_key=time_key, labeling_times=labeling_times
    )

    obs_dist_argsort = get_obs_dist_argsort(
        adata=adata, labeling_time_mask=labeling_time_mask
    )

    var_names = adata.var_names.tolist()
    obs_names = adata.obs_names.tolist()

    def _fit(obs_id: int):
        alpha = []
        gamma = []
        r0 = []
        success = []
        res = []
        for var_id, var_name in enumerate(var_names):
            neighbors = n_neighbors[var_name].iloc[obs_id, :]

            _counts = get_counts(
                gex=X[:, var_id],
                labeling_time_mask=labeling_time_mask,
                obs_dist_argsort=obs_dist_argsort,
                obs_id=obs_id,
                neighbors=neighbors,
            )

            _labeling_times = np.hstack(
                [[idx] * (val + 1) for idx, val in neighbors.iteritems()]
            )

            experiment = (
                pd.Series(_labeling_times)
                .replace(time_to_experiment)
                .str.lower()
                .values
            )

            params = _get_parameters(
                measured_labeled=_counts,
                labeling_times=_labeling_times,
                experiment=experiment,
                x0=x0,
            )
            alpha.append(params[0])
            gamma.append(params[1])
            r0.append(params[2])
            success.append(params[3])

            # TODO: Remove
            res.append(params[4])

        return alpha, gamma, r0, success, res

    res = Parallel(n_jobs)(
        delayed(_fit)(obs_id) for obs_id in tqdm(range(len(obs_names)))
    )

    alpha, gamma, r0, success, res = zip(*res)
    alpha = pd.DataFrame(np.array(alpha), index=obs_names, columns=var_names)
    gamma = pd.DataFrame(np.array(gamma), index=obs_names, columns=var_names)
    r0 = pd.DataFrame(np.array(r0), index=obs_names, columns=var_names)
    success = pd.DataFrame(np.array(success), index=obs_names, columns=var_names)

    return alpha, gamma, r0, success, res
