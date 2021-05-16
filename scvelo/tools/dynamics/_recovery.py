import os
from typing import Dict, List, Optional, Tuple, Union

from typing_extensions import Literal

import numpy as np
import pandas as pd

from anndata import AnnData

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import get_modality, get_n_jobs, parallelize, set_modality
from scvelo.preprocessing.neighbors import get_connectivities


# TODO: Add docstrings
# TODO: Finish type hints
def recover_dynamics(
    adata: AnnData,
    modalities: List,
    dynamics: Union[str, List, Tuple],
    model_parameters: List,
    initial_parameter_fit: Dict,
    pretrain: List[List[str]],
    train: List[List[str]],
    max_iter: int = 10,
    var_names: str = "velocity_genes",
    fit_connected_states: bool = True,
    fit_scaling: bool = True,
    assignment_mode: Optional[
        Literal["full_projection", "partial_projection", "projection"]
    ] = "projection",
    align: bool = True,
    t_max: Optional[float] = None,
    inplace: bool = True,
    n_jobs=None,
    backend="loky",
    **kwargs,
) -> Optional[AnnData]:
    """Recover dynamical system for a given set of variables.

    Arguments
    ---------
    adata
        Annotated data containing data on which dynamical system is fitted.
    modalites
        Relevant modalities for fitting dynamical system. For example, in case of RNA
        velocity `modalities=["unspliced", "spliced"]`.
    dynamics
        Dynamics to recover. If not specified by a string, the first entry gives the
        dynamical system, the second entry the class to recover its model parameters.
    model_parameters
        Name of model parameters.
    initial_paremter_fit
            Specification of initial parameter fit through a grid. Each fit is specified
            by the parameter names `parameter_names` to fit, the width of the grid
            `sight` and the numbe of points `num` in it. For example, consider `alpha`
            being estimated initially as `alpha=1`. Using
            `{"parameter_names": ["alpha"], "sight": 0.5, "num": 5}`, `alpha` is updated
            by the value of `array([0.5 , 0.75, 1.  , 1.25, 1.5])` resulting in the
            minimal loss.
    pretrain
            Parameter pairing to fit with time assignment.
    train
        Parameter pairings to fit with fixed time assignment.
    max_iter
        Maximal iterations in the optimization algorithm.
    var_names
        Name of variables for which dynamical system to fit.
    fit_connected_states
        Boolean flag to restrict fit to neighbors given by connectivies.
    fit_scaling
        Boolean flag to fit scaling of counts or not.
    assignment_mode
        Mode used to assign time points.
    align
        Boolean flag to align inferred time points.
    t_max
        Time point of last observations (in a temporal sense). Will be estimated if not
        specified.
    inplace
        Whether or not to manipulate `adata` inplace.
    n_jobs
        Number of cores used during parallelized parameter inference.
    backend
        Backend used by `joblib.Parallel`.
    kwargs
        Additional keyword arguments for recovering dynamics.

    Returns
    -------
    Optional[AnnData]
        Returns updated annotated data matrix if `inplace=False`, `None` otherwise.
    """

    if not isinstance(dynamics, str):
        DynamicsRecovery = dynamics[1]
        dynamics = dynamics[0]

    if not inplace:
        adata = adata.copy()

    if isinstance(var_names, str) and var_names not in adata.var_names:
        if var_names in adata.var.columns:
            var_names = adata.var_names[adata.var[var_names].values]
        elif var_names == "all":
            var_names = adata.var_names

    connectivities = get_connectivities(adata) if fit_connected_states else None

    fitted_time = np.empty(adata.shape)
    fitted_time.fill(np.nan)

    fitted_tau = np.empty(adata.shape)
    fitted_tau.fill(np.nan)

    fitted_tau_ = np.empty(adata.shape)
    fitted_tau_.fill(np.nan)

    fitted_parameters = pd.DataFrame(
        np.zeros((adata.n_vars, len(model_parameters) + 2 + len(modalities) + 1))
        * np.nan,
        index=adata.var_names,
        columns=model_parameters
        + ["t_", "scaling"]
        + [f"steady_{modality}" for modality in modalities]
        + ["likelihood"],
    )
    fitted_parameters = fitted_parameters.add_prefix("fit_")

    loss = pd.DataFrame()
    likelihood = pd.DataFrame()

    var_idx = []

    n_jobs = get_n_jobs(n_jobs=n_jobs)
    logg.info(f"recovering dynamics (using {n_jobs}/{os.cpu_count()} cores)", r=True)

    res = parallelize(
        _fit_recovery,
        var_names,
        n_jobs,
        unit="gene",
        as_array=False,
        backend=backend,
        show_progress_bar=len(var_names) > 9,
    )(
        adata=adata,
        modalities=modalities,
        DynamicsRecovery=DynamicsRecovery,
        dynamics=dynamics,
        model_parameters=model_parameters,
        connectivities=connectivities,
        max_iter=max_iter,
        fit_scaling=fit_scaling,
        assignment_mode=assignment_mode,
        initial_parameter_fit=initial_parameter_fit,
        pretrain=pretrain,
        train=train,
        **kwargs,
    )
    var_names, dms = map(_flatten, zip(*res))

    for var_name, dm in zip(var_names, dms):
        loss = loss.append(pd.DataFrame(dm.loss_, columns=[var_name]).T)
        likelihood = likelihood.append(
            pd.DataFrame(dm.likelihood_, columns=[var_name]).T
        )

        var_id = adata.var_names.get_loc(var_name)
        var_idx.append(var_id)

        fitted_time[:, var_id] = dm.t
        fitted_tau[:, var_id] = dm.tau
        fitted_tau_[:, var_id] = dm.tau_

        dm.steady_state *= dm.scaling

        fitted_parameters.loc[var_name, :] = np.array(
            np.concatenate(
                (dm.params.iloc[-1, :], dm.steady_state, [dm.likelihood_[-1]])
            )
        )

    adata.layers["fit_t"] = (
        fitted_time if connectivities is None else connectivities.dot(fitted_time)
    )
    adata.layers["fit_tau"] = fitted_tau
    adata.layers["fit_tau_"] = fitted_tau_
    adata.var = adata.var.merge(fitted_parameters, left_index=True, right_index=True)

    if "loss" not in adata.varm:
        adata.varm["loss"] = pd.DataFrame(index=adata.var_names)
    adata.varm["loss"] = (
        adata.varm["loss"]
        .merge(loss, left_index=True, right_index=True, how="outer")
        .add_prefix("update_")
    )
    if "likelihood" not in adata.varm:
        adata.varm["likelihood"] = pd.DataFrame(index=adata.var_names)
    adata.varm["likelihood"] = (
        adata.varm["likelihood"]
        .merge(likelihood, left_index=True, right_index=True, how="outer")
        .add_prefix("update_")
    )

    if align:
        align_dynamics(
            adata, model_parameters=model_parameters, t_max=t_max, var_idx=var_idx
        )

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")

    if not inplace:
        return adata


# TODO: Add docstrings
# TODO: Finish type hints
def align_dynamics(
    adata: AnnData,
    model_parameters: List,
    t_max: Optional[float] = None,
    var_idx: Optional[List] = None,
    mode: Literal["align_total_time"] = "align_total_time",
):
    """Align inferred time points.

    Arguments
    ---------
    adata
        Annotated data matrix.
    model_parameters
        Parameter names of dynamical system.
    t_max
        Time point of last observations (in a temporal sense). Will be estimated if not
        specified.
    var_idx
        Indices of variables for which dynamical systems have been fitted.
    mode
        Type of alignment.
    """

    if t_max is None:
        t_max = 20

    if var_idx is None:
        var_idx = np.arange(0, adata.n_vars)

    m = np.ones(adata.n_vars)

    t_ = adata.var["fit_t_"].iloc[var_idx].values
    fitted_time = adata[:, var_idx].layers["fit_t"]

    if mode == "align_total_time" and t_max is not False:
        T_max = np.max(fitted_time * (fitted_time < t_), axis=0)
        T_max += np.max((fitted_time - t_) * (fitted_time > t_), axis=0)

        denom = (
            1 - np.sum((fitted_time == t_) | (fitted_time == 0), axis=0) / adata.n_obs
        )
        denom += denom == 0

        T_max = T_max / denom
        T_max += T_max == 0

        m[var_idx] = t_max / T_max

    m = m[var_idx]

    for parameter in model_parameters:
        adata.var[f"fit_{parameter}"].iloc[var_idx] = (
            adata.var[f"fit_{parameter}"].iloc[var_idx] / m
        )
    adata.var["fit_t_"].iloc[var_idx] = adata.var["fit_t_"].iloc[var_idx] * m

    for layer_name in ["fit_t", "fit_tau", "fit_tau_"]:
        set_modality(adata, adata[:, var_idx].layers[layer_name] * m, layer_name)


def _fit_recovery(
    var_names,
    adata,
    modalities,
    DynamicsRecovery,
    dynamics,
    model_parameters,
    connectivities,
    max_iter,
    fit_scaling,
    assignment_mode,
    initial_parameter_fit,
    pretrain,
    train,
    queue,
    **kwargs,
):

    _var_names, dms = [], []
    for var_name in var_names:
        dm = DynamicsRecovery(
            dynamics=dynamics,
            model_parameters=model_parameters,
            connectivities=connectivities,
            max_iter=max_iter,
            fit_scaling=fit_scaling,
        )
        dm.fit(
            counts=np.hstack(
                [get_modality(adata[:, var_name], modality) for modality in modalities]
            ),
            initial_parameter_fit=initial_parameter_fit,
            pretrain=pretrain,
            train=train,
            assignment_mode=assignment_mode,
        )

        dms.append(dm)
        _var_names.append(var_name)

        if queue is not None:
            queue.put(1)

    if queue is not None:
        queue.put(None)

    return _var_names, dms


def _flatten(iterable):
    return [i for it in iterable for i in it]
