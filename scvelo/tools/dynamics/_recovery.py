import os
import warnings
from typing import Dict, List, Optional, Tuple, Union

from typing_extensions import Literal

import numpy as np
import pandas as pd
from numpy import ndarray

from anndata import AnnData

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import (
    clipped_log,
    get_modality,
    get_n_jobs,
    invert,
    LinearRegression,
    parallelize,
    set_modality,
    SplicingDynamics,
)
from scvelo.preprocessing.neighbors import get_connectivities
from scvelo.tools.dynamical_model_utils import convolve
from ._dynamics_base import DynamicsRecoveryBase


# TODO: Remove argument `model_parameters` and infer them from the class `dynamics`.
class SplicingDynamicsRecovery(DynamicsRecoveryBase):
    def _initialize_parameters(self, subsetted_counts: ndarray):
        """Initialize model parameters.

        Arguments
        ---------
        obs_subset_counts
            Filter mask to subset observations.
        """

        # TODO: Remove hard coded percentile
        _obs_subset = subsetted_counts >= np.percentile(subsetted_counts, 98, axis=0)

        u_inf = subsetted_counts[:, 0][_obs_subset[:, 0] | _obs_subset[:, 1]].mean()
        s_inf = subsetted_counts[:, 1][_obs_subset[:, 1]].mean()

        obs_subset_g = _obs_subset[:, 1] | self.steady_state_prior[self.obs_subset_]

        self.gamma = (
            LinearRegression()
            .fit(
                convolve(subsetted_counts[:, 1], obs_subset_g),
                convolve(subsetted_counts[:, 0], obs_subset_g),
            )
            .coef_
            + 1e-6
        )

        if self.gamma < 0.05 / self.scaling[0]:
            self.gamma *= 1.2
        elif self.gamma > 1.5 / self.scaling[0]:
            self.gamma /= 1.2

        if self.pval_steady < 1e-3:
            u_inf = np.mean([u_inf, self.steady_state[0]])
            self.alpha = self.gamma * s_inf
            self.beta = self.alpha / u_inf
            self.initial_state_ = np.array([u_inf, s_inf])
        else:
            self.initial_state_ = np.array([u_inf, s_inf])
            self.alpha = u_inf
            self.beta = 1

        self.alpha_ = 0

    def _set_steady_state_ratio(self, **model_params):
        """Set ratio between steady states of dynamical system.

        Arguments
        ---------
        model_params
            Parameters of dynamical system and their values.
        """

        self.steady_state_ratio = model_params["gamma"] / model_params["beta"]

    # TODO: Finish docstrings
    # TODO: Find better name
    def _check_projection(self, **model_parameters) -> bool:
        """

        Arguments
        ---------
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        bool
        """

        if model_parameters["beta"] < model_parameters["gamma"]:
            return True
        else:
            return False

    def get_approx_time_assignment(
        self,
        state: ndarray,
        alpha: Union[float, ndarray],
        beta: Union[float, ndarray],
        gamma: Union[float, ndarray],
        initial_state: Union[ndarray, List] = [0, 0],
        full_projection: bool = False,
    ) -> ndarray:
        """Get approximate time assignment.

        Arguments
        ---------
        state
            Observed values considered.
        initial_state
            Initial state of dynamical system.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        ndarray
            Approximate time assignment.
        """

        u = state[:, 0]
        s = state[:, 1]

        if full_projection:
            u = np.min(state[s > 0, 0])

        u0 = initial_state[0]
        s0 = initial_state[1]

        if (gamma >= beta) or full_projection:
            uinf = alpha / beta
            tau = -1 / beta * clipped_log((u - uinf) / (u0 - uinf))
        elif not full_projection:
            beta_ = beta * invert(gamma - beta)
            xinf = alpha / gamma - beta_ * (alpha / beta)
            tau = (
                -1
                / gamma
                * clipped_log((s - beta_ * u - xinf) / (s0 - beta_ * u0 - xinf))
            )

        if tau.size == 1 and not np.isscalar(tau):
            return tau[0]
        else:
            return tau

    def _set_parameters_after_switch(self, model_parameters):
        """Set model parameters during repression phase.

        Arguments
        ---------
        model_parameters
            Dictionary of model parameters and their current values.

        Returns
        -------
        Dict
            Model parameters during repression.
        """

        model_parameters["alpha"] = 0

        return model_parameters

    def get_residuals(
        self,
        t: Optional[ndarray] = None,
        t_: Optional[float] = None,
        scaling: Optional[Union[float, ndarray]] = None,
        initial_state_: Optional[ndarray] = None,
        refit_time: Optional[bool] = None,
        subsetted: bool = True,
        obs_subset_cluster: Optional[List] = None,
        return_model_kwargs: bool = False,
        **model_parameters,
    ) -> ndarray:
        """Get residuals of estimated trajectory and measurements.

        Arguments
        ---------
        t
            Time assigned to observations.
        t_
            Time when system switches states.
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        initial_state_
            State of system at switching point.
        refit_time
            Boolean flag to refit time assignment or not.
        subsetted
            Boolean flag to subset observations or not.
        obs_subset_cluster
            TODO: Add description.
        return_model_kwargs
            Boolean flag to return model parameters in addition to residuals.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        ndarray
            Residuals.
        """

        model_parameters = self.get_model_parameters(**model_parameters)
        scaling, t_ = self.get_vars(scaling, t_, initial_state_, **model_parameters)

        subsetted_counts = self.get_counts(
            scaling, subsetted=subsetted, obs_subset_cluster=obs_subset_cluster
        )

        t, tau, o = self.get_time_assignment(
            scaling,
            t_,
            initial_state_,
            t,
            refit_time,
            subsetted=subsetted,
            obs_subset_cluster=obs_subset_cluster,
            **model_parameters,
        )

        tau, alpha, initial_state = self._vectorize(t, t_, **model_parameters)
        model_parameters.update({"alpha": alpha})

        sol = self.dynamics(
            initial_state=initial_state, **model_parameters
        ).get_solution(tau)

        if return_model_kwargs:
            return model_parameters, (sol - subsetted_counts) / self.std_ * scaling
        else:
            return (sol - subsetted_counts) / self.std_ * scaling

    def _get_regularization(
        self, subsetted_counts: ndarray, **model_parameters
    ) -> ndarray:
        """Calculate regularization term.

        Arguments
        ---------
        subsetted_counts
            Count matrix of subsetted observations.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        ndarray
            Regularization term.
        """

        beta = model_parameters["beta"]
        gamma = model_parameters["gamma"]
        return (
            (gamma / beta - self.steady_state_ratio)
            * subsetted_counts[:, 1]
            / self.std_[1]
        )

    # TODO: Add as method to a base class for dynamical models
    def _vectorize(
        self,
        t: ndarray,
        t_: float,
        alpha_: float = 0,
        initial_state: Union[ndarray, List] = [0, 0],
        sorted: bool = False,
        **model_parameters,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        """Vectorize parameters and initial states.

        Arguments
        ---------
        t
            Time assigned to observations.
        t_
            Time when system switches states.
        initial_state
            Initial state of system.
        sorted
            Boolean flag to sort vectorized variables according to time assigment.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        Tuple[ndarray, ndarray, ndarray]
            Time assignments w.r.t. state of system, vectorized value of transcription
            and initial state.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            o = np.array(t < t_, dtype=int)
        tau = t * o + (t - t_) * (1 - o)

        state_after_switch = self.dynamics(
            initial_state=initial_state, **model_parameters
        ).get_solution(t_)

        # vectorize u0, s0 and alpha
        initial_state = (
            np.array(initial_state) * o[:, None] + state_after_switch * (1 - o)[:, None]
        )
        alpha = model_parameters["alpha"] * o + alpha_ * (1 - o)

        if sorted:
            idx = np.argsort(t)
            tau, alpha = tau[idx], alpha[idx]
            initial_state = initial_state[idx, :]
        return tau, alpha, initial_state


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
    elif dynamics.lower() == "splicing":
        dynamics = SplicingDynamics
        DynamicsRecovery = SplicingDynamicsRecovery

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
