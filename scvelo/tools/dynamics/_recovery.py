import warnings

import numpy as np
import pandas as pd

from scvelo.core import (
    clipped_log,
    get_modality,
    invert,
    LinearRegression,
    set_modality,
    SplicingDynamics,
)
from scvelo.preprocessing.neighbors import get_connectivities
from scvelo.tools.dynamical_model_utils import convolve
from ._recovery_base import DynamicsRecoveryBase


# TODO: Add docstrings
# TODO: Finish type hints
# TODO: Remove argument `model_parameters` and infer them from the class `dynamics`.
class SplicingDynamicsRecovery(DynamicsRecoveryBase):
    def _initialize_parameters(self, weighted_counts):
        # TODO: Remove hard coded percentile
        _weights = weighted_counts >= np.percentile(weighted_counts, 98, axis=0)

        u_inf = weighted_counts[:, 0][_weights[:, 0] | _weights[:, 1]].mean()
        s_inf = weighted_counts[:, 1][_weights[:, 1]].mean()

        weights_g = _weights[:, 1] | self.steady_state_prior[self.weights_]

        self.gamma = (
            LinearRegression()
            .fit(
                convolve(weighted_counts[:, 1], weights_g),
                convolve(weighted_counts[:, 0], weights_g),
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
        self.steady_state_ratio = model_params["gamma"] / model_params["beta"]

    # TODO: Find better name
    def _check_projection(self, **model_parameters):
        if model_parameters["beta"] < model_parameters["gamma"]:
            return True
        else:
            return False

    def tau_inv(
        self,
        state,
        alpha,
        beta,
        gamma,
        initial_state=[0, 0],
        full_projection=False,
    ):
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
        model_parameters["alpha"] = 0

        return model_parameters

    def get_residuals(
        self,
        t=None,
        t_=None,
        scaling=None,
        initial_state_=None,
        refit_time=None,
        weighted=True,
        weights_cluster=None,
        return_model_kwargs=False,
        **model_parameters,
    ):
        model_parameters = self.get_model_parameters(**model_parameters)
        scaling, t_ = self.get_vars(scaling, t_, initial_state_, **model_parameters)

        weighted_counts = self.get_counts(
            scaling, weighted=weighted, weights_cluster=weights_cluster
        )

        t, tau, o = self.get_time_assignment(
            scaling,
            t_,
            initial_state_,
            t,
            refit_time,
            weighted=weighted,
            weights_cluster=weights_cluster,
            **model_parameters,
        )

        tau, alpha, initial_state = self._vectorize(t, t_, **model_parameters)
        model_parameters.update({"alpha": alpha})

        sol = self.dynamics(
            initial_state=initial_state, **model_parameters
        ).get_solution(tau)

        if return_model_kwargs:
            return model_parameters, (sol - weighted_counts) / self.std_ * scaling
        else:
            return (sol - weighted_counts) / self.std_ * scaling

    def _get_regularization(self, weighted_counts, **model_parameters):
        beta = model_parameters["beta"]
        gamma = model_parameters["gamma"]
        return (
            (gamma / beta - self.steady_state_ratio)
            * weighted_counts[:, 1]
            / self.std_[1]
        )

    # TODO: Add as method to a base class for dynamical models
    def _vectorize(
        self, t, t_, alpha_=0, initial_state=[0, 0], sorted=False, **model_parameters
    ):
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
    adata,
    modalities,
    dynamics,
    model_parameters,
    initial_parameter_fit,
    pretrain,
    train,
    max_iter: int = 10,
    var_names="velocity_genes",
    fit_connected_states: bool = True,
    fit_scaling: bool = True,
    assignment_mode="projection",
    align=True,
    t_max=None,
    inplace=True,
):
    """"""

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

    var_idx = []

    dm = DynamicsRecovery(
        dynamics=dynamics,
        model_parameters=model_parameters,
        connectivities=connectivities,
        max_iter=max_iter,
        fit_scaling=fit_scaling,
    )

    # TODO: Check if better solution exists than setting dm.assignment_mode to None
    for var in var_names:
        dm.assignment_mode = None
        dm.initialize(
            np.hstack(
                [get_modality(adata[:, var], modality) for modality in modalities]
            ),
            initial_parameter_fit=initial_parameter_fit,
        ).fit(pretrain=pretrain, train=train, assignment_mode=assignment_mode)
        loss = loss.append(pd.DataFrame(dm.loss_, columns=[var]).T)

        var_id = adata.var_names.get_loc(var)
        var_idx.append(var_id)

        fitted_time[:, var_id] = dm.t
        fitted_tau[:, var_id] = dm.tau
        fitted_tau_[:, var_id] = dm.tau_

        dm.steady_state *= dm.scaling

        fitted_parameters.loc[var, :] = np.array(
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

    if align:
        align_dynamics(
            adata, model_parameters=model_parameters, t_max=t_max, var_idx=var_idx
        )

    if not inplace:
        return adata


# TODO: Add docstrings
# TODO: Finish type hints
def align_dynamics(
    adata, model_parameters, t_max=None, var_idx=None, mode="align_total_time"
):
    """"""

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
