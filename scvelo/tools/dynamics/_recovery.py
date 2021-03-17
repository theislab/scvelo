import numpy as np
import pandas as pd

from scvelo.core import get_modality, set_modality
from scvelo.preprocessing.neighbors import get_connectivities


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
