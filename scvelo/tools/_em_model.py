import os
from dataclasses import asdict, dataclass, field, fields
from typing import List, Optional, Sequence, Union

import numpy as np

from anndata import AnnData

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import get_n_jobs, parallelize
from scvelo.preprocessing.moments import get_connectivities
from ._core import BaseInference
from ._em_model_core import DynamicsRecovery
from ._steady_state_model import SteadyStateModel
from .utils import make_unique_list


@dataclass
class EMParams:
    """EM parameters."""

    r2: np.ndarray = field(metadata={"is_matrix": False})
    alpha: np.ndarray = field(metadata={"is_matrix": False})
    beta: np.ndarray = field(metadata={"is_matrix": False})
    gamma: np.ndarray = field(metadata={"is_matrix": False})
    t_: np.ndarray = field(metadata={"is_matrix": False})
    scaling: np.ndarray = field(metadata={"is_matrix": False})
    std_u: np.ndarray = field(metadata={"is_matrix": False})
    std_s: np.ndarray = field(metadata={"is_matrix": False})
    likelihood: np.ndarray = field(metadata={"is_matrix": False})
    u0: np.ndarray = field(metadata={"is_matrix": False})
    s0: np.ndarray = field(metadata={"is_matrix": False})
    pval_steady: np.ndarray = field(metadata={"is_matrix": False})
    steady_u: np.ndarray = field(metadata={"is_matrix": False})
    steady_s: np.ndarray = field(metadata={"is_matrix": False})
    variance: np.ndarray = field(metadata={"is_matrix": False})
    alignment_scaling: np.ndarray = field(metadata={"is_matrix": False})
    T: np.ndarray = field(metadata={"is_matrix": True})
    Tau: np.ndarray = field(metadata={"is_matrix": True})
    Tau_: np.ndarray = field(metadata={"is_matrix": True})

    @classmethod
    def from_adata(cls, adata: AnnData, key: str = "fit"):
        parameter_dict = {}
        for parameter in fields(cls):
            para_name = parameter.name
            if parameter.metadata["is_matrix"]:
                if f"{key}_{para_name.lower()}" in adata.layers.keys():
                    parameter_dict[para_name] = adata.layers[
                        f"{key}_{para_name.lower()}"
                    ]
                else:
                    _vals = np.empty(adata.shape)
                    _vals.fill(np.nan)
                    parameter_dict[para_name] = _vals
            else:
                if f"{key}_{para_name.lower()}" in adata.var.keys():
                    parameter_dict[para_name] = adata.var[
                        f"{key}_{para_name.lower()}"
                    ].values
                else:
                    _vals = np.empty(adata.n_vars)
                    _vals.fill(np.nan)
                    parameter_dict[para_name] = _vals
        return cls(**parameter_dict)

    def export_to_adata(self, adata: AnnData, key: str = "fit"):
        for parameter in fields(self):
            para_name = parameter.name
            value = getattr(self, para_name)
            # The parameter is only written if not all entries are nan.
            if not np.all(np.isnan(value)):
                if parameter.metadata["is_matrix"]:
                    adata.layers[f"{key}_{para_name.lower()}"] = value
                else:
                    adata.var[f"{key}_{para_name.lower()}"] = value
        return adata


# TODO: Implement abstract methods
class ExpectationMaximizationModel(BaseInference):
    """EM 'Dynamical' model for velocity estimation.

    Parameters
    ----------
    adata
        Annotated data matrix.
    spliced_layer
        Key for spliced layer.
    unspliced_layer
        Key for unspliced layer.
    var_names
        Names of variables/genes to use for the fitting. If `var_names='velocity_genes'`
        but there is no column `'velocity_genes'` in `adata.var`, velocity genes are
        estimated using the steady state model.
    n_top_genes
        Number of top velocity genes to use for the dynamical model.
    max_iter
        Maximum number of iterations.
    projection_assignment
        Whether to use projection assignment (default), or inverse approximation.
    t_max
        Maximum time point to use for the fitting.
    fit_time
        Whether to fit time or keep initially given time fixed.
    fit_scaling
        Whether to fit scaling between unspliced and spliced.
    fit_steady_states
        Whether to explicitly model and fit steady states (next to induction/repression).
    fit_connected_states
        Restricts fitting to neighbors given by connectivities.
    fit_basal_transcription
        Enables model to incorporate basal transcriptions.
    steady_state_prior
        Mask for indices used for steady state regression.
    n_jobs
        Number of jobs for parallelization.
    backend
        Backend used for multiprocessing. See :class:`joblib.Parallel` for valid
        options.
    """

    def __init__(
        self,
        adata: AnnData,
        spliced_layer: str = "Ms",
        unspliced_layer: str = "Mu",
        var_names_key: str = "velocity_genes",
        fit_key: str = "fit",
        n_top_genes: Optional[int] = None,
        max_iter: int = 10,
        projection_assignment: bool = True,
        t_max: int = 20,
        fit_time: bool = True,
        fit_scaling: bool = True,
        fit_steady_states: bool = True,
        fit_connected_states: bool = None,
        fit_basal_transcription: bool = None,
        steady_state_prior: Optional[Sequence[bool]] = None,
        n_jobs: Optional[int] = None,
        backend: str = "loky",
    ):
        super().__init__(
            adata,
        )
        self._spliced_layer = spliced_layer
        self._unspliced_layer = unspliced_layer
        self._var_names_key = var_names_key
        self._fit_key = fit_key
        self._n_top_genes = n_top_genes
        self._max_iter = max_iter
        self._projection_assignment = projection_assignment
        self._t_max = t_max
        self._fit_time = fit_time
        self._fit_scaling = fit_scaling
        self._fit_steady_states = fit_steady_states
        self._fit_connected_states = fit_connected_states
        self._fit_basal_transcription = fit_basal_transcription
        self._steady_state_prior = steady_state_prior
        self._n_jobs = get_n_jobs(n_jobs=n_jobs)
        self._backend = backend

        if len(set(self._adata.var_names)) != len(self._adata.var_names):
            logg.warn("Duplicate var_names found. Making them unique.")
            self._adata.var_names_make_unique()
        self._state_dict = EMParams.from_adata(adata)
        self._use_raw = False

    def _prepare_genes(self):
        """Initialize genes to use for the fitting."""
        var_names = self._var_names_key
        adata = self._adata
        if isinstance(var_names, str) and var_names not in adata.var_names:
            if var_names in adata.var.keys():
                var_names = adata.var_names[adata.var[var_names].values]
            elif var_names == "all":
                var_names = adata.var_names
            elif "_genes" in var_names:
                velo = SteadyStateModel(
                    adata,
                    spliced_layer=self._spliced_layer,
                    unspliced_layer=self._unspliced_layer,
                )
                velo.fit()
                var_names = adata.var_names[velo.state_dict()["velocity_genes"]]
                self._state_dict.r2 = velo.state_dict()["r2"]
            else:
                raise ValueError("Variable name not found in var keys.")
        if not isinstance(var_names, str):
            var_names = list(np.ravel(var_names))

        var_names = make_unique_list(var_names, allow_array=True)
        var_names = np.array([name for name in var_names if name in adata.var_names])
        if len(var_names) == 0:
            raise ValueError("Variable name not found in var keys.")
        if self._n_top_genes is not None and len(var_names) > self._n_top_genes:
            X = adata.layers[self._spliced_layer][:, var_names]
            var_names = var_names[np.argsort(np.sum(X, 0))[::-1][: self._n_top_genes]]
        self._var_names = var_names

    def state_dict(self):
        """Return the state of the model."""
        return asdict(self._state_dict)

    def export_results_adata(self, copy: bool = True, add_key: str = "fit"):
        """Export the results to the AnnData object and return it."""
        adata = self._adata.copy() if copy else self._adata
        self._state_dict.export_to_adata(adata, add_key)
        adata.uns["recover_dynamics"] = {
            "fit_connected_states": self._fit_connected_states,
            "fit_basal_transcription": self._fit_basal_transcription,
            "use_raw": self._use_raw,
        }
        # loss is only written after the execution of fit()
        if hasattr(self, "_loss"):
            adata.varm["loss"] = self._loss
        return adata

    # TODO: Remove `use_raw` argument
    # TODO: Remove `return_model` argument
    def fit(
        self,
        copy: bool = True,
        use_raw: bool = False,
        return_model: Optional[bool] = None,
        load_pars: bool = False,
        steady_state_prior: Optional[List[bool]] = None,
        assignment_mode: str = "projection",
        show_progress_bar: bool = True,
        **kwargs,
    ):
        """Fit the model."""
        logg.info(
            f"recovering dynamics (using {self._n_jobs}/{os.cpu_count()} cores)", r=True
        )

        if (
            "Ms" not in self._adata.layers.keys()
            or "Mu" not in self._adata.layers.keys()
        ):
            use_raw = True
        # TODO: Refactor the definition of 'use_raw'; Move to init?
        self._use_raw = use_raw
        if self._fit_connected_states is None:
            self._fit_connected_states = not use_raw

        self._prepare_genes()
        if return_model is None:
            return_model = len(self._var_names) < 5

        sd = self._state_dict
        idx, L = [], []
        conn = get_connectivities(self._adata) if self._fit_connected_states else None

        res = parallelize(
            self._fit,
            self._var_names,
            self._n_jobs,
            unit="gene",
            as_array=False,
            backend=self._backend,
            show_progress_bar=show_progress_bar,
        )(
            use_raw=use_raw,
            load_pars=load_pars,
            steady_state_prior=steady_state_prior,
            conn=conn,
            assignment_mode=assignment_mode,
            **kwargs,
        )
        idx, dms = map(_flatten, zip(*res))

        for ix, dm in zip(idx, dms):
            sd.T[:, ix], sd.Tau[:, ix], sd.Tau_[:, ix] = dm.t, dm.tau, dm.tau_
            (
                sd.alpha[ix],
                sd.beta[ix],
                sd.gamma[ix],
                sd.t_[ix],
                sd.scaling[ix],
            ) = dm.pars[:, -1]
            sd.u0[ix], sd.s0[ix], sd.pval_steady[ix] = dm.u0, dm.s0, dm.pval_steady
            sd.steady_u[ix], sd.steady_s[ix] = dm.steady_u, dm.steady_s
            sd.beta[ix] /= sd.scaling[ix]
            sd.steady_u[ix] *= sd.scaling[ix]

            sd.std_u[ix], sd.std_s[ix] = dm.std_u, dm.std_s
            sd.likelihood[ix], sd.variance[ix] = dm.likelihood, dm.varx
            L.append(dm.loss)

        adata = self._adata

        if conn is not None:
            if f"{self._fit_key}_t" in adata.layers.keys():
                sd.T[:, idx] = conn.dot(sd.T[:, idx])
            else:
                sd.T = conn.dot(sd.T)

        # is False if only one invalid / irrecoverable gene was given in var_names
        if L:
            cur_len = adata.varm["loss"].shape[1] if "loss" in adata.varm.keys() else 2
            max_len = max(np.max([len(loss) for loss in L]), cur_len) if L else cur_len
            self._loss = np.empty((adata.n_vars, max_len))
            self._loss.fill(np.nan)

            if "loss" in adata.varm.keys():
                self._loss[:, :cur_len] = adata.varm["loss"]

            self._loss[idx] = np.vstack(
                [
                    np.concatenate([loss, np.ones(max_len - len(loss)) * np.nan])
                    for loss in L
                ]
            )

        # TODO: Fix s.t. `self._t_max` is only integer
        if self._t_max is not False:
            dm = self._align_dynamics(t_max=self._t_max, dm=dm, idx=idx)

        logg.info(
            "    finished", time=True, end=" " if settings.verbosity > 2 else "\n"
        )
        logg.hint(
            "added \n"
            f"    '{self._fit_key}_pars', "
            f"fitted parameters for splicing dynamics (adata.var)"
        )

        if return_model:
            logg.info("\noutputs model fit of gene:", dm.gene)

        return dm if return_model else adata if copy else None

    def _align_dynamics(
        self,
        t_max: Optional[Union[float, bool]] = None,
        dm: Optional[DynamicsRecovery] = None,
        idx: Optional[List[bool]] = None,
        mode: Optional[str] = None,
        remove_outliers: bool = False,
    ):
        """Align dynamics to a common set of parameters.

        Arguments:
        ---------
        t_max: `float`, `False` or `None` (default: `None`)
            Total range for time assignments.
        dm: :class:`~DynamicsRecovery`
            DynamicsRecovery object to perform alignment on.
        idx: list of `bool` or `None` (default: `None`)
            Mask for indices used for alignment.
        mode: `str` or None (default: `'align_total_time`)
            What to align. Takes the following arguments:
            common_splicing_rate, common_scaling, align_increments, align_total_time
        remove_outliers: `bool` or `None` (default: `None`)
            Whether to remove outliers.
        copy: `bool` (default: `False`)
            Return a copy instead of writing to `adata`.

        Returns
        -------
        `alpha`, `beta`, `gamma`, `t_`, `alignment_scaling`: `.var`
            Aligned parameters
        `fit_t`, `fit_tau`, `fit_tau_`: `.layer`
            Aligned time
        """
        sd = self._state_dict
        if idx is None:
            idx = ~np.isnan(np.sum(sd.T, axis=0))
        if np.all(np.isnan(sd.alignment_scaling)):
            sd.alignment_scaling.fill(1)
        if mode is None:
            mode = "align_total_time"

        m = np.ones(self._adata.n_vars)
        mz = sd.alignment_scaling
        mz_prev = np.array(mz)

        if dm is not None:  # newly fitted
            mz[idx] = 1

        if mode == "align_total_time" and t_max is not False:
            # transient 'on'
            T_max = np.max(sd.T[:, idx] * (sd.T[:, idx] < sd.t_[idx]), axis=0)
            # 'off'
            T_max += np.max(
                (sd.T[:, idx] - sd.t_[idx]) * (sd.T[:, idx] > sd.t_[idx]), axis=0
            )

            denom = 1 - np.sum(
                (sd.T[:, idx] == sd.t_[idx]) | (sd.T[:, idx] == 0), axis=0
            ) / len(sd.T)
            denom += denom == 0

            T_max = T_max / denom
            T_max += T_max == 0

            if t_max is None:
                t_max = 20
            m[idx] = t_max / T_max
            mz *= m
        else:
            m = 1 / mz
            mz = np.ones(self.adata.n_vars)

        if remove_outliers:
            mu, std = np.nanmean(mz), np.nanstd(mz)
            mz = np.clip(mz, mu - 3 * std, mu + 3 * std)
            m = mz / mz_prev

        if idx is None:
            sd.alpha, sd.beta, sd.gamma = sd.alpha / m, sd.beta / m, sd.gamma / m
            sd.T, sd.t_, sd.Tau, sd.Tau_ = sd.T * m, sd.t_ * m, sd.Tau * m, sd.Tau_ * m
        else:
            m_ = m[idx]
            sd.alpha[idx] = sd.alpha[idx] / m_
            sd.beta[idx] = sd.beta[idx] / m_
            sd.gamma[idx] = sd.gamma[idx] / m_
            sd.T[:, idx], sd.t_[idx] = sd.T[:, idx] * m_, sd.t_[idx] * m_
            sd.Tau[:, idx], sd.Tau_[:, idx] = sd.Tau[:, idx] * m_, sd.Tau_[:, idx] * m_

        mz[mz == 1] = np.nan

        if dm is not None and dm.recoverable:
            dm.m = m[idx]
            dm.alpha = dm.alpha / dm.m[-1]
            dm.beta = dm.beta / dm.m[-1]
            dm.gamma = dm.gamma / dm.m[-1]
            dm.pars[:3] = dm.pars[:3] / dm.m[-1]

            dm.t = dm.t * dm.m[-1]
            dm.tau = dm.tau * dm.m[-1]
            dm.t_ = dm.t_ * dm.m[-1]
            dm.pars[4] = dm.pars[4] * dm.m[-1]
            return dm

    # TODO: Add docstrings
    def _fit(
        self,
        var_names,
        conn,
        queue,
        use_raw: bool = False,
        load_pars: bool = False,
        steady_state_prior: Optional[List[bool]] = None,
        assignment_mode: str = "projection",
        **kwargs,
    ):
        """TODO."""
        idx, dms = [], []
        for gene in var_names:
            dm = DynamicsRecovery(
                self._adata,
                gene,
                use_raw=use_raw,
                load_pars=load_pars,
                max_iter=self._max_iter,
                fit_time=self._fit_time,
                fit_steady_states=self._fit_steady_states,
                fit_connected_states=conn,
                fit_scaling=self._fit_scaling,
                fit_basal_transcription=self._fit_basal_transcription,
                steady_state_prior=steady_state_prior,
                **kwargs,
            )
            if dm.recoverable:
                dm.fit(assignment_mode=assignment_mode)

                ix = np.where(self._adata.var_names == gene)[0][0]
                idx.append(ix)
                dms.append(dm)
            else:
                logg.warn(dm.gene, "not recoverable due to insufficient samples.")

            if queue is not None:
                queue.put(1)

        if queue is not None:
            queue.put(None)

        return idx, dms


def _flatten(iterable):
    return [i for it in iterable for i in it]
