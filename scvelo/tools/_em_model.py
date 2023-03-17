import os
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from anndata import AnnData

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import get_n_jobs, parallelize
from scvelo.preprocessing.moments import get_connectivities
from ._core import BaseInference
from ._em_model_core import (
    _flatten,
    _read_pars,
    _write_pars,
    align_dynamics,
    DynamicsRecovery,
)
from ._steady_state_model import SteadyStateModel
from .utils import make_unique_list


@dataclass
class EMParams:
    """EM parameters."""

    alpha: np.ndarray
    beta: np.ndarray
    gamma: np.ndarray
    t_: np.ndarray
    scaling: np.ndarray
    std_u: np.ndarray
    std_s: np.ndarray
    likelihood: np.ndarray
    u0: np.ndarray
    s0: np.ndarray
    pval_steady: np.ndarray
    steady_u: np.ndarray
    steady_s: np.ndarray
    variance: np.ndarray


# TODO: Refactor to use `EMParams`
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
                self.r2 = velo.state_dict()["r2"]
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

    # TODO: Implement
    def state_dict(self):
        """Return the state of the model."""
        raise NotImplementedError

    # TODO: Implement
    def export_results_adata(self):
        """Export the results to the AnnData object."""
        raise NotImplementedError

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
        **kwargs,
    ):
        """Fit the model."""
        logg.info(
            f"recovering dynamics (using {self._n_jobs}/{os.cpu_count()} cores)", r=True
        )

        # TODO: Remove or move to `__init__`
        if len(set(self._adata.var_names)) != len(self._adata.var_names):
            logg.warn("Duplicate var_names found. Making them unique.")
            self._adata.var_names_make_unique()

        if (
            "Ms" not in self._adata.layers.keys()
            or "Mu" not in self._adata.layers.keys()
        ):
            use_raw = True
        if self._fit_connected_states is None:
            fit_connected_states = not use_raw

        self._prepare_genes()
        if return_model is None:
            return_model = len(self._var_names) < 5

        pars = _read_pars(self._adata)
        alpha, beta, gamma, t_, scaling, std_u, std_s, likelihood = pars[:8]
        u0, s0, pval, steady_u, steady_s, varx = pars[8:]
        # likelihood[np.isnan(likelihood)] = 0
        idx, L = [], []
        T = np.zeros(self._adata.shape) * np.nan
        Tau = np.zeros(self._adata.shape) * np.nan
        Tau_ = np.zeros(self._adata.shape) * np.nan
        if f"{self._fit_key}_t" in self._adata.layers.keys():
            T = self._adata.layers[f"{self._fit_key}_t"]
        if f"{self._fit_key}_tau" in self._adata.layers.keys():
            Tau = self._adata.layers[f"{self._fit_key}_tau"]
        if f"{self._fit_key}_tau_" in self._adata.layers.keys():
            Tau_ = self._adata.layers[f"{self._fit_key}_tau_"]

        conn = get_connectivities(self._adata) if fit_connected_states else None

        res = parallelize(
            self._fit,
            self._var_names,
            self._n_jobs,
            unit="gene",
            as_array=False,
            backend=self._backend,
            show_progress_bar=len(self._var_names) > 9,
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
            T[:, ix], Tau[:, ix], Tau_[:, ix] = dm.t, dm.tau, dm.tau_
            alpha[ix], beta[ix], gamma[ix], t_[ix], scaling[ix] = dm.pars[:, -1]
            u0[ix], s0[ix], pval[ix] = dm.u0, dm.s0, dm.pval_steady
            steady_u[ix], steady_s[ix] = dm.steady_u, dm.steady_s
            beta[ix] /= scaling[ix]
            steady_u[ix] *= scaling[ix]

            std_u[ix], std_s[ix] = dm.std_u, dm.std_s
            likelihood[ix], varx[ix] = dm.likelihood, dm.varx
            L.append(dm.loss)

        _pars = [
            alpha,
            beta,
            gamma,
            t_,
            scaling,
            std_u,
            std_s,
            likelihood,
            u0,
            s0,
            pval,
            steady_u,
            steady_s,
            varx,
        ]

        adata = self._adata.copy() if copy else self._adata

        adata.uns["recover_dynamics"] = {
            "fit_connected_states": fit_connected_states,
            "fit_basal_transcription": self._fit_basal_transcription,
            "use_raw": use_raw,
        }

        _write_pars(adata, _pars, add_key=self._fit_key)
        if f"{self._fit_key}_t" in adata.layers.keys():
            adata.layers[f"{self._fit_key}_t"][:, idx] = (
                T[:, idx] if conn is None else conn.dot(T[:, idx])
            )
        else:
            adata.layers[f"{self._fit_key}_t"] = T if conn is None else conn.dot(T)
        adata.layers[f"{self._fit_key}_tau"] = Tau
        adata.layers[f"{self._fit_key}_tau_"] = Tau_

        # is False if only one invalid / irrecoverable gene was given in var_names
        if L:
            cur_len = adata.varm["loss"].shape[1] if "loss" in adata.varm.keys() else 2
            max_len = max(np.max([len(loss) for loss in L]), cur_len) if L else cur_len
            loss = np.ones((adata.n_vars, max_len)) * np.nan

            if "loss" in adata.varm.keys():
                loss[:, :cur_len] = adata.varm["loss"]

            loss[idx] = np.vstack(
                [
                    np.concatenate([loss, np.ones(max_len - len(loss)) * np.nan])
                    for loss in L
                ]
            )
            adata.varm["loss"] = loss

        # TODO: Fix s.t. `self._t_max` is only integer
        if self._t_max is not False:
            dm = align_dynamics(adata, t_max=self._t_max, dm=dm, idx=idx)

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
