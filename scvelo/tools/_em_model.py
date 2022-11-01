from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from anndata import AnnData

from scvelo.core import BaseInference, get_n_jobs
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


class ExpectationMaximizationModel(BaseInference):
    """
    EM 'Dynamical' model for velocity estimation.

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

    def fit(self):
        """Fit the model."""
        # TODO: wrap what recover dynamics does
