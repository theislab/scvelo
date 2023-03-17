from dataclasses import asdict, dataclass
from typing import Dict, Literal, Optional, Tuple, Union

import numpy as np

from anndata import AnnData

from scvelo import logging as logg
from scvelo.core import LinearRegression
from scvelo.preprocessing.moments import second_order_moments
from ._core import BaseInference
from .optimization import leastsq_generalized, maximum_likelihood
from .utils import make_dense, R_squared


@dataclass
class SteadyStateParams:
    """Steady state parameters."""

    residual: Tuple[np.ndarray, None]
    offset: np.ndarray
    gamma: np.ndarray
    qreg_ratio: np.ndarray
    r2: np.ndarray
    beta: np.ndarray
    velocity_genes: np.ndarray


@dataclass
class SecondOrderSteadyStateParams(SteadyStateParams):
    """Second order steady state parameters."""

    offset2: np.ndarray
    gamma2: np.ndarray
    residual2: Tuple[np.ndarray, None]


# TODO: Make equivalent to calling `scvelo.tools.velocity(adata, mode="deterministic")`
# TODO: Save velocity (residual) to `state_dict`
class SteadyStateModel(BaseInference):
    """Steady-state model for velocity estimation.

    Parameters
    ----------
    adata
        Annotated data matrix.
    spliced_layer
        Key for spliced layer.
    unspliced_layer
        Key for unspliced layer.
    groups_for_fit
        Subset of cells to fit the model, e.g. ['g1', 'g2', 'g3'].
    groupby
        Key of observations grouping to consider in adata.obs.
    constrain_ratio
        Bounds for the steady-state ratio.
    min_r2
        Minimum threshold for coefficient of determination.
    min_ratio
        Minimum threshold for steady-state ratio
    use_highly_variable
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    r2_adjusted
        Whether to compute coefficient of determination
        on full data fit (adjusted) or extreme quantile fit.
    fit_offset
        Whether to fit with offset for first order moment dynamics.
    perc
        Percentile for extreme quantile fit.
    """

    def __init__(
        self,
        adata: AnnData,
        spliced_layer: str = "Ms",
        unspliced_layer: str = "Mu",
        constrain_ratio: Optional[Union[float, Tuple[float, float]]] = None,
        min_r2: float = 0.01,
        min_ratio: float = 0.01,
        use_highly_variable: bool = True,
        r2_adjusted: bool = True,
        fit_offset: bool = False,
        perc: Tuple[int, int] = (5, 95),
    ):
        super().__init__(adata)
        self._Ms = adata.layers[spliced_layer]
        self._Mu = adata.layers[unspliced_layer]
        self._Ms, self._Mu = make_dense(self._Ms), make_dense(self._Mu)
        self._constrain_ratio = constrain_ratio
        self._r2_adjusted = r2_adjusted
        self._min_r2 = min_r2
        self._min_ratio = min_ratio
        self._highly_variable = None
        if use_highly_variable is not None and adata is not None:
            if "highly_variable" in adata.var.keys():
                self._highly_variable = adata.var["highly_variable"]
        self._fit_offset = fit_offset
        self._perc = perc
        _, self.n_vars = self._Ms.shape

    def _initialize_state_dict(self):
        self._state_dict = SteadyStateParams(
            residual=None,
            offset=np.zeros(self.n_vars),
            gamma=np.zeros(self.n_vars),
            qreg_ratio=np.zeros(self.n_vars),
            r2=np.zeros(self.n_vars),
            beta=np.zeros(self.n_vars),
            velocity_genes=np.zeros(self.n_vars, dtype=bool),
        )

    def state_dict(self) -> Dict[str, np.ndarray]:
        """Return the state of the model."""
        return asdict(self._state_dict)

    def load_state_dict(self, state_dict: Dict[str, np.ndarray]):
        """Load the state of the model."""
        self._state_dict = SteadyStateParams(**state_dict)

    def fit(self) -> None:
        """Fit the model."""
        if self._state_dict is None:
            self._initialize_state_dict()
        sd = self._state_dict
        Ms = self._Ms
        Mu = self._Mu

        lr = LinearRegression(
            fit_intercept=self._fit_offset, percentile=list(self._perc)
        )
        lr.fit(Ms, Mu)
        sd.offset = lr.intercept_
        sd.gamma = lr.coef_

        if self._constrain_ratio is not None:
            if np.size(self._constrain_ratio) < 2:
                self._constrain_ratio = [None, self._constrain_ratio]
            cr = self._constrain_ratio
            sd.gamma = np.clip(sd.gamma, cr[0], cr[1])

        self.residual = self._Mu - sd.gamma * self._Ms
        if self._fit_offset:
            self.residual -= sd.offset
        residual = self.residual

        # velocity genes
        if self._r2_adjusted:
            lr = LinearRegression(fit_intercept=self._fit_offset)
            lr.fit(Ms, Mu)
            offset = lr.intercept_
            gamma = lr.coef_

            residual = self._Mu - gamma * self._Ms
            if self._fit_offset:
                residual -= offset

        self.qreg_ratio = np.array(sd.gamma)  # quantile regression ratio

        sd.r2 = R_squared(residual, total=self._Mu - self._Mu.mean(0))
        sd.velocity_genes = (
            (sd.r2 > self._min_r2)
            & (sd.gamma > self._min_ratio)
            & (np.max(self._Ms > 0, 0) > 0)
            & (np.max(self._Mu > 0, 0) > 0)
        )

        if self._highly_variable is not None:
            sd.velocity_genes &= self._highly_variable

        if np.sum(sd.velocity_genes) < 2:
            min_r2 = np.percentile(sd.r2, 80)
            sd.velocity_genes = sd.r2 > min_r2
            min_r2 = np.round(min_r2, 4)
            logg.warn(
                f"You seem to have very low signal in splicing dynamics.\n"
                f"The correlation threshold has been reduced to {min_r2}.\n"
                f"Please be cautious when interpreting results."
            )

    def get_velocity(self) -> None:
        """Get velocity from steady-state model."""
        return self._state_dict.residual

    # TODO: Actually implement result export
    def export_results_adata(self) -> None:
        """Export results to adata."""
        return super().export_results_adata()


class SecondOrderSteadyStateModel(SteadyStateModel):
    """Second-order steady-state model for velocity estimation.

    Also called the 'stochatic' steady-state model.

    Parameters
    ----------
    adata
        Annotated data matrix.
    spliced_layer
        Key for spliced layer.
    unspliced_layer
        Key for unspliced layer.
    groups_for_fit
        Subset of cells to fit the model, e.g. ['g1', 'g2', 'g3'].
    groupby
        Key of observations grouping to consider in adata.obs.
    constrain_ratio
        Bounds for the steady-state ratio.
    min_r2
        Minimum threshold for coefficient of determination.
    min_ratio
        Minimum threshold for steady-state ratio
    use_highly_variable
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    r2_adjusted
        Whether to compute coefficient of determination
        on full data fit (adjusted) or extreme quantile fit.
    fit_offset
        Whether to fit with offset for first order moment dynamics.
    perc
        Percentile for extreme quantile fit.
    fit_offset2
        Whether to fit with offset for second order moment dynamics.
    mode
        Whether to use a maximimum likelihood or least squares fit.
    """

    def __init__(
        self,
        adata: AnnData,
        spliced_layer: str = "Ms",
        unspliced_layer: str = "Mu",
        constrain_ratio: Optional[Union[float, Tuple[float, float]]] = None,
        min_r2: float = 0.01,
        min_ratio: float = 0.01,
        use_highly_variable: bool = True,
        r2_adjusted: bool = False,
        fit_offset: bool = False,
        perc: Tuple[int, int] = (5, 95),
        fit_offset2: bool = False,
        mode: Literal["bayes", "ols"] = "ols",
    ):
        super().__init__(
            adata,
            spliced_layer,
            unspliced_layer,
            constrain_ratio,
            min_r2,
            min_ratio,
            use_highly_variable,
            r2_adjusted,
            fit_offset,
            perc,
        )
        self._fit_offset2 = fit_offset2
        self._mode = mode

    def _initialize_state_dict(self):
        super()._initialize_state_dict()
        self._state_dict = SecondOrderSteadyStateParams(
            offset2=np.zeros(self.n_vars),
            gamma2=np.zeros(self.n_vars),
            residual2=None,
            **self.state_dict(),
        )

    def fit(self):
        """Fit the model."""
        super().fit()
        sd = self._state_dict
        idx = np.ones(sd.velocity_genes.shape, dtype=bool)
        if np.any(sd.velocity_genes):
            idx = sd.velocity_genes
        is_subset = len(set(idx)) > 1

        _adata = self._adata[:, idx] if is_subset else self._adata
        _Ms = self._Ms[:, idx] if is_subset else self._Ms
        _Mu = self._Mu[:, idx] if is_subset else self._Mu
        _residual = self.residual[:, idx] if is_subset else sd.residual

        _Mss, _Mus = second_order_moments(_adata)

        var_ss = 2 * _Mss - _Ms
        cov_us = 2 * _Mus + _Mu

        lr = LinearRegression(fit_intercept=self._fit_offset2)
        lr.fit(var_ss, cov_us)
        _offset2 = lr.intercept_
        _gamma2 = lr.coef_

        # initialize covariance matrix
        res_std = _residual.std(0)
        res2_std = (cov_us - _gamma2 * var_ss - _offset2).std(0)

        # solve multiple regression
        sd.offset[idx], sd.offset2[idx], sd.gamma[idx] = (
            maximum_likelihood(
                _Ms, _Mu, _Mus, _Mss, self._fit_offset, self._fit_offset2
            )
            if self._mode == "bayes"
            else leastsq_generalized(
                _Ms,
                _Mu,
                var_ss,
                cov_us,
                res_std,
                res2_std,
                self._fit_offset,
                self._fit_offset2,
                list(self._perc),
            )
        )

        sd.residual = self._Mu - sd.gamma * self._Ms
        if self._fit_offset:
            sd.residual -= sd.offset

        _residual2 = (cov_us - 2 * _Ms * _Mu) - sd.gamma[idx] * (var_ss - 2 * _Ms**2)
        if self._fit_offset:
            _residual2 += 2 * sd.offset[idx] * _Ms
        if self._fit_offset2:
            _residual2 -= sd.offset2[idx]
        if is_subset:
            sd.residual2 = np.zeros(self._Ms.shape, dtype=np.float32)
            sd.residual2[:, idx] = _residual2
        else:
            sd.residual2 = _residual2

    def get_variance_velocity(self):
        """Get variance of velocity from steady-state model."""
        return self._state_dict.residual2
