from typing import List, Optional, Tuple, Union

import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix, issparse

from ._arithmetic import prod_sum, sum


class LinearRegression:
    """Extreme quantile and constraint least square linear regression.

    Arguments:
    ---------
    percentile
        Percentile of data on which linear regression line is fit. If `None`, all data
        is used, if a single value is given, it is interpreted as the upper quantile.
        Defaults to `None`.
    fit_intercept
        Whether to calculate the intercept for model. Defaults to `False`.
    positive_intercept
        Whether the intercept it constraint to positive values. Only plays a role when
        `fit_intercept=True`. Defaults to `True`.
    constrain_ratio
        Ratio to which coefficients are clipped. If `None`, the coefficients are not
        constraint. Defaults to `None`.

    Attributes
    ----------
    `coef_`
        Estimated coefficients of the linear regression line.

    `intercept_`
        Fitted intercept of linear model. Set to `0.0` if `fit_intercept=False`.

    """

    def __init__(
        self,
        percentile: Optional[Union[Tuple, int, float]] = None,
        fit_intercept: bool = False,
        positive_intercept: bool = True,
        constrain_ratio: Optional[Union[Tuple, float]] = None,
    ):
        if not fit_intercept and isinstance(percentile, (list, tuple)):
            self.percentile = percentile[1]
        else:
            self.percentile = percentile
        self.fit_intercept = fit_intercept
        self.positive_intercept = positive_intercept

        if constrain_ratio is None:
            self.constrain_ratio = [-np.inf, np.inf]
        elif len(constrain_ratio) == 1:
            self.constrain_ratio = [-np.inf, constrain_ratio]
        else:
            self.constrain_ratio = constrain_ratio

    def _trim_data(self, data: List) -> List:
        """Trim data to extreme values.

        Arguments:
        ---------
        data
            Data to be trimmed to extreme quantiles.

        Returns
        -------
        List
            Number of non-trivial entries per column and trimmed data.
        """
        if not isinstance(data, List):
            data = [data]

        data = np.array(
            [data_mat.A if issparse(data_mat) else data_mat for data_mat in data]
        )

        # TODO: Add explanatory comment
        normalized_data = np.sum(
            data / data.max(axis=1, keepdims=True).clip(1e-3, None), axis=0
        )

        bound = np.percentile(normalized_data, self.percentile, axis=0)

        if bound.ndim == 1:
            trimmer = csr_matrix(normalized_data >= bound).astype(bool)
        else:
            trimmer = csr_matrix(
                (normalized_data <= bound[0]) | (normalized_data >= bound[1])
            ).astype(bool)

        return [trimmer.getnnz(axis=0)] + [
            trimmer.multiply(data_mat).tocsr() for data_mat in data
        ]

    def fit(self, x: ndarray, y: ndarray):
        """Fit linear model per column.

        Arguments:
        ---------
        x
            Training data of shape `(n_obs, n_vars)`.
        y
            Target values of shape `(n_obs, n_vars)`.

        Returns
        -------
        self
            Returns: an instance of self.
        """
        n_obs = x.shape[0]

        if self.percentile is not None:
            n_obs, x, y = self._trim_data(data=[x, y])

        _xx = prod_sum(x, x, axis=0)
        _xy = prod_sum(x, y, axis=0)

        if self.fit_intercept:
            _x = sum(x, axis=0) / n_obs
            _y = sum(y, axis=0) / n_obs
            self.coef_ = (_xy / n_obs - _x * _y) / (_xx / n_obs - _x**2)
            self.intercept_ = _y - self.coef_ * _x

            if self.positive_intercept:
                idx = self.intercept_ < 0
                if self.coef_.ndim > 0:
                    self.coef_[idx] = _xy[idx] / _xx[idx]
                else:
                    self.coef_ = _xy / _xx
                self.intercept_ = np.clip(self.intercept_, 0, None)
        else:
            self.coef_ = _xy / _xx
            self.intercept_ = np.zeros(x.shape[1]) if x.ndim > 1 else 0

        if not np.isscalar(self.coef_):
            self.coef_[np.isnan(self.coef_)] = 0
            self.intercept_[np.isnan(self.intercept_)] = 0
        else:
            if np.isnan(self.coef_):
                self.coef_ = 0
            if np.isnan(self.intercept_):
                self.intercept_ = 0

        self.coef_ = np.clip(self.coef_, *self.constrain_ratio)

        return self
