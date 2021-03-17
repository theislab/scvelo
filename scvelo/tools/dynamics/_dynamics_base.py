from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from numpy import ndarray
from scipy.optimize import minimize

from scvelo.tools.dynamical_model_utils import adjust_increments
from scvelo.tools.utils import test_bimodality


# TODO: Add docstrings
# TODO: Finish type hints
# TODO: Remove argument `model_parameters` and infer them from the class `dynamics`.
# TODO: Combine `initialize` and `fit` in single method
class DynamicsRecoveryBase(ABC):
    def __init__(
        self,
        dynamics,
        model_parameters,
        connectivities=None,
        percentile: Union[int, float] = 99,
        max_iter: int = 10,
        fit_time: bool = True,
        fit_scaling: Union[bool, int, float] = True,
        fit_steady_states: bool = True,
        fit_connected_states: bool = True,
        fit_basal_transcription=None,
        high_pars_resolution: bool = False,
        steady_state_prior: Optional[ndarray] = None,
        optimization_method: str = "Nelder-Mead",
    ):
        # TODO: Check how to remove this
        self.assignment_mode = None

        self.dynamics = dynamics
        self.model_parameters = model_parameters
        self.connectivities = connectivities

        self.fit_basal_transcription = fit_basal_transcription
        self.fit_connected_states = fit_connected_states
        self.fit_time = fit_time
        self.fit_scaling = fit_scaling
        self.fit_steady_states = fit_steady_states
        self.high_pars_resolution = high_pars_resolution
        self.max_iter = max_iter
        self.percentile = percentile
        self.steady_state_prior = steady_state_prior

        self.optimization_kwargs = {
            "method": optimization_method,
            "options": {"maxiter": int(self.max_iter / 5)},
        }

    def initialize(
        self,
        counts,
        weighted: bool = True,
        initial_parameter_fit: Optional[Union[Dict, List[Dict]]] = None,
    ):
        if initial_parameter_fit is None:
            initial_parameter_fit = []
        elif isinstance(initial_parameter_fit, dict):
            initial_parameter_fit = [initial_parameter_fit]

        counts = counts.copy().astype("float")

        # TODO: Remove attribute counts
        self.raw_counts = counts.copy()

        weights = (counts > 0).all(axis=1)
        self.recoverable = weights.sum() > 2

        if not self.recoverable:
            return

        if weighted:
            ub = np.percentile(counts[weights, :], self.percentile, axis=0)
            weights &= (counts[:, ub > 0] <= ub[ub > 0]).all(axis=1)

        self.weights_ = weights.copy()
        self.std_ = counts[weights, :].std(axis=0)
        if any(self.std_ == 0):
            self.std_ = np.ones(len(self.std_))

        self.weights_upper_ = weights
        if weights.any():
            self.weights_upper_ &= (counts > (counts[weights, :].max(axis=0) / 3)).all(
                axis=1
            )

        weighted_counts = counts[self.weights_, :]

        if isinstance(self.fit_scaling, bool) and self.fit_scaling:
            self.scaling = self.std_ / self.std_[-1]
        elif isinstance(self.fit_scaling, bool) and not self.fit_scaling:
            self.scaling = np.ones(len(self.std_))
        elif isinstance(self.fit_scaling, (int, float)):
            self.scaling = np.array([self.fit_scaling, 1])
            self.fit_scaling = True
        else:
            self.scaling = np.array(self.fit_scaling)
            self.fit_scaling = True

        counts /= self.scaling
        weighted_counts /= self.scaling

        if self.steady_state_prior is None:
            self.steady_state_prior = np.array([False] * len(self.weights_))

        # TODO: Check why try-except is needed and which error is thrown
        try:
            _, pvals, means = zip(
                *[
                    test_bimodality(weighted_counts[:, col_id], kde=True)
                    for col_id in range(counts.shape[1])
                ]
            )
            means = np.array(means).T
        except ValueError:
            pvals = np.ones(counts.shape[1])
            means = np.zeros((counts.shape[1], counts.shape[1]))

        self.pval_steady = max(pvals)
        self.steady_state = means[1, :]

        self._initialize_parameters(weighted_counts)

        model_params = self.get_model_parameters()
        self.t_ = self.tau_inv(state=self.initial_state_[None, :], **model_params)

        model_params.update({"t_": self.t_, "scaling": self.scaling[0]})
        self.params = pd.DataFrame(model_params, index=[0])

        self.t, self.tau, self.o = self.get_time_assignment()

        self.loss_ = [self.get_loss()]

        for kwargs in initial_parameter_fit:
            self.initial_parameter_fit(**kwargs)

        self._set_steady_state_ratio(**self.get_model_parameters())

        return self

    # TODO: Add case `weighted == "outer"`
    def get_weights(self, weighted=False, weights_cluster=None):
        if not weighted:
            weights = np.ones(len(self.weights_), bool)
        elif weighted == "upper":
            weights = self.weights_upper_
        else:
            weights = self.weights_

        if weights_cluster is not None and len(weights) == len(weights_cluster):
            weights &= weights_cluster
        return weights

    # TODO: Try to remove function
    def get_counts(self, scaling=None, weighted=False, weights_cluster=None):
        scaling = self.scaling if scaling is None else scaling

        counts = self.raw_counts / scaling

        if weighted or weights_cluster is not None:
            weights = self.get_weights(
                weighted=weighted, weights_cluster=weights_cluster
            )
            counts = counts[weights, :]

        return counts

    def get_vars(self, scaling=None, t_=None, initial_state_=None, **model_parameters):
        scaling = self.scaling if scaling is None else scaling
        if t_ is None or t_ == 0:
            t_ = (
                self.t_
                if initial_state_ is None
                else self.tau_inv(
                    state=initial_state_,
                    **self.get_model_parameters(**model_parameters),
                )
            )
        return scaling, t_

    def get_model_parameters(self, after_switch=False, **kwargs):
        model_parameters = {
            model_parameter: kwargs.get(model_parameter, getattr(self, model_parameter))
            for model_parameter in self.model_parameters
        }

        if after_switch:
            self._set_parameters_after_switch(model_parameters)

        return model_parameters

    @abstractmethod
    def _set_parameters_after_switch(self, model_parameters):
        return model_parameters

    def set_model_parameters(self, **model_parameters):
        for parameter, value in model_parameters.items():
            setattr(self, parameter, value)

    def _get_parameter_dict(self, parameter_names, parameter_values):
        parameter_dict = dict(zip(parameter_names, parameter_values))
        if "scaling" in parameter_dict:
            parameter_dict["scaling"] = np.append(parameter_dict["scaling"], 1)

        return parameter_dict

    def get_variance(self, regularize=False, **kwargs):
        if "weighted" not in kwargs:
            kwargs.update({"weighted": "upper"})
        residuals = self.get_residuals(**kwargs)
        sum_squares = (residuals ** 2).sum(axis=1)
        if regularize:
            sum_squares += self.get_regularization(**kwargs) ** 2
        return (
            np.mean(sum_squares)
            - np.mean(np.sign(residuals[:, -1]) * np.sqrt(sum_squares)) ** 2
        )

    def get_loglikelihood(self, varx=None, noise_model="normal", **kwargs):
        dist = self.get_distances(**kwargs)
        n = np.clip(len(dist) - self.raw_counts.shape[0] * 0.01, 2, None)

        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        if varx is None:
            varx = self.get_variance(regularize=True, **kwargs)
        varx += varx == 0  # edge case of mRNAs levels to be the same across all cells

        if noise_model == "normal":
            return -1 / 2 / n * np.sum(dist ** 2) / varx - 1 / 2 * np.log(
                2 * np.pi * varx
            )
        elif noise_model == "laplace":
            return -1 / np.sqrt(2) / n * np.sum(dist) / np.sqrt(varx) - 1 / 2 * np.log(
                2 * varx
            )

    def get_likelihood(self, **kwargs):
        if "weighted" not in kwargs:
            kwargs.update({"weighted": "upper"})

        return np.exp(self.get_loglikelihood(**kwargs))

    @abstractmethod
    def _initialize_parameters(self, weighted_counts):
        pass

    @abstractmethod
    def _set_steady_state_ratio(self, **model_params):
        pass

    def initial_parameter_fit(self, parameter_names, sight=0.5, num=4):
        if isinstance(parameter_names, str):
            parameter_names = [parameter_names]

        if "scaling" in parameter_names:
            # TODO: Check for better implementation / solution
            scaling_values = (
                self.scaling
                + np.hstack(
                    [
                        np.linspace(
                            [-1] * (len(self.initial_state_) - 1),
                            [1] * (len(self.initial_state_) - 1),
                            num=num,
                        ),
                        np.zeros((num, 1)),
                    ]
                )
                * self.scaling
                * sight
            )
            for scaling_value in scaling_values:
                params = {
                    parameter: getattr(self, parameter)
                    / self.scaling[0]
                    * scaling_value[0]
                    for parameter in parameter_names
                    if parameter != "scaling"
                }
                self.update(scaling=scaling_value, **params)
        else:
            parameter_values = [
                getattr(self, parameter)
                + np.linspace(-1, 1, num=num) * getattr(self, parameter) * sight
                for parameter in parameter_names
            ]
            for itr in [
                dict(zip(parameter_names, value)) for value in zip(*parameter_values)
            ]:
                self.update(**itr)

    # TODO: Add rescale_factor
    def get_time_assignment(
        self,
        scaling=None,
        t_=None,
        initial_state_=None,
        t=None,
        refit_time=None,
        weighted=None,
        weights_cluster=None,
        **model_parameters,
    ):
        model_parameters = self.get_model_parameters(**model_parameters)

        if refit_time is None:
            refit_time = self.fit_time

        if t is not None:
            t_ = self.t_ if t_ is None else t_
            o = np.array(t < t_, dtype=int)
            tau = t * o + (t - t_) * (1 - o)
        elif refit_time:
            t, tau, o = self.get_divergence(
                scaling,
                t_,
                initial_state_=initial_state_,
                mode="assign_timepoints",
                **model_parameters,
            )
        else:
            t, tau, o = self.t, self.tau, self.o

        if weighted or weights_cluster is not None:
            weights = self.get_weights(
                weighted=weighted, weights_cluster=weights_cluster
            )
            t, tau, o = t[weights], tau[weights], o[weights]
        return t, tau, o

    # TODO: Find better name
    @abstractmethod
    def _check_projection(self, **model_parameters):
        return True

    def assign_tau(
        self, state, initial_state_, t_, assignment_mode=None, **model_parameters
    ):
        model_parameters_ = self.get_model_parameters(
            after_switch=True, **model_parameters
        )

        if assignment_mode in {"full_projection", "partial_projection"} or (
            assignment_mode == "projection"
            and self._check_projection(**model_parameters)
        ):
            t0 = self.tau_inv(
                state=state,
                initial_state=initial_state_,
                full_projection=True,
                **model_parameters_,
            )

            num = np.clip(int(state.shape[0] / 5), 200, 500)
            tpoints = np.linspace(0, t_, num=num)
            tpoints_ = np.linspace(0, t0, num=num)[1:]

            xt = self.dynamics(**model_parameters).get_solution(tpoints)
            xt_ = self.dynamics(
                initial_state=initial_state_, **model_parameters_
            ).get_solution(tpoints_)

            # assign time points (oth. projection onto 'on' and 'off' curve)
            tau = tpoints[
                ((xt[None, :, :] - state[:, None, :]) ** 2).sum(axis=2).argmin(axis=1)
            ]
            tau_ = tpoints_[
                ((xt_[None, :, :] - state[:, None, :]) ** 2).sum(axis=2).argmin(axis=1)
            ]
        else:
            tau = self.tau_inv(state=state, **model_parameters)
            tau = np.clip(tau, 0, t_)

            # TODO: Check if clipping is correct in general
            tau_ = self.tau_inv(
                state=state, initial_state=initial_state_, **model_parameters_
            )
            tau_ = np.clip(tau_, 0, np.max(tau_[state[:, -1] > 0]))

        return tau, tau_

    @abstractmethod
    def tau_inv(self, state, initial_state, **model_parameters):
        pass

    def get_divergence(
        self,
        scaling=None,
        t_=None,
        initial_state_=None,
        mode=None,
        **kwargs,
    ):
        model_parameters = self.get_model_parameters(**kwargs)
        scaling, t_ = self.get_vars(scaling, t_, initial_state_, **model_parameters)

        scaled_counts = self.raw_counts / scaling

        return self.compute_divergence(
            state=scaled_counts,
            t_=t_,
            initial_state_=initial_state_,
            std=self.std_,
            mode=mode,
            scaling=scaling,
            **model_parameters,
        )

    def compute_divergence(
        self,
        state,
        scaling=1,
        t_=None,
        initial_state_=None,
        tau=None,
        tau_=None,
        std=1,
        mode="distance",
        assignment_mode=None,
        var_scale=False,
        kernel_width=None,
        constraint_time_increments=True,
        noise_model="chi",
        **model_parameters,
    ):
        """"""

        if assignment_mode is None:
            assignment_mode = self.assignment_mode

        model_parameters_ = self.get_model_parameters(
            after_switch=True, **model_parameters
        )

        # set tau, tau_
        if initial_state_ is None:
            initial_state_ = (
                self.dynamics(**model_parameters).get_solution(t_).flatten()
            )
        if tau is None or tau_ is None:
            tau, tau_ = self.assign_tau(
                state, initial_state_, t_, assignment_mode, **model_parameters
            )

        std /= scaling

        # adjust increments of tau, tau_ to avoid meaningless jumps
        if constraint_time_increments:
            sol = self.dynamics(**model_parameters).get_solution(tau)
            sol_ = self.dynamics(
                initial_state=initial_state_, **model_parameters_
            ).get_solution(tau_)

            scaled_diff = (state - sol) / std
            scaled_diff_ = (state - sol_) / std

            res = np.array(
                [(scaled_diff_ ** 2).sum(axis=1), (scaled_diff ** 2).sum(axis=1)]
            )
            if self.connectivities is not None and self.connectivities is not False:
                res = (
                    np.array([self.connectivities.dot(r) for r in res])
                    if res.ndim > 2
                    else self.connectivities.dot(res.T).T
                )

            o = np.argmin(res, axis=0)

            off, on = o == 0, o == 1
            if np.any(on) and np.any(off):
                tau[on], tau_[off] = adjust_increments(tau[on], tau_[off])
            elif np.any(on):
                tau[on] = adjust_increments(tau[on])
            elif np.any(off):
                tau_[off] = adjust_increments(tau_[off])

        if mode == "tau":
            return [tau, tau_]

        # compute induction/repression state distances
        sol = self.dynamics(**model_parameters).get_solution(tau)
        sol_ = self.dynamics(
            initial_state=initial_state_, **model_parameters_
        ).get_solution(tau_)

        scaled_diff = (state - sol) / std
        scaled_dist_ = (state - sol_) / std

        sum_of_squares = (scaled_diff ** 2).sum(axis=1)
        sum_of_squares_ = (scaled_dist_ ** 2).sum(axis=1)

        res = np.array([sum_of_squares_, sum_of_squares])
        varx = 1

        # compute steady state distances
        if self.fit_steady_states:
            sum_of_squares_steady = np.sum(
                ((state - self.dynamics(**model_parameters).get_steady_states()) / std)
                ** 2,
                axis=1,
            )
            sum_of_squares_steady_ = ((state / std) ** 2).sum(axis=1)

            res = np.array(
                [
                    sum_of_squares_,
                    sum_of_squares,
                    sum_of_squares_steady_,
                    sum_of_squares_steady,
                ]
            )

        if self.connectivities is not None and self.connectivities is not False:
            res = (
                np.array([self.connectivities.dot(r) for r in res])
                if res.ndim > 2
                else self.connectivities.dot(res.T).T
            )

        # compute variances
        if noise_model == "chi":
            if var_scale:
                o = np.argmin([sum_of_squares_, sum_of_squares], axis=0)
                sign = np.sign(
                    scaled_diff * o.reshape(-1, 1)
                    + scaled_diff_ * (1 - o).reshape(-1, 1)
                )[:, 1]
                dist = sum_of_squares * o + sum_of_squares_ * (1 - o)
                varx = (
                    np.mean(dist, axis=0) - np.mean(sign * np.sqrt(dist), axis=0) ** 2
                )
                if kernel_width is not None:
                    varx *= kernel_width ** 2
                res /= varx
            elif kernel_width is not None:
                res /= kernel_width ** 2

        if mode in {"assign_timepoints", "time"}:
            o = np.argmin(res, axis=0)

            tau_ *= o == 0
            tau *= o == 1

            if 2 in o:
                o[o == 2] = 1
            if 3 in o:
                o[o == 3] = 0

            t = tau * (o == 1) + (tau_ + t_) * (o == 0)
            res = [t, tau, o] if mode == "assign_timepoints" else t

        return res

    @abstractmethod
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
        pass

    @abstractmethod
    def _get_regularization(self, weighted_counts, **model_parameters):
        pass

    def get_regularization(
        self, scaling=None, weighted=True, weights_cluster=None, **model_parameters
    ):
        model_parameters = self.get_model_parameters(**model_parameters)

        weighted_counts = self.get_counts(
            scaling, weighted=weighted, weights_cluster=weights_cluster
        )

        if getattr(self, "steady_state_ratio", None) is not None:
            return self._get_regularization(weighted_counts, **model_parameters)
        else:
            return 0

    def get_distances(
        self,
        t=None,
        t_=None,
        scaling=None,
        initial_state_=None,
        refit_time=None,
        weighted=True,
        weights_cluster=None,
        regularize=True,
        **model_parameters,
    ):
        model_parameters, diff = self.get_residuals(
            t=t,
            t_=t_,
            scaling=scaling,
            refit_time=refit_time,
            weighted=weighted,
            weights_cluster=weights_cluster,
            return_model_kwargs=True,
            **model_parameters,
        )

        if regularize:
            reg = self.get_regularization(
                scaling=scaling,
                weighted=weighted,
                weights_cluster=weights_cluster,
                **model_parameters,
            )
        else:
            reg = 0

        return np.sqrt((diff ** 2).sum(axis=1) + reg ** 2)

    def get_loss(
        self,
        kind="se",
        regularize=True,
        t=None,
        t_=None,
        scaling=None,
        initial_state_=None,
        refit_time=None,
        **model_parameters,
    ):
        kwargs = {
            "t": t,
            "t_": t_,
            "scaling": scaling,
            "initial_state_": initial_state_,
            "refit_time": refit_time,
        }

        if kind == "se":
            return (
                self.get_distances(regularize=regularize, **kwargs, **model_parameters)
                ** 2
            ).sum()
        elif kind == "mse":
            return (
                self.get_distances(regularize=regularize, **kwargs, **model_parameters)
                ** 2
            ).mean()

    # TODO: Add adjust_t_ part from original implementation
    def update(
        self, t=None, t_=None, scaling=None, initial_state_=None, **model_parameters
    ):
        loss_prev = self.loss_[-1] if len(self.loss_) > 0 else 1e6

        model_parameters = self.get_model_parameters(**model_parameters)
        scaling, t_ = self.get_vars(scaling, t_, initial_state_, **model_parameters)

        t, tau, o = self.get_time_assignment(
            scaling=scaling,
            t_=t_,
            initial_state_=initial_state_,
            t=t,
            **model_parameters,
        )

        loss = self.get_loss(t=t, t_=t_, scaling=scaling, **model_parameters)
        perform_update = loss < loss_prev

        if perform_update:
            if scaling is not None:
                self.steady_state *= self.scaling / scaling
                self.initial_state_ *= self.scaling / scaling
            if initial_state_ is not None:
                self.initial_state_ = initial_state_

            self.t, self.tau, self.o = t, tau, o
            self.set_model_parameters(**model_parameters)

            # BUG: Only correct for 2D systems
            if self.fit_scaling and np.isscalar(scaling):
                self.scaling = np.array([scaling, 1])
            elif self.fit_scaling:
                self.scaling = scaling
            self.t_ = t_

            if not hasattr(self, "likelihood_"):
                self.likelihood_ = [self.get_likelihood(refit_time=False)]
            self.likelihood_.append(self.get_likelihood(refit_time=False))

            new_params = {"scaling": self.scaling[0]}
            # TODO: Find better fix, i.e. make sure `self.t_` is always a scalar
            if np.isscalar(self.t_):
                new_params["t_"] = self.t_
            else:
                new_params["t_"] = self.t_[0]
            new_params.update(model_parameters)
            self.params = self.params.append(new_params, ignore_index=True)
            self.loss_.append(loss)

    # TODO: Generalize to more general scaling attributes, e.g. for a 3D model.
    def fit_parameters(self, parameters, **kwargs):
        def mse(x):
            return self.get_loss(
                kind="mse",
                **self._get_parameter_dict(
                    parameter_names=parameters, parameter_values=x
                ),
                **kwargs,
            )

        def callback(x):
            if self.high_pars_resolution:
                return self.update(**dict(zip(parameters, x)))

        # A for loop is chosen for better readability. Given the small number of
        # parameters, the speed advantage of a list comprehension is negligible
        x0 = []
        for parameter in parameters:
            if "scaling" == parameter:
                x0.append(getattr(self, "scaling")[0])
            else:
                x0.append(getattr(self, parameter))

        res = minimize(mse, x0, callback=callback, **self.optimization_kwargs)

        self.update(
            **self._get_parameter_dict(
                parameter_names=parameters, parameter_values=res.x
            )
        )

    # TODO: Use better argument names: pretrain runs everything with non-optimal time
    # assignment
    def fit(self, pretrain, train, assignment_mode=None):
        if self.max_iter > 0:
            # for comparison with exact time assignment
            if assignment_mode == "full_projection":
                self.assignment_mode = assignment_mode

            for parameters in pretrain:
                self.fit_parameters(parameters=parameters)

            self.assignment_mode = assignment_mode
            self.update(adjust_t_=False)
            for parameters in train:
                self.fit_parameters(parameters=parameters, refit_time=False)

        self.update()
        self.tau, self.tau_ = self.get_divergence(mode="tau")

        self.likelihood_.append(self.get_likelihood(refit_time=False))
        self.varx = self.get_variance()

        return self

    # TODO: Add as method to a base class for dynamical models
    @abstractmethod
    def _vectorize(self, t, t_, initial_state, sorted=False, **model_paramers):
        pass
