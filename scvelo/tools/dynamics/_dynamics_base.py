from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Union

from typing_extensions import Literal

import numpy as np
import pandas as pd
from numpy import ndarray
from scipy.optimize import minimize

from scvelo.core._base import DynamicsBase
from scvelo.tools.dynamical_model_utils import adjust_increments
from scvelo.tools.utils import test_bimodality


# TODO: Remove argument `model_parameters` and infer them from the class `dynamics`.
# TODO: Check if refit_time needs to be `None` or simply boolean
class DynamicsRecoveryBase(ABC):
    """Parameter inference of dynamical systems.

    Arguments
    ---------
    dynamics
        Dynamic system for which model parameters are inferred.
    model_parameters
        Model parmeters to infer.
    connectivities
        Connectivity matrix.
    percentile
        Percentile of cells considered when subsetting.
    max_iter
        Maximal iterations in the optimization algorithm.
    fit_time
        Boolean flag to fit time or keep time assignments or not.
    fit_scaling
        Boolean flag to fit scaling of counts or not.
    fit_steady_states
        Boolean flag to fit steady states of induction and repression phase or not.
    fit_connected_states
        Boolean flag to restrict fit to neighbors given by connectivies.
    fit_basal_transcription
        Boolean flag to fit a basal transcription or not.
    high_pars_resolution
        TODO: Add description
    steady_state_prior
        Prior on steady state of system.
    optimization_method
        Algorithm used for solving optimization problems.
    """

    def __init__(
        self,
        dynamics: DynamicsBase,
        model_parameters: List,
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

    def _initialize(
        self,
        counts: ndarray,
        subsetted: bool = True,
        initial_parameter_fit: Optional[Union[Dict, List[Dict]]] = None,
    ):
        """Initialize parameter fitting.

        Arguments
        ---------
        counts
            Count matrix used to fit parameters. Each column corresponds to a variable
            of the dynamical system. The order must be the same as returned by the
            model class `self.dynamics`.
        subsetted
            Indicator whether counts should be subsetted (based on percentile) or not.
        initial_paremter_fit
            Specification of initial parameter fit through a grid. Each fit is specified
            by the parameter names `parameter_names` to fit, the width of the grid
            `sight` and the numbe of points `num` in it. For example, consider `alpha`
            being estimated initially as `alpha=1`. Using
            `{"parameter_names": ["alpha"], "sight": 0.5, "num": 5}`, `alpha` is updated
            by the value of `array([0.5 , 0.75, 1.  , 1.25, 1.5])` resulting in the
            minimal loss.
        """

        if initial_parameter_fit is None:
            initial_parameter_fit = []
        elif isinstance(initial_parameter_fit, dict):
            initial_parameter_fit = [initial_parameter_fit]

        counts = counts.copy().astype("float")

        # TODO: Remove attribute counts
        self.raw_counts = counts.copy()

        obs_subset = (counts > 0).all(axis=1)
        self.recoverable = obs_subset.sum() > 2

        if not self.recoverable:
            return

        if subsetted:
            ub = np.percentile(counts[obs_subset, :], self.percentile, axis=0)
            obs_subset &= (counts[:, ub > 0] <= ub[ub > 0]).all(axis=1)

        self.obs_subset_ = obs_subset.copy()

        # self.std_ = counts[obs_subset, :].std(axis=0) results in different estimates
        # in simulated data
        subsetted_counts = counts[obs_subset, :]
        self.std_ = np.array(
            [subsetted_counts[:, col].std() for col in range(subsetted_counts.shape[1])]
        )
        if any(self.std_ == 0):
            self.std_ = np.ones(len(self.std_))

        self.obs_subset_upper_ = obs_subset
        if obs_subset.any():
            self.obs_subset_upper_ &= (
                counts > (counts[obs_subset, :].max(axis=0) / 3)
            ).all(axis=1)

        subsetted_counts = counts[self.obs_subset_, :]

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
        subsetted_counts /= self.scaling

        if self.steady_state_prior is None:
            self.steady_state_prior = np.array([False] * len(self.obs_subset_))

        # TODO: Check why try-except is needed and which error is thrown
        try:
            _, pvals, means = zip(
                *[
                    test_bimodality(subsetted_counts[:, col_id], kde=True)
                    for col_id in range(counts.shape[1])
                ]
            )
            means = np.array(means).T
        except ValueError:
            pvals = np.ones(counts.shape[1])
            means = np.zeros((counts.shape[1], counts.shape[1]))

        self.pval_steady = max(pvals)
        self.steady_state = means[1, :]

        self._initialize_parameters(subsetted_counts)

        model_params = self.get_model_parameters()
        self.t_ = self.get_approx_time_assignment(
            state=self.initial_state_[None, :], **model_params
        )

        model_params.update({"t_": self.t_, "scaling": self.scaling[0]})
        self.params = pd.DataFrame(model_params, index=[0])

        self.t, self.tau, self.o = self.get_time_assignment()

        self.loss_ = [self.get_loss()]

        for kwargs in initial_parameter_fit:
            self.initial_parameter_fit(**kwargs)

        self._set_steady_state_ratio(**self.get_model_parameters())

    # TODO: Add case `subsetted == "outer"`
    def get_obs_subset(
        self, subsetted: bool = False, obs_subset_cluster: Optional[List] = None
    ) -> ndarray:
        """Get subset of observations in form of a filter mask.

        Arguments
        ---------
        subsetted
            Boolean flag to indicate whether observations should be subsetted or not.
        obs_subset_cluster
            TODO: Add description.

        Returns
        -------
        ndarray
            Filter mask of observation subset.
        """

        if not subsetted:
            obs_subset = np.ones(len(self.obs_subset_), bool)
        elif subsetted == "upper":
            obs_subset = self.obs_subset_upper_
        else:
            obs_subset = self.obs_subset_

        if obs_subset_cluster is not None and len(obs_subset) == len(
            obs_subset_cluster
        ):
            obs_subset &= obs_subset_cluster
        return obs_subset

    # TODO: Try to remove function
    def get_counts(
        self,
        scaling: Optional[ndarray] = None,
        subsetted: bool = False,
        obs_subset_cluster: Optional[List] = None,
    ) -> ndarray:
        """Retrieve relevant (scaled) counts.

        Arguments
        ---------
        scaling
            Scaling factors of counts. If not specified, `self.scaling` will be used.
        subsetted
            Boolean flag to indicate whether observations should be subsetted or not.
        obs_subset_cluster
            TODO: Add description.

        Returns
        -------
        ndarray
            Scaled and subsetted observations.
        """

        scaling = self.scaling if scaling is None else scaling

        counts = self.raw_counts / scaling

        if subsetted or obs_subset_cluster is not None:
            obs_subset = self.get_obs_subset(
                subsetted=subsetted, obs_subset_cluster=obs_subset_cluster
            )
            counts = counts[obs_subset, :]

        return counts

    def get_vars(
        self,
        scaling: Optional[ndarray] = None,
        t_: Optional[ndarray] = None,
        initial_state_: Optional[ndarray] = None,
        **model_parameters,
    ) -> Tuple[ndarray, float]:
        """Retrieve relevant (scaled) counts.

        Arguments
        ---------
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        t_
            Time where system switches from induction to repression. If neither `t_` nor
            `initial_state_` are specified, `self.t_` is used.
        initial_state_
            State of system when switching between induction and repression.

        Returns
        -------
        Tuple[ndarray, float]
            Scaling and time of switching.
        """

        scaling = self.scaling if scaling is None else scaling
        if t_ is None or t_ == 0:
            t_ = (
                self.t_
                if initial_state_ is None
                else self.get_approx_time_assignment(
                    state=initial_state_,
                    **self.get_model_parameters(**model_parameters),
                )
            )
        return scaling, t_

    def get_model_parameters(self, after_switch: bool = False, **kwargs) -> Dict:
        """Retrieve model parameters from dictionary of keyword arguments.

        If a model parameter `m` is not specified in `kwargs`, it is set to `self.m`.

        Arguments
        ---------
        after_switch
            Flag to indicate if model parameters during repression should be returned.
        kwargs
            Dictionary of keyword arguments to retrieve model parameters from.

        Returns
        -------
        Dict
            Model parameters.
        """

        model_parameters = {
            model_parameter: kwargs.get(model_parameter, getattr(self, model_parameter))
            for model_parameter in self.model_parameters
        }

        if after_switch:
            self._set_parameters_after_switch(model_parameters)

        return model_parameters

    @abstractmethod
    def _set_parameters_after_switch(self, model_parameters: List) -> Dict:
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

        return model_parameters

    def set_model_parameters(self, **model_parameters):
        """Set attributes of model parameters to new value.

        Arguments
        ---------
        model_parameters
            Dictionary mapping model parameter to its new value.
        """

        for parameter, value in model_parameters.items():
            setattr(self, parameter, value)

    def _get_parameter_dict(
        self, parameter_names: List, parameter_values: List
    ) -> Dict:
        """Get dictionary mapping parameter to its value.

        Arguments
        ---------
        parameter_names
            Names of parameters.
        parameter_values
            Value of parameters.

        Returns
        -------
        Dict
            Dictionary mapping parameter to its value.
        """

        parameter_dict = dict(zip(parameter_names, parameter_values))
        if "scaling" in parameter_dict:
            parameter_dict["scaling"] = np.append(parameter_dict["scaling"], 1)

        return parameter_dict

    def get_variance(self, regularize: bool = False, **kwargs) -> float:
        """Get (regularized) variance of residuals.

        Arguments
        ---------
        regularize
            Boolean flag to regularize variance or not.
        kwargs
            Keyword arguments needed to calculate residuals.

        Returns
        -------
        float
            Variance of residuals.
        """

        if "susetted" not in kwargs:
            kwargs.update({"subsetted": "upper"})
        residuals = self.get_residuals(**kwargs)
        sum_squares = (residuals ** 2).sum(axis=1)
        if regularize:
            sum_squares += self.get_regularization(**kwargs) ** 2
        return (
            np.mean(sum_squares)
            - np.mean(np.sign(residuals[:, -1]) * np.sqrt(sum_squares)) ** 2
        )

    def get_loglikelihood(
        self,
        variance: Optional[ndarray] = None,
        noise_model: Literal["laplace", "normal"] = "normal",
        **kwargs,
    ) -> float:
        """Calculate log-likelihood.

        Arguments
        ---------
        variance
            Variance used to calculate log-likelihood.
        noise_model
            Noise model used in likelihood calculation.
        kwargs
            Keyword arguments used to a) calculate distance between prediction and
            observations, b) variance.

        Returns
        -------
        float
            Log-likelihood.
        """

        dist = self.get_distances(**kwargs)
        n = np.clip(len(dist) - self.raw_counts.shape[0] * 0.01, 2, None)

        # compute variance / equivalent to np.var(np.sign(sdiff) * np.sqrt(distx))
        if variance is None:
            variance = self.get_variance(regularize=True, **kwargs)
        variance += variance == 0

        if noise_model == "normal":
            return -1 / 2 / n * np.sum(dist ** 2) / variance - 1 / 2 * np.log(
                2 * np.pi * variance
            )
        elif noise_model == "laplace":
            return -1 / np.sqrt(2) / n * np.sum(dist) / np.sqrt(
                variance
            ) - 1 / 2 * np.log(2 * variance)

    def get_likelihood(self, **kwargs) -> float:
        """Calculate likelihood.

        Arguments
        ---------
        kwargs
            Keyword arguments passed to `self.get_loglikelihood`.

        Returns
        -------
        float
            Likelihood.
        """

        if "subsetted" not in kwargs:
            kwargs.update({"subsetted": "upper"})

        return np.exp(self.get_loglikelihood(**kwargs))

    @abstractmethod
    def _initialize_parameters(self, obs_subset_counts: ndarray):
        """Initialize model parameters.

        Arguments
        ---------
        obs_subset_counts
            Filter mask to subset observations.
        """

        pass

    @abstractmethod
    def _set_steady_state_ratio(self, **model_params):
        """Set ratio between steady states of dynamical system.

        Arguments
        ---------
        model_params
            Parameters of dynamical system and their values.
        """

        pass

    def initial_parameter_fit(
        self, parameter_names: Union[str, List], sight: float = 0.5, num: int = 4
    ):
        """Estimate parameters through grid search.

        Arguments
        ---------
        parameter_names
            Names of parameters to fit.
        sight
            Sight of grid search, i.e., width of interval.
        num
            Number of points in the grid.
        """

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
        t_: Optional[Union[float, ndarray]] = None,
        initial_state_: Optional[ndarray] = None,
        t: Optional[ndarray] = None,
        refit_time: Optional[bool] = None,
        subsetted: Optional[bool] = None,
        obs_subset_cluster: Optional[List] = None,
        **model_parameters,
    ) -> Tuple[ndarray, ndarray, ndarray]:
        """Assing time points to observations

        Arguments
        ---------
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        t_
            Time when system switches states.
        initial_state_
            State of system at switching point.
        t
            Time points assigned to observations.
        refit_time
            Boolean flag to control whether the time assignment should be refitted or
            not.
        subsetted
            Flag to indicate whether observations should be subsetted or not.
        obs_subset_cluster
            TODO: Add description.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        Tuple[ndarray, ndarray, ndarray]
            Time assigned to each observation w.r.t. global start (`t`), w.r.t. start of
            state (`tau`) and mask assigning each observation to a state.
        """

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

        if subsetted or obs_subset_cluster is not None:
            obs_subset = self.get_obs_subset(
                subsetted=subsetted, obs_subset_cluster=obs_subset_cluster
            )
            t, tau, o = t[obs_subset], tau[obs_subset], o[obs_subset]

        return t, tau, o

    # TODO: Finish docstrings
    # TODO: Find better name
    @abstractmethod
    def _check_projection(self, **model_parameters) -> bool:
        """

        Arguments
        ---------
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        """

        return True

    # TODO: More detailed description for `assignment_mode`
    def assign_tau(
        self,
        state: ndarray,
        initial_state_: ndarray,
        t_: float,
        assignment_mode: Optional[
            Literal["full_projection", "partial_projection", "projection"]
        ] = None,
        **model_parameters,
    ) -> Tuple[ndarray, ndarray]:
        """Assing time (w.r.t. state) to observations.

        Arguments
        ---------
        state
            Relevant observations of the system.
        initial_state_
            State of system at switching point.
        t_
            Time point at which system switches states.
        assignment_mode
            Mode used to assign time points.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        Tuple[ndarray, ndarray]
            Pre- and post-switch time w.r.t. system state assigned to observations.
        """

        model_parameters_ = self.get_model_parameters(
            after_switch=True, **model_parameters
        )

        if assignment_mode in {"full_projection", "partial_projection"} or (
            assignment_mode == "projection"
            and self._check_projection(**model_parameters)
        ):
            t0 = self.get_approx_time_assignment(
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
            tau = self.get_approx_time_assignment(state=state, **model_parameters)
            tau = np.clip(tau, 0, t_)

            # TODO: Check if clipping is correct in general
            tau_ = self.get_approx_time_assignment(
                state=state, initial_state=initial_state_, **model_parameters_
            )
            tau_ = np.clip(tau_, 0, np.max(tau_[state[:, -1] > 0]))

        return tau, tau_

    @abstractmethod
    def get_approx_time_assignment(
        self, state: ndarray, initial_state: ndarray, **model_parameters
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

        pass

    def get_divergence(
        self,
        scaling: Optional[ndarray] = None,
        t_: float = None,
        initial_state_: ndarray = None,
        mode: Literal["tau", "assign_timepoints", "time"] = None,
        **kwargs,
    ) -> Tuple:
        """Get divergence of system.

        This is simply a helper function to call `compute_divergence`.

        Arguments
        ---------
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        t_
            Time when system switches states.
        initial_state_
            State of system at switching point.
        mode
            Mode of divergence. See `compute_divergence` for more detailed info.
        kwargs
            Keyword arguments from which model parameters and their values are
            extracted from.

        Returns
        -------
        Tuple
            Result of calling `compute_divergence`.
        """

        model_parameters = self.get_model_parameters(**kwargs)
        scaling, t_ = self.get_vars(scaling, t_, initial_state_, **model_parameters)

        scaled_counts = self.raw_counts / scaling

        return self.compute_divergence(
            state=scaled_counts,
            t_=t_,
            initial_state_=initial_state_,
            std=self.std_.copy(),
            mode=mode,
            scaling=scaling,
            **model_parameters,
        )

    # TODO: Add more details w.r.t. `mode`
    # TODO: Add more details w.r.t. `assignment_mode`
    # TODO: Add docstrings for returned value depending on chosen mode
    def compute_divergence(
        self,
        state: ndarray,
        scaling: Union[float, ndarray] = 1,
        t_: Optional[float] = None,
        initial_state_: Optional[ndarray] = None,
        tau: Optional[ndarray] = None,
        tau_: Optional[ndarray] = None,
        std: Union[float, ndarray] = 1,
        mode: Literal["tau", "assign_timepoints", "time"] = "assign_timepoints",
        assignment_mode: Optional[
            Literal["full_projection", "partial_projection", "projection"]
        ] = None,
        var_scale: bool = False,
        kernel_width: Optional[ndarray] = None,
        constrain_time_increments: bool = True,
        noise_model: Literal["chi"] = "chi",
        **model_parameters,
    ) -> Tuple:
        """Compute divergence.

        Arguments
        ---------
        state
            Observations (subsetted and scaled) to consider.
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        t_
            Time when system switches states.
        initial_state_
            State of system at switching point.
        tau
            Time per observation w.r.t. starting time of system state.
        tau_
            Clipped time per observation w.r.t. to time of switching.
        std
            Standard deviation of observations.
        mode
            Mode to calculate divergence for.
        assignment_mode
            Assignment mode for calculating `tau` and `tau_` if not specified.
        var_scale
            TODO: Add description.
        kernel_width
            Kernel width when using chi-distribution as noise model.
        constrain_time_increments
            Boolean flag to adjust increments of tau, tau_ to avoid meaningless jumps.
        noise_model
            Noise model to use.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        Tuple
            Divergence.
        """

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
        if constrain_time_increments:
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
        variance = 1

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
                variance = (
                    np.mean(dist, axis=0) - np.mean(sign * np.sqrt(dist), axis=0) ** 2
                )
                if kernel_width is not None:
                    variance *= kernel_width ** 2
                res /= variance
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

        pass

    @abstractmethod
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

        pass

    def get_regularization(
        self,
        scaling: Optional[ndarray] = None,
        subsetted: bool = True,
        obs_subset_cluster: Optional[List] = None,
        **model_parameters,
    ) -> ndarray:
        """Get regularization term.

        Arguments
        ---------
        scaling
            Scaling for counts. If not specified, `self.scaling` will be used.
        subsetted
            Boolean flag to subset observations or not.
        obs_subset_cluster
            TODO: Add description.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        ndarray
            Regularization term.
        """

        model_parameters = self.get_model_parameters(**model_parameters)

        subsetted_counts = self.get_counts(
            scaling, subsetted=subsetted, obs_subset_cluster=obs_subset_cluster
        )

        if getattr(self, "steady_state_ratio", None) is not None:
            return self._get_regularization(subsetted_counts, **model_parameters)
        else:
            return 0

    def get_distances(
        self,
        t: Optional[ndarray] = None,
        t_: Optional[float] = None,
        scaling: Optional[ndarray] = None,
        initial_state_: Optional[ndarray] = None,
        refit_time: Optional[bool] = None,
        subsetted: bool = True,
        obs_subset_cluster: Optional[List] = None,
        regularize: bool = True,
        **model_parameters,
    ) -> ndarray:
        """Calculate distances between measured observations and their predictions.

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
        regularize
            Boolean flag to regularize distances.
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        ndarray
            (Regularized) distances between measurement and prediction.
        """

        model_parameters, diff = self.get_residuals(
            t=t,
            t_=t_,
            scaling=scaling,
            refit_time=refit_time,
            subsetted=subsetted,
            obs_subset_cluster=obs_subset_cluster,
            return_model_kwargs=True,
            **model_parameters,
        )

        if regularize:
            reg = self.get_regularization(
                scaling=scaling,
                subsetted=subsetted,
                obs_subset_cluster=obs_subset_cluster,
                **model_parameters,
            )
        else:
            reg = 0

        return np.sqrt((diff ** 2).sum(axis=1) + reg ** 2)

    def get_loss(
        self,
        kind: Literal["mse", "se"] = "se",
        regularize: bool = True,
        t: Optional[ndarray] = None,
        t_: Optional[float] = None,
        scaling: Optional[ndarray] = None,
        initial_state_: Optional[ndarray] = None,
        refit_time: Optional[bool] = None,
        **model_parameters,
    ) -> float:
        """Calculate value of loss function.

        Arguments
        ---------
        kind
            Type of loss. Supported are the mean squared error `mse` and standard
            error `se`.
        regularize
            Boolean flag to regularize loss calculation.
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
        model_parameters
            Parameters of dynamical system and their values.

        Returns
        -------
        float
            Loss for a given set of parameter estimates.
        """

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
        self,
        t: Optional[ndarray] = None,
        t_: Optional[float] = None,
        scaling: Optional[ndarray] = None,
        initial_state_: Optional[ndarray] = None,
        **model_parameters,
    ):
        """Update parameter estimates.

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
        model_parameters
            Parameters of dynamical system and their values.
        """

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
    def fit_parameters(self, parameters: List, **kwargs):
        """Find parameters minizing loss.

        Arguments
        ---------
        parameters
            Parameters to fit.
        kwargs
            Keyword arguments passed to `self.get_loss`.
        """

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
    def fit(
        self,
        counts: ndarray,
        pretrain: List[List[str]],
        train: List[List[str]],
        subsetted: bool = True,
        initial_parameter_fit: Optional[Union[Dict, List[Dict]]] = None,
        assignment_mode: Optional[
            Literal["full_projection", "partial_projection", "projection"]
        ] = None,
    ):
        """Fit dynamical model to observed data.

        Arguments
        ---------
        counts
            Counts matrix of observed values. Each column corresponds to a variable of
            the dynamical system. The order must be the same as returned by the model
            class `self.dynamics`.
        pretrain
            Parameter pairing to fit with time assignment.
        train
            Parameter pairings to fit with fixed time assignment.
        subsetted
            Indicator whether counts should be subsetted (based on percentile) or not.
        initial_paremter_fit
            Specification of initial parameter fit through a grid. Each fit is specified
            by the parameter names `parameter_names` to fit, the width of the grid
            `sight` and the numbe of points `num` in it. For example, consider `alpha`
            being estimated initially as `alpha=1`. Using
            `{"parameter_names": ["alpha"], "sight": 0.5, "num": 5}`, `alpha` is updated
            by the value of `array([0.5 , 0.75, 1.  , 1.25, 1.5])` resulting in the
            minimal loss.
        assignment_mode
            Mode used to assign time points.

        Returns
        -------
        self
        """

        self._initialize(
            counts=counts,
            subsetted=subsetted,
            initial_parameter_fit=initial_parameter_fit,
        )

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
        self.variance_ = self.get_variance()

        return self
