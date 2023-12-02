from typing import Dict, List, Tuple, Union

import numpy as np
from numpy import ndarray

from ._arithmetic import invert
from ._base import DynamicsBase


# TODO: Improve parameter names: alpha -> transcription_rate; beta -> splicing_rate;
# gamma -> degradation_rate
# TODO: Handle cases beta = 0, gamma == 0, beta == gamma
class SplicingDynamics(DynamicsBase):
    """Splicing dynamics.

    Arguments:
    ---------
    alpha
        Transcription rate.
    beta
        Translation rate.
    gamma
        Splicing degradation rate.
    initial_state
        Initial state of system. Defaults to `[0, 0]`.

    Attributes
    ----------
    alpha
        Transcription rate.
    beta
        Translation rate.
    gamma
        Splicing degradation rate.
    initial_state
        Initial state of system. Defaults to `[0, 0]`.
    u0
        Initial abundance of unspliced RNA.
    s0
        Initial abundance of spliced RNA.

    """

    def __init__(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: Union[List, ndarray] = None,
    ):
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        if initial_state is None:
            self.initial_state = [0, 0]
        else:
            self.initial_state = initial_state

    # TODO: Add docstrings
    @property
    def initial_state(self):
        """TODO."""
        return self._initial_state

    # TODO: Add docstrings
    @initial_state.setter
    def initial_state(self, val):
        """TODO."""
        if isinstance(val, list) or (isinstance(val, ndarray) and (val.ndim == 1)):
            self.u0 = val[0]
            self.s0 = val[1]
        else:
            self.u0 = val[:, 0]
            self.s0 = val[:, 1]
        self._initial_state = val

    def get_solution(
        self, t: ndarray, stacked: bool = True, with_keys: bool = False
    ) -> Union[Dict, ndarray]:
        """Calculate solution of dynamics.

        Arguments:
        ---------
        t
            Time steps at which to evaluate solution.
        stacked
            Whether to stack states or return them individually. Defaults to `True`.
        with_keys
            Whether to return solution labelled by variables in form of a dictionary.
            Defaults to `False`.

        Returns
        -------
        Union[Dict, ndarray]
            Solution of system. If `with_keys=True`, the solution is returned in form of
            a dictionary with variables as keys. Otherwise, the solution is given as
            a `numpy.ndarray` of form `(n_steps, 2)`.
        """
        expu = np.exp(-self.beta * t)
        exps = np.exp(-self.gamma * t)

        unspliced = self.u0 * expu + self.alpha / self.beta * (1 - expu)
        c = (self.alpha - self.u0 * self.beta) * invert(self.gamma - self.beta)
        spliced = (
            self.s0 * exps + self.alpha / self.gamma * (1 - exps) + c * (exps - expu)
        )

        if with_keys:
            return {"u": unspliced, "s": spliced}
        elif not stacked:
            return unspliced, spliced
        else:
            if isinstance(t, np.ndarray) and t.ndim == 2:
                return np.stack([unspliced, spliced], axis=2)
            else:
                return np.column_stack([unspliced, spliced])

    # TODO: Handle cases `beta = 0`, `gamma = 0`
    def get_steady_states(
        self, stacked: bool = True, with_keys: bool = False
    ) -> Union[Dict[str, ndarray], Tuple[ndarray], ndarray]:
        """Return steady state of system.

        Arguments:
        ---------
        stacked
            Whether to stack states or return them individually. Defaults to `True`.
        with_keys
            Whether to return solution labelled by variables in form of a dictionary.
            Defaults to `False`.

        Returns
        -------
        Union[Dict[str, ndarray], Tuple[ndarray], ndarray]
            Steady state of system.
        """
        if (self.beta <= 0) or (self.gamma <= 0):
            raise ValueError("Both `beta` and `gamma` need to be strictly positive.")
        else:
            unspliced = self.alpha / self.beta
            spliced = self.alpha / self.gamma

        if with_keys:
            return {"u": unspliced, "s": spliced}
        elif not stacked:
            return unspliced, spliced
        else:
            return np.array([unspliced, spliced])
