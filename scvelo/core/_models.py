from typing import Dict, List, Union

import numpy as np
from numpy import ndarray

from ._base import DynamicsBase


# TODO: Improve parameter names: alpha -> transcription_rate; beta -> splicing_rate;
# gamma -> degradation_rate
class SplicingDynamics(DynamicsBase):
    """Splicing dynamics.

    Arguments
    ---------
    alpha:
        Transcription rate.

    beta:
        Translation rate.

    gamma:
        Splicing degradation rate.

    initial_state:
        Initial state of system. Defaults to `[0, 0]`.

    Attributes
    ----------

    """

    def __init__(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: Union[List, ndarray] = [0, 0],
    ):

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        self.initial_state = initial_state

    @property
    def initial_state(self):
        return self._initial_state

    @initial_state.setter
    def initial_state(self, val):
        if isinstance(val, list) or (isinstance(val, ndarray) and (val.ndim == 1)):
            self.u0 = val[0]
            self.s0 = val[1]
        else:
            self.u0 = val[:, 0]
            self.s0 = val[:, 1]
        self._initial_state = val

    # TODO: If parameters are not arrays, broadcast them to right dimension and work
    # with unspliced[beta == 0 & self.gamma == 0] = ..., etc
    def get_solution(self, t: ndarray, with_keys: bool = False) -> Union[Dict, ndarray]:
        """Calculate solution of dynamics.

        Arguments
        ---------
        t:
            Time steps at which to evaluate solution.

        with_keys:
            Whether to return solution labelled by variables in form of a dictionary.
            Defaults to `False`.

        Returns
        -------
        Union[Dict, ndarray]
            Solution of system. If `with_keys=True`, the solution is returned in form of
            a dictionary with variables as keys. Otherwise, the solution is given as
            a `numpy.ndarray` of form `(n_steps, 2)`.

        """

        beta = self.beta if not isinstance(self.beta, np.ndarray) else self.beta[0]
        gamma = self.gamma if not isinstance(self.gamma, np.ndarray) else self.gamma[0]

        expu = np.exp(-self.beta * t)
        exps = np.exp(-self.gamma * t)

        if beta == 0:
            unspliced = self.u0 + self.alpha * t
            spliced = self.s0 * exps
        elif gamma == 0:
            unspliced = self.u0 * expu + self.alpha / self.beta * (1 - expu)
            spliced = (
                self.s0
                + self.alpha * t
                + (self.u0 - self.alpha / self.beta) * (1 - expu)
            )
        elif beta == gamma:
            unspliced = self.u0 * expu + self.alpha / self.beta * (1 - expu)
            spliced = (
                self.s0 * exps
                + (self.beta * self.u0 - self.alpha) * t * exps
                + self.alpha / self.gamma * (1 - exps)
            )
        else:
            unspliced = self.u0 * expu + self.alpha / self.beta * (1 - expu)
            c = (self.alpha - self.u0 * self.beta) / (self.gamma - self.beta)
            spliced = (
                self.s0 * exps
                + self.alpha / self.gamma * (1 - exps)
                + c * (exps - expu)
            )

        if with_keys:
            return {"u": unspliced, "s": spliced}
        else:
            return np.column_stack([unspliced, spliced])

    # TODO: Handle cases `beta = 0`, `gamma = 0`
    def get_steady_states(self) -> ndarray:
        """Return steady state of system.

        Arguments
        ---------

        Returns
        -------
        ndarray:
            Steady state of system.

        """

        if (self.beta > 0) and (self.gamma > 0):
            return np.array([self.alpha / self.beta, self.alpha / self.gamma])
