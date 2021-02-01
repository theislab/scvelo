from abc import ABC, abstractmethod
from typing import Dict, Union

from numpy import ndarray


class DynamicsBase(ABC):
    @abstractmethod
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
            a `numpy.ndarray` of form `(n_steps, n_vars)`.

        """

        return

    @abstractmethod
    def get_steady_states(self) -> ndarray:
        """Return steady state of system.

        Arguments
        ---------

        Returns
        -------
        ndarray:
            Steady state of system.

        """

        return
