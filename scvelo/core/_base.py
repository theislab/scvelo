from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union

from numpy import ndarray


class DynamicsBase(ABC):
    @abstractmethod
    def get_solution(
        self, t: ndarray, stacked: True, with_keys: bool = False
    ) -> Union[Dict, Tuple[ndarray], ndarray]:
        """Calculate solution of dynamics.

        Arguments
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
        Union[Dict, Tuple[ndarray], ndarray]
            Solution of system. If `with_keys=True`, the solution is returned in form of
            a dictionary with variables as keys. Otherwise, the solution is given as
            a `numpy.ndarray` of form `(n_steps, n_vars)`.
        """

        return

    @abstractmethod
    def get_steady_states(
        self, stacked: True, with_keys: False
    ) -> Union[Dict[str, ndarray], Tuple[ndarray], ndarray]:
        """Return steady state of system.

        Arguments
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

        return
