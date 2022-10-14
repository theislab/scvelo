from typing import List

import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import ndarray
from scipy.integrate import odeint

from scvelo.core import SplicingDynamics


class TestSplicingDynamics:
    @given(
        alpha=st.floats(min_value=0, allow_infinity=False),
        beta=st.floats(min_value=0, max_value=1, exclude_min=True),
        gamma=st.floats(min_value=0, max_value=1, exclude_min=True),
        initial_state=st.lists(
            st.floats(min_value=0, allow_infinity=False), min_size=2, max_size=2
        ),
        t=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        with_keys=st.booleans(),
    )
    def test_output_form(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: List[float],
        t: ndarray,
        with_keys: bool,
    ):
        if beta == gamma:
            gamma = gamma + 1e-6

        splicing_dynamics = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )
        solution = splicing_dynamics.get_solution(t=t, with_keys=with_keys)

        if not with_keys:
            assert type(solution) == ndarray
            assert solution.shape == (len(t), 2)
        else:
            assert len(solution) == 2
            assert type(solution) == dict
            assert list(solution.keys()) == ["u", "s"]
            assert all([len(var) == len(t) for var in solution.values()])

    # TODO: Check how / if hypothesis can be used instead.
    @pytest.mark.parametrize(
        "alpha, beta, gamma, initial_state",
        [
            (5, 0.5, 0.4, [0, 1]),
        ],
    )
    def test_solution(self, alpha, beta, gamma, initial_state):
        def model(y, t, alpha, beta, gamma):
            dydt = np.zeros(2)
            dydt[0] = alpha - beta * y[0]
            dydt[1] = beta * y[0] - gamma * y[1]

            return dydt

        t = np.linspace(0, 20, 10000)
        splicing_dynamics = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )
        exact_solution = splicing_dynamics.get_solution(t=t)

        numerical_solution = odeint(
            model,
            np.array(initial_state),
            t,
            args=(
                alpha,
                beta,
                gamma,
            ),
        )

        assert np.allclose(numerical_solution, exact_solution)

    @pytest.mark.parametrize(
        "alpha, beta, gamma, initial_state",
        [
            (5, 0.5, 0.4, [0, 1]),
        ],
    )
    def test_2d_time(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: List[float],
    ):
        def model(y, t, alpha, beta, gamma):
            dydt = np.zeros(2)
            dydt[0] = alpha - beta * y[0]
            dydt[1] = beta * y[0] - gamma * y[1]

            return dydt

        t = np.linspace(0, 20, 10000)
        t = np.vstack([t, t, t]).T

        splicing_dynamics = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )
        exact_solution = splicing_dynamics.get_solution(t=t)

        assert exact_solution.shape == (*t.shape, 2)

        numerical_solution = np.stack(
            [
                odeint(
                    model,
                    np.array(initial_state),
                    t[:, col_id],
                    args=(
                        alpha,
                        beta,
                        gamma,
                    ),
                )
                for col_id in range(t.shape[1])
            ],
            axis=1,
        )

        assert np.allclose(numerical_solution, exact_solution)

    @given(
        alpha=st.floats(),
        beta=st.floats(),
        gamma=st.floats(),
        initial_state=arrays(
            float,
            shape=st.integers(min_value=2, max_value=2),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    def test_intitial_state_1d(self, alpha, beta, gamma, initial_state):
        dm = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )

        np.testing.assert_array_equal(dm.initial_state, initial_state)

        dm.initial_state = np.array([0, 0])
        np.testing.assert_array_equal(dm.initial_state, np.array([0, 0]))

    @given(
        alpha=st.floats(),
        beta=st.floats(),
        gamma=st.floats(),
        initial_state=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=2, max_value=2),
            ),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    def test_intitial_state_2d(self, alpha, beta, gamma, initial_state):
        dm = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )

        np.testing.assert_array_equal(dm.initial_state, initial_state)

        dm.initial_state = np.zeros(2)
        np.testing.assert_array_equal(dm.initial_state, np.zeros(2))

    @given(
        alpha=st.floats(allow_infinity=False, allow_nan=False),
        beta=st.floats(min_value=1e-10, allow_infinity=False, allow_nan=False),
        gamma=st.floats(min_value=1e-10, allow_infinity=False, allow_nan=False),
        initial_state=arrays(
            float,
            shape=st.integers(min_value=2, max_value=2),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        with_keys=st.booleans(),
        stacked=st.booleans(),
    )
    def test_steady_state_1d(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: ndarray,
        with_keys: bool,
        stacked: bool,
    ):
        dm = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )
        steady_states = dm.get_steady_states(stacked=stacked, with_keys=with_keys)

        if with_keys:
            assert isinstance(steady_states, dict)
            assert [*steady_states] == ["u", "s"]
            assert steady_states["u"] == alpha / beta
            assert steady_states["s"] == alpha / gamma
        elif not stacked:
            assert isinstance(steady_states, tuple)
            assert len(steady_states) == 2
            assert steady_states[0] == alpha / beta
            assert steady_states[1] == alpha / gamma
        else:
            assert isinstance(steady_states, np.ndarray)
            assert steady_states.shape == (2,)
            assert steady_states[0] == alpha / beta
            assert steady_states[1] == alpha / gamma

    @given(
        alpha=st.floats(),
        beta=st.floats(max_value=0, allow_infinity=False, allow_nan=False),
        gamma=st.floats(max_value=0, allow_infinity=False, allow_nan=False),
        initial_state=arrays(
            float,
            shape=st.integers(min_value=2, max_value=2),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        with_keys=st.booleans(),
        stacked=st.booleans(),
    )
    def test_steady_state_not_defined(
        self,
        alpha: float,
        beta: float,
        gamma: float,
        initial_state: ndarray,
        with_keys: bool,
        stacked: bool,
    ):
        dm = SplicingDynamics(
            alpha=alpha, beta=beta, gamma=gamma, initial_state=initial_state
        )

        with pytest.raises(
            ValueError, match=("Both `beta` and `gamma` need to be strictly positive.")
        ):
            _ = dm.get_steady_states(stacked=stacked, with_keys=with_keys)
