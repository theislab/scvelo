from typing import List

import pytest
from hypothesis import given
from hypothesis import strategies as st

import numpy as np

from anndata import AnnData

from scvelo.datasets import simulation


class TestSimulation:
    @given(
        n_obs=st.integers(min_value=5, max_value=300),
        n_vars=st.integers(min_value=5, max_value=300),
        t_max=st.floats(
            min_value=1, max_value=50, allow_nan=False, allow_infinity=False
        ),
        alpha=st.floats(
            min_value=0, max_value=10, allow_nan=False, allow_infinity=False
        ),
        beta=st.floats(
            min_value=0.1, max_value=10, allow_nan=False, allow_infinity=False
        ),
        gamma=st.floats(
            min_value=0.1, max_value=10, allow_nan=False, allow_infinity=False
        ),
        noise_level=st.floats(
            min_value=0, max_value=5, allow_nan=False, allow_infinity=False
        ),
    )
    def test_normal_noise(
        self,
        n_obs: int,
        n_vars: int,
        t_max: float,
        alpha: float,
        beta: float,
        gamma: float,
        noise_level: float,
    ):
        if beta == gamma:
            beta += 1e-3
        adata = simulation(
            n_obs=n_obs,
            n_vars=n_vars,
            t_max=t_max,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            noise_level=noise_level,
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (n_obs, n_vars)
        assert len(adata.layers) == 2
        assert set(adata.layers) == {"unspliced", "spliced"}

        assert len(adata.obs.columns) == 1
        assert adata.obs.columns.isin(["true_t"]).all()
        assert adata.obs["true_t"].max() == np.round(t_max, 2)

        assert len(adata.var.columns) == 5
        assert adata.var.columns.isin(
            ["true_t_", "true_alpha", "true_beta", "true_gamma", "true_scaling"]
        ).all()
        assert (adata.var["true_alpha"] == alpha).all()
        assert (adata.var["true_beta"] == beta).all()
        assert (adata.var["true_gamma"] == gamma).all()
        assert (adata.var["true_scaling"] == 1).all()

    def test_time_dependent_parameters(self):
        adata = simulation(
            n_obs=5,
            alpha=np.array([5, 4, 0, 0, 0]),
            beta=np.array([0.5, 0.3, 0.6, 0.4, 0.7]),
            gamma=np.array([0.25, 0.4, 0.5, 0.2, 0.4]),
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (5, 4)

        assert len(adata.var.columns) == 5
        assert (
            adata.var.columns
            == ["true_t_", "true_alpha", "true_beta", "true_gamma", "true_scaling"]
        ).all()
        assert adata.var["true_alpha"].isna().sum() == 4
        assert adata.var["true_beta"].isna().sum() == 4
        assert adata.var["true_gamma"].isna().sum() == 4

    @pytest.mark.parametrize("noise_model", ("gillespie", "normal"))
    @pytest.mark.parametrize(
        "switches",
        ([0.25, 0.5, 0.75, 1], [0.3, 0.4, 0.5], [0.01, 0.2, 0.4, 0.3, 0.61, 0.7]),
    )
    def test_switch(self, noise_model, switches: List[float]):
        def ceil(x, precision=0):
            return np.true_divide(np.ceil(x * 10**precision), 10**precision)

        adata = simulation(t_max=1, switches=switches, noise_model=noise_model)

        assert isinstance(adata, AnnData)
        np.testing.assert_equal(
            ceil(adata.var["true_t_"].values, precision=2), switches
        )

    @pytest.mark.parametrize("n_obs", (5, 10, 100))
    @pytest.mark.parametrize("n_vars", (5, 10, 100))
    @pytest.mark.parametrize("t_max", (1, 10))
    @pytest.mark.parametrize("alpha", (5, 7))
    @pytest.mark.parametrize("beta", (0.5, 0.3))
    @pytest.mark.parametrize("gamma", (0.4, 0.2))
    def test_gillespie(
        self,
        n_obs: int,
        n_vars: int,
        t_max: float,
        alpha: float,
        beta: float,
        gamma: float,
    ):
        adata = simulation(
            n_obs=n_obs,
            n_vars=n_vars,
            t_max=t_max,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            noise_model="gillespie",
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (n_obs, n_vars)
        assert len(adata.layers) == 2
        assert set(adata.layers) == {"unspliced", "spliced"}

        np.testing.assert_equal(adata.X % 1, 0)
        np.testing.assert_equal(adata.layers["unspliced"] % 1, 0)
        np.testing.assert_equal(adata.layers["spliced"] % 1, 0)

        assert len(adata.obs.columns) == 1
        assert adata.obs.columns.isin(["true_t"]).all()
        assert adata.obs["true_t"].max() == np.round(t_max, 2)

        assert len(adata.var.columns) == 5
        assert adata.var.columns.isin(
            ["true_t_", "true_alpha", "true_beta", "true_gamma", "true_scaling"]
        ).all()
        assert (adata.var["true_alpha"] == alpha).all()
        assert (adata.var["true_beta"] == beta).all()
        assert (adata.var["true_gamma"] == gamma).all()
        assert (adata.var["true_scaling"] == 1).all()
