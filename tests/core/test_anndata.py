from typing import Dict, List, Optional

import hypothesis.strategies as st
import pytest
from hypothesis import given

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal
from scipy.sparse import issparse

from anndata import AnnData

from scvelo.core import (
    clean_obs_names,
    cleanup,
    get_initial_size,
    get_modality,
    get_size,
    make_dense,
    make_sparse,
    set_initial_size,
    set_modality,
    show_proportions,
    sum,
)
from scvelo.core._anndata import obs_df, var_df
from .test_base import get_adata, TestBase


# TODO: Make more sophisticated
class TestCleanObsNames:
    @pytest.mark.parametrize(
        "obs_names, obs_names_cleaned",
        [
            (
                ["sample1_ABCD", "sample2_ABCD", "sample3_DCBA"],
                ["ABCD", "ABCD-1", "DCBA"],
            ),
            (
                ["sample1_ABCD0815", "sample2_AMNC0707", "sample3_AAAA0902"],
                ["ABCD", "AMNC", "AAAA"],
            ),
        ],
    )
    @pytest.mark.parametrize("inplace", (True, False))
    def test_equal_obs_id_length(
        self, obs_names: List[str], obs_names_cleaned: List[str], inplace: bool
    ):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        _adata = clean_obs_names(adata, inplace=inplace)

        if inplace:
            assert _adata is None
        else:
            assert isinstance(_adata, AnnData)
            adata = _adata

        assert (adata.obs_names == obs_names_cleaned).all()
        assert "sample_batch" in adata.obs
        assert adata.obs["sample_batch"].str.startswith("sample").all()

    @pytest.mark.parametrize(
        "obs_names, obs_names_cleaned",
        [
            (
                ["sample1_ABCDE0815", "sample2_AMNC0707", "sample3_AAAA0902"],
                ["ABCD", "AMNC", "AAAA"],
            )
        ],
    )
    @pytest.mark.parametrize("inplace", (True, False))
    def test_different_obs_id_length(
        self, obs_names: List[str], obs_names_cleaned: List[str], inplace: bool
    ):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        _adata = clean_obs_names(adata, inplace=inplace)

        if inplace:
            assert _adata is None
        else:
            assert isinstance(_adata, AnnData)
            adata = _adata

        assert (adata.obs_names == obs_names_cleaned).all()
        assert "sample_batch" in adata.obs
        assert adata.obs["sample_batch"].str.startswith("sample").all()


class TestCleanup(TestBase):
    @given(adata=get_adata(max_obs=5, max_vars=5), inplace=st.booleans())
    def test_cleanup_all(self, adata: AnnData, inplace: bool):
        returned_adata = cleanup(adata=adata, clean="all", inplace=inplace)

        if not inplace:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        assert len(adata.layers) == 0
        assert len(adata.uns) == 0
        assert len(adata.obs.columns) == 0
        assert len(adata.var.columns) == 0

    @given(adata=get_adata(max_obs=5, max_vars=5), inplace=st.booleans())
    def test_cleanup_default_clean_w_random_adata(self, adata: AnnData, inplace: bool):
        n_obs_cols = len(adata.obs.columns)
        n_var_cols = len(adata.var.columns)
        n_uns_slots = len(adata.uns)

        returned_adata = cleanup(adata=adata, inplace=inplace)
        if not inplace:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        assert len(adata.layers) == 0
        assert len(adata.uns) == n_uns_slots
        assert len(adata.obs.columns) == n_obs_cols
        assert len(adata.var.columns) == n_var_cols

    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced", "Ms", "Mu", "random"],
        ),
        inplace=st.booleans(),
    )
    def test_cleanup_default_clean(self, adata: AnnData, inplace: bool):
        n_obs_cols = len(adata.obs.columns)
        n_var_cols = len(adata.var.columns)
        n_uns_slots = len(adata.uns)

        returned_adata = cleanup(adata=adata, inplace=inplace)

        if not inplace:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        assert len(adata.layers) == 4
        assert len(adata.uns) == n_uns_slots
        assert len(adata.obs.columns) == n_obs_cols
        assert len(adata.var.columns) == n_var_cols

    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        inplace=st.booleans(),
        n_modalities=st.integers(min_value=0),
        n_cols=st.integers(min_value=0),
    )
    def test_cleanup_some(
        self, adata: AnnData, inplace: bool, n_modalities: int, n_cols: int
    ):
        layers_to_keep = self._subset_modalities(
            adata,
            n_modalities,
            from_obsm=False,
        )
        obs_cols_to_keep = self._subset_columns(adata, n_cols=n_cols, from_var=False)
        var_cols_to_keep = self._subset_columns(adata, n_cols=n_cols, from_obs=False)

        # Update in case adata.layers, adata.obs, adata.var share same keys
        layers_to_keep += set(adata.layers).intersection(obs_cols_to_keep)
        layers_to_keep += set(adata.layers).intersection(var_cols_to_keep)

        obs_cols_to_keep += set(adata.obs.columns).intersection(var_cols_to_keep)
        obs_cols_to_keep += set(adata.obs.columns).intersection(layers_to_keep)

        var_cols_to_keep += set(adata.var.columns).intersection(obs_cols_to_keep)
        var_cols_to_keep += set(adata.var.columns).intersection(layers_to_keep)

        returned_adata = cleanup(
            adata=adata,
            keep=layers_to_keep + obs_cols_to_keep + var_cols_to_keep,
            clean="all",
            inplace=inplace,
        )

        if not inplace:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        # Distinction is needed since `self._subset_modalities` always includes `"X"`
        if "X" in adata.layers:
            assert set(adata.layers.keys()) == set(layers_to_keep)
        else:
            assert set(adata.layers.keys()) == set(layers_to_keep).difference({"X"})
        assert set(adata.obs.columns) == set(obs_cols_to_keep)
        assert set(adata.var.columns) == set(var_cols_to_keep)


class TestGetInitialSize(TestBase):
    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced", "ambiguous"],
            obs_col_names=[
                "initial_size",
                "initial_size_unspliced",
                "initial_size_spliced",
                "initial_size_ambiguous",
            ],
        ),
        by_total_size=st.booleans(),
        layer=st.sampled_from([None, "X", "unspliced", "spliced", "ambiguous"]),
    )
    def test_get_initial_size(
        self, adata: AnnData, layer: Optional[None], by_total_size: bool
    ):
        initial_size = get_initial_size(
            adata=adata, layer=layer, by_total_size=by_total_size
        )

        if by_total_size:
            assert np.allclose(
                initial_size,
                adata.obs["initial_size_unspliced"] + adata.obs["initial_size_spliced"],
            )
        elif layer in adata.layers:
            assert np.allclose(initial_size, adata.obs[f"initial_size_{layer}"])
        else:
            assert np.allclose(initial_size, adata.obs["initial_size"])

    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced", "ambiguous"],
        ),
        layer=st.text(min_size=2, max_size=5),
    )
    def test_not_existing_modality(self, adata: AnnData, layer: str):
        initial_size = get_initial_size(adata=adata, layer=layer)

        assert initial_size is None

    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced", "ambiguous"],
        ),
        layer=st.sampled_from([None, "X", "unspliced", "spliced", "ambiguous"]),
    )
    def test_initial_size_not_in_adata_obs(self, adata: AnnData, layer: Optional[str]):
        initial_size = get_initial_size(adata=adata, layer=layer)

        if layer in [None, "X"]:
            np.testing.assert_allclose(initial_size, get_size(adata=adata))
        else:
            np.testing.assert_allclose(initial_size, get_size(adata=adata, layer=layer))


class TestGetModality(TestBase):
    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_get_modality(self, adata: AnnData):
        modality_to_get = self._subset_modalities(adata, 1)[0]
        modality_retrieved = get_modality(adata=adata, modality=modality_to_get)

        if modality_to_get == "X":
            assert_array_equal(adata.X, modality_retrieved)
        elif modality_to_get in adata.layers:
            assert_array_equal(adata.layers[modality_to_get], modality_retrieved)
        else:
            assert_array_equal(adata.obsm[modality_to_get], modality_retrieved)

    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_modality_equals_none(self, adata: AnnData):
        modality_retrieved = get_modality(adata=adata, modality=None)

        assert_array_equal(adata.X, modality_retrieved)


class TestGetSize(TestBase):
    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_get_size(self, adata: AnnData):
        modality = self._subset_modalities(adata, n_modalities=1)[0]

        np.testing.assert_allclose(
            sum(get_modality(adata=adata, modality=modality), axis=1),
            get_size(adata=adata, modality=modality),
        )

    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_modality_set_to_none(self, adata: AnnData):
        np.testing.assert_allclose(
            sum(adata.X, axis=1),
            get_size(adata=adata, modality=None),
        )


class TestMakeDense(TestBase):
    @given(
        adata=get_adata(max_obs=5, max_vars=5, sparse_entries=True),
        inplace=st.booleans(),
        n_modalities=st.integers(min_value=0),
    )
    def test_make_dense(self, adata: AnnData, inplace: bool, n_modalities: int):
        modalities_to_densify = self._subset_modalities(adata, n_modalities)

        returned_adata = make_dense(
            adata=adata, modalities=modalities_to_densify, inplace=inplace
        )

        if inplace:
            assert returned_adata is None
            assert np.all(
                [
                    not issparse(get_modality(adata=adata, modality=modality))
                    for modality in modalities_to_densify
                ]
            )
        else:
            assert isinstance(returned_adata, AnnData)
            assert np.all(
                [
                    not issparse(get_modality(adata=returned_adata, modality=modality))
                    for modality in modalities_to_densify
                ]
            )
            assert np.all(
                [
                    issparse(get_modality(adata=adata, modality=modality))
                    for modality in modalities_to_densify
                ]
            )

    @given(
        adata=get_adata(max_obs=5, max_vars=5, sparse_entries=True),
        inplace=st.booleans(),
    )
    def test_modalities_passed_as_string(self, adata: AnnData, inplace: bool):
        modality_to_densify = self._subset_modalities(adata, n_modalities=1)[0]

        returned_adata = make_dense(
            adata=adata, modalities=modality_to_densify, inplace=inplace
        )

        if inplace:
            assert returned_adata is None
            assert not issparse(get_modality(adata=adata, modality=modality_to_densify))
        else:
            assert isinstance(returned_adata, AnnData)
            assert not issparse(
                get_modality(adata=returned_adata, modality=modality_to_densify)
            )
            assert issparse(get_modality(adata=adata, modality=modality_to_densify))


class TestMakeSparse(TestBase):
    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        inplace=st.booleans(),
        n_modalities=st.integers(min_value=0),
    )
    def test_make_sparse(self, adata: AnnData, inplace: bool, n_modalities: int):
        modalities_to_make_sparse = self._subset_modalities(adata, n_modalities)

        returned_adata = make_sparse(
            adata=adata, modalities=modalities_to_make_sparse, inplace=inplace
        )

        if inplace:
            assert returned_adata is None
            assert np.all(
                [
                    issparse(get_modality(adata=adata, modality=modality))
                    for modality in modalities_to_make_sparse
                    if modality != "X"
                ]
            )
        else:
            assert isinstance(returned_adata, AnnData)
            assert np.all(
                [
                    issparse(get_modality(adata=returned_adata, modality=modality))
                    for modality in modalities_to_make_sparse
                    if modality != "X"
                ]
            )
            assert np.all(
                [
                    not issparse(get_modality(adata=adata, modality=modality))
                    for modality in modalities_to_make_sparse
                    if modality != "X"
                ]
            )

    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        inplace=st.booleans(),
    )
    def test_modalities_passed_as_string(self, adata: AnnData, inplace: bool):
        modality_to_make_sparse = self._subset_modalities(adata, n_modalities=1)[0]

        returned_adata = make_sparse(
            adata=adata, modalities=modality_to_make_sparse, inplace=inplace
        )

        if inplace:
            assert returned_adata is None
            if modality_to_make_sparse != "X":
                assert issparse(
                    get_modality(adata=adata, modality=modality_to_make_sparse)
                )
        else:
            assert isinstance(returned_adata, AnnData)
            if modality_to_make_sparse != "X":
                assert issparse(
                    get_modality(adata=returned_adata, modality=modality_to_make_sparse)
                )
                assert not issparse(
                    get_modality(adata=adata, modality=modality_to_make_sparse)
                )


class TestObsDf(TestBase):
    @given(data=st.data(), adata=get_adata(max_obs=5, max_vars=5))
    def test_obs_df(self, data, adata: AnnData):
        adata.var_names = "var_" + adata.var_names

        modality = self._subset_modalities(adata, n_modalities=1, from_obsm=False)[0]

        var_names = data.draw(
            st.lists(
                st.sampled_from(adata.var_names.to_list()),
                max_size=len(adata.var_names),
                unique=True,
            )
        )

        if modality == "X":
            df = obs_df(adata=adata, keys=var_names)
        else:
            df = obs_df(adata=adata, keys=var_names, layer=modality)

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == var_names).all()
        if len(var_names) == 0:
            assert df.shape == (adata.n_obs, 0)
        else:
            np.testing.assert_equal(
                df.values, get_modality(adata[:, var_names], modality=modality)
            )
        assert (df.index == adata.obs_names).all()

    @pytest.mark.parametrize(
        "var_names", (["var_1", "var_2"], ["var_0", "Var_1", "var_2"])
    )
    def test_warning_for_nonexisting_var_names(self, capfd, var_names):
        adata = AnnData(np.eye(len(var_names)), var=pd.DataFrame(index=var_names))

        df = obs_df(adata=adata, keys=var_names + ["VAR_1", "VAR_2"])

        actual_warning, _ = capfd.readouterr()
        expected_warning = (
            "WARNING: Keys ['VAR_1', 'VAR_2'] were not found in `adata.var_names`.\n"
        )

        assert actual_warning == expected_warning
        assert isinstance(df, pd.DataFrame)
        assert (df.index == adata.obs_names).all()


class TestSetInitialSize(TestBase):
    @given(
        adata=get_adata(max_obs=5, max_vars=5), n_modalities=st.integers(min_value=0)
    )
    def test_added_columns(self, adata: AnnData, n_modalities: int):
        layers = self._subset_modalities(
            adata=adata, n_modalities=n_modalities, from_obsm=False
        )

        set_initial_size(adata=adata, layers=layers)

        if "X" in layers and "X" not in adata.layers:
            assert (
                sum(
                    adata.obs.columns.isin(
                        [f"initial_size_{layer}" for layer in layers]
                    )
                )
                == len(layers) - 1
            )
        else:
            assert sum(
                adata.obs.columns.isin([f"initial_size_{layer}" for layer in layers])
            ) == len(layers)

        assert "initial_size" in adata.obs.columns

    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_non_existing_columns_specified(self, adata: AnnData):
        layers = "_" + adata.obs.columns
        set_initial_size(adata=adata, layers=layers)

        assert "initial_size" in adata.obs.columns
        assert len(adata.obs.columns) == 3

    @given(adata=get_adata(max_obs=5, max_vars=5, layer_keys=["unspliced", "spliced"]))
    def test_layers_not_specified(self, adata: AnnData):
        set_initial_size(adata=adata)

        assert "initial_size" in adata.obs.columns
        assert "initial_size_unspliced" in adata.obs.columns
        assert "initial_size_spliced" in adata.obs.columns
        assert adata.obs.columns.str.startswith("initial_size").sum() == 3

    @pytest.mark.parametrize(
        "X, layers, initial_size",
        [
            (
                np.eye(2),
                {"unspliced": np.ones((2, 2)), "spliced": np.array([[1, 2], [3, 3]])},
                {
                    "X": np.ones(2),
                    "unspliced": 2 * np.ones(2),
                    "spliced": np.array([3, 6]),
                },
            )
        ],
    )
    def test_calculated_initial_size(
        self, X: np.ndarray, layers: np.ndarray, initial_size: np.ndarray
    ):
        adata = AnnData(X=X, layers=layers)
        set_initial_size(adata=adata, layers=["unspliced", "spliced"])

        np.testing.assert_equal(adata.obs["initial_size"], initial_size["X"])
        np.testing.assert_equal(
            adata.obs["initial_size_unspliced"], initial_size["unspliced"]
        )
        np.testing.assert_equal(
            adata.obs["initial_size_spliced"], initial_size["spliced"]
        )


class TestSetModality(TestBase):
    @given(adata=get_adata(max_obs=5, max_vars=5), inplace=st.booleans())
    def test_set_modality(self, adata: AnnData, inplace: bool):
        modality_to_set = self._subset_modalities(adata, 1)[0]

        if (modality_to_set == "X") or (modality_to_set in adata.layers):
            new_value = np.random.randn(adata.n_obs, adata.n_vars)
        else:
            new_value = np.random.randn(
                adata.n_obs, np.random.randint(low=1, high=10000)
            )

        returned_adata = set_modality(
            adata=adata, new_value=new_value, modality=modality_to_set, inplace=inplace
        )

        if inplace:
            assert returned_adata is None
            if modality_to_set == "X":
                assert_array_equal(adata.X, new_value)
            elif modality_to_set in adata.layers:
                assert_array_equal(adata.layers[modality_to_set], new_value)
            else:
                assert_array_equal(adata.obsm[modality_to_set], new_value)
        else:
            assert isinstance(returned_adata, AnnData)
            if modality_to_set == "X":
                assert_array_equal(returned_adata.X, new_value)
            elif modality_to_set in adata.layers:
                assert_array_equal(returned_adata.layers[modality_to_set], new_value)
            else:
                assert_array_equal(returned_adata.obsm[modality_to_set], new_value)


class TestShowProportions(TestBase):
    @pytest.mark.parametrize(
        "layers",
        (
            {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
            {
                "unspliced": np.eye(2),
                "spliced": 2 * np.eye(2),
                "ambiguous": 3 * np.eye(2),
            },
            {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
            {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
        ),
    )
    @pytest.mark.parametrize("use_raw", (True, False))
    def test_layers_not_specified(self, capfd, layers: Dict, use_raw: bool):
        adata = AnnData(X=np.eye(2), layers=layers)

        show_proportions(adata=adata, layers=None, use_raw=use_raw)
        actual_output, _ = capfd.readouterr()

        if len(layers) == 2:
            expected_output = f"Abundance of {[*layers]}: [0.33 0.67]\n"
        else:
            expected_output = f"Abundance of {[*layers]}: [0.17 0.33 0.5 ]\n"

        assert actual_output == expected_output

    @pytest.mark.parametrize(
        "layers",
        (
            {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
            {
                "unspliced": np.eye(2),
                "spliced": 2 * np.eye(2),
                "ambiguous": 3 * np.eye(2),
            },
            {"layer_1": np.eye(2), "layer_2": 2 * np.eye(2)},
        ),
    )
    @pytest.mark.parametrize("use_raw", (True, False))
    def test_layers_specified(self, capfd, layers: Dict, use_raw: bool):
        adata = AnnData(X=np.eye(2), layers=layers)

        show_proportions(adata=adata, layers=layers.keys(), use_raw=use_raw)
        actual_output, _ = capfd.readouterr()

        if len(layers) == 2:
            expected_output = f"Abundance of {[*layers]}: [0.33 0.67]\n"
        else:
            expected_output = f"Abundance of {[*layers]}: [0.17 0.33 0.5 ]\n"

        assert actual_output == expected_output

    @pytest.mark.parametrize(
        "layers",
        (
            {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
            {
                "unspliced": np.eye(2),
                "spliced": 2 * np.eye(2),
                "ambiguous": 3 * np.eye(2),
            },
            {"layer_1": np.eye(2), "layer_2": 2 * np.eye(2)},
        ),
    )
    @pytest.mark.parametrize("use_raw", (True, False))
    def test_passing_nonexisting_layers(self, capfd, layers: Dict, use_raw: bool):
        adata = AnnData(X=np.eye(2), layers=layers)

        show_proportions(
            adata=adata, layers=[*layers] + ["random_1", "random_2"], use_raw=use_raw
        )
        actual_output, _ = capfd.readouterr()

        if len(layers) == 2:
            expected_output = f"Abundance of {[*layers]}: [0.33 0.67]\n"
        else:
            expected_output = f"Abundance of {[*layers]}: [0.17 0.33 0.5 ]\n"

        assert actual_output == expected_output

    @pytest.mark.parametrize(
        "layers, obs",
        (
            (
                {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
                {
                    "initial_size_unspliced": np.ones(2),
                    "initial_size_spliced": np.ones(2),
                },
            ),
            (
                {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
                {
                    "initial_size_unspliced": np.ones(2),
                    "initial_size_spliced": np.ones(2),
                },
            ),
            (
                {"unspliced": np.eye(2), "spliced": 2 * np.eye(2)},
                {"initial_size_unspliced": np.ones(2)},
            ),
        ),
    )
    @pytest.mark.parametrize("use_raw", (True, False))
    def test_initial_size_specified(
        self, capfd, layers: Dict, obs: Dict, use_raw: bool
    ):
        adata = AnnData(X=np.eye(2), layers=layers, obs=obs)

        show_proportions(adata=adata, layers=[*layers], use_raw=use_raw)
        actual_output, _ = capfd.readouterr()

        if len(adata.obs.columns) == 2:
            if use_raw:
                expected_output = f"Abundance of {[*layers]}: [0.5 0.5]\n"
            else:
                expected_output = f"Abundance of {[*layers]}: [0.33 0.67]\n"
        else:
            expected_output = f"Abundance of {[*layers]}: [0.33 0.67]\n"

        assert actual_output == expected_output


class TestVarDf(TestBase):
    @given(data=st.data(), adata=get_adata())
    def test_var_df(self, data, adata: AnnData):
        adata.obs_names = "obs_" + adata.obs_names

        modality = self._subset_modalities(adata, n_modalities=1, from_obsm=False)[0]

        obs_names = data.draw(
            st.lists(
                st.sampled_from(adata.obs_names.to_list()),
                max_size=len(adata.obs_names),
                unique=True,
            )
        )

        if modality == "X":
            df = var_df(adata=adata, keys=obs_names)
        else:
            df = var_df(adata=adata, keys=obs_names, layer=modality)

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == obs_names).all()
        if len(obs_names) == 0:
            assert df.shape == (adata.n_vars, 0)
        else:
            np.testing.assert_equal(
                df.values, get_modality(adata[obs_names, :], modality=modality).T
            )
        assert (df.index == adata.var_names).all()

    @pytest.mark.parametrize(
        "obs_names", (["obs_1", "obs_2"], ["obs_0", "Obs_1", "obs_2"])
    )
    def test_warning_for_nonexisting_obs_names(self, capfd, obs_names):
        adata = AnnData(np.eye(len(obs_names)), obs=pd.DataFrame(index=obs_names))

        df = var_df(adata=adata, keys=obs_names + ["OBS_1", "OBS_2"])

        actual_warning, _ = capfd.readouterr()
        expected_warning = (
            "WARNING: Keys ['OBS_1', 'OBS_2'] were not found in `adata.obs_names`.\n"
        )

        assert actual_warning == expected_warning
        assert isinstance(df, pd.DataFrame)
        assert (df.index == adata.var_names).all()
