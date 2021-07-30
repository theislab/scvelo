from typing import Optional

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
    set_modality,
    sum,
)
from scvelo.core._anndata import obs_df
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
    def test_equal_obs_id_length(self, obs_names, obs_names_cleaned):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        clean_obs_names(adata)

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
    def test_different_obs_id_length(self, obs_names, obs_names_cleaned):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        clean_obs_names(adata)

        assert (adata.obs_names == obs_names_cleaned).all()
        assert "sample_batch" in adata.obs
        assert adata.obs["sample_batch"].str.startswith("sample").all()


class TestCleanup(TestBase):
    @given(adata=get_adata(), inplace=st.booleans())
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

    @given(adata=get_adata(), inplace=st.booleans())
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
        adata=get_adata(layer_keys=["unspliced", "spliced", "Ms", "Mu", "random"]),
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
        adata=get_adata(),
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
            layer_keys=["unspliced", "spliced", "ambiguous"],
        ),
        layer=st.text(min_size=1, max_size=5),
    )
    def test_not_existing_modality(self, adata: AnnData, layer: str):
        initial_size = get_initial_size(adata=adata, layer=layer)

        assert initial_size is None

    @given(
        adata=get_adata(
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
    @given(adata=get_adata())
    def test_get_modality(self, adata: AnnData):
        modality_to_get = self._subset_modalities(adata, 1)[0]
        modality_retrieved = get_modality(adata=adata, modality=modality_to_get)

        if modality_to_get == "X":
            assert_array_equal(adata.X, modality_retrieved)
        elif modality_to_get in adata.layers:
            assert_array_equal(adata.layers[modality_to_get], modality_retrieved)
        else:
            assert_array_equal(adata.obsm[modality_to_get], modality_retrieved)

    @given(adata=get_adata())
    def test_modality_equals_none(self, adata: AnnData):
        modality_retrieved = get_modality(adata=adata, modality=None)

        assert_array_equal(adata.X, modality_retrieved)


class TestGetSize(TestBase):
    @given(adata=get_adata())
    def test_get_size(self, adata: AnnData):
        modality = self._subset_modalities(adata, n_modalities=1)[0]

        np.testing.assert_allclose(
            sum(get_modality(adata=adata, modality=modality), axis=1),
            get_size(adata=adata, modality=modality),
        )

    @given(adata=get_adata())
    def test_modality_set_to_none(self, adata: AnnData):
        np.testing.assert_allclose(
            sum(adata.X, axis=1),
            get_size(adata=adata, modality=None),
        )


class TestMakeDense(TestBase):
    @given(
        adata=get_adata(sparse_entries=True),
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


class TestMakeSparse(TestBase):
    @given(
        adata=get_adata(),
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


class TestObsDf(TestBase):
    @given(data=st.data(), adata=get_adata())
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


class TestSetModality(TestBase):
    @given(adata=get_adata(), inplace=st.booleans())
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
