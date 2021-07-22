import hypothesis.strategies as st
import pytest
from hypothesis import given

import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse import issparse

from anndata import AnnData

from scvelo.core import (
    clean_obs_names,
    cleanup,
    get_modality,
    make_dense,
    make_sparse,
    set_modality,
)
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
    @given(adata=get_adata(), copy=st.booleans())
    def test_cleanup_all(self, adata: AnnData, copy: bool):
        returned_adata = cleanup(adata, clean="all", copy=copy)

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        assert len(adata.layers) == 0
        assert len(adata.uns) == 0
        assert len(adata.obs.columns) == 0
        assert len(adata.var.columns) == 0

    @given(adata=get_adata(), copy=st.booleans())
    def test_cleanup_default_clean_w_random_adata(self, adata: AnnData, copy: bool):
        n_obs_cols = len(adata.obs.columns)
        n_var_cols = len(adata.var.columns)
        n_uns_slots = len(adata.uns)

        returned_adata = cleanup(adata)
        assert returned_adata is None

        assert len(adata.layers) == 0
        assert len(adata.uns) == n_uns_slots
        assert len(adata.obs.columns) == n_obs_cols
        assert len(adata.var.columns) == n_var_cols

    @given(
        adata=get_adata(layer_keys=["unspliced", "spliced", "Ms", "Mu", "random"]),
        copy=st.booleans(),
    )
    def test_cleanup_default_clean(self, adata: AnnData, copy: bool):
        n_obs_cols = len(adata.obs.columns)
        n_var_cols = len(adata.var.columns)
        n_uns_slots = len(adata.uns)

        returned_adata = cleanup(adata, copy=copy)

        if copy:
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
        copy=st.booleans(),
        n_modalities=st.integers(min_value=0),
        n_cols=st.integers(min_value=0),
    )
    def test_cleanup_some(
        self, adata: AnnData, copy: bool, n_modalities: int, n_cols: int
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
        obs_cols_to_keep += set(adata.var.columns).intersection(layers_to_keep)

        returned_adata = cleanup(
            adata,
            keep=layers_to_keep + obs_cols_to_keep + var_cols_to_keep,
            clean="all",
            copy=copy,
        )

        if copy:
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
