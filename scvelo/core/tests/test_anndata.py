import hypothesis.strategies as st
from hypothesis import given

import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse import issparse

from anndata import AnnData

from scvelo.core import get_modality, make_dense, make_sparse, set_modality
from .test_base import get_adata, TestBase


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
