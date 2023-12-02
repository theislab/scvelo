from typing import Dict, List, Optional, Tuple, Union

import hypothesis.strategies as st
import pytest
from hypothesis import given, settings

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal
from scipy.sparse import csr_matrix, issparse

from anndata import AnnData

from scvelo.core import (
    clean_obs_names,
    cleanup,
    get_df,
    get_initial_size,
    get_modality,
    get_size,
    make_dense,
    make_sparse,
    merge,
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
        "obs_names, obs_names_cleaned, id_length",
        [
            (
                ["sample1_ABCD", "sample2_ABCD", "sample3_DCBA"],
                ["ABCD", "ABCD-1", "DCBA"],
                4,
            ),
            (
                ["sample1_ABCD0815", "sample2_AMNC0707", "sample3_AAAA0902"],
                ["ABCD", "AMNC", "AAAA"],
                4,
            ),
            (
                [
                    "possorted_genome_bam_H66NQ:AAAGAACAGACATATGx",
                    "possorted_genome_bam_7YCS2:AAACCCAGTCGGCTACx",
                    "possorted_genome_bam_10UXK:AAAGGATGTGAATGATx",
                ],
                ["AAAGAACAGACATATG", "AAACCCAGTCGGCTAC", "AAAGGATGTGAATGAT"],
                16,
            ),
        ],
    )
    @pytest.mark.parametrize("inplace", (True, False))
    def test_equal_obs_id_length(
        self,
        obs_names: List[str],
        obs_names_cleaned: List[str],
        id_length: int,
        inplace: bool,
    ):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        _adata = clean_obs_names(adata, inplace=inplace, id_length=id_length)

        if inplace:
            assert _adata is None
        else:
            assert isinstance(_adata, AnnData)
            adata = _adata

        assert (adata.obs_names == obs_names_cleaned).all()
        assert "sample_batch" in adata.obs
        assert adata.obs["sample_batch"].str.startswith(("sample", "possorted")).all()

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
        self,
        obs_names: List[str],
        obs_names_cleaned: List[str],
        inplace: bool,
    ):
        adata = AnnData(np.eye(3))
        adata.obs_names = obs_names

        _adata = clean_obs_names(adata, inplace=inplace, id_length=4)

        if inplace:
            assert _adata is None
        else:
            assert isinstance(_adata, AnnData)
            adata = _adata

        assert (adata.obs_names == obs_names_cleaned).all()
        assert "sample_batch" in adata.obs
        assert adata.obs["sample_batch"].str.startswith("sample").all()


class TestCleanup(TestBase):
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("layer", [None, "unspliced", "spliced"])
    @pytest.mark.parametrize("dense", [True, False])
    @pytest.mark.parametrize("inplace", [True, False])
    def test_cleanup_all(
        self, adata, dataset: str, n_obs: int, layer: bool, dense: bool, inplace: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.layers["dummy_layer"] = csr_matrix(np.eye(adata.n_obs, adata.n_vars))
        adata.uns["dummy_entry"] = {"key": csr_matrix(np.eye(5, 7))}

        if dense:
            if layer is None:
                adata.X = adata.X.A
            else:
                adata.layers[layer] = adata.layers[layer].A
        returned_adata = cleanup(adata=adata, clean="all", inplace=inplace)

        if not inplace:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        assert list(adata.layers.keys()) == ["Ms", "Mu", "spliced", "unspliced"]
        assert list(adata.uns.keys()) == ["neighbors"]
        assert len(adata.obs.columns) == 0
        assert len(adata.var.columns) == 0

    @given(adata=get_adata(max_obs=5, max_vars=5), inplace=st.booleans())
    @settings(max_examples=10, deadline=1000)
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

        assert len(adata.layers) == len(
            set(adata.layers).intersection(["unspliced", "spliced", "Mu", "Ms"])
        )
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

        assert set(adata.layers.keys()) == set(layers_to_keep).difference({"X"})
        assert set(adata.obs.columns) == set(obs_cols_to_keep)
        assert set(adata.var.columns) == set(var_cols_to_keep)


class TestGetDf:
    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1", "layer_2"],
        ),
        modality=st.sampled_from([None, "X", "layer_1", "layer_2"]),
    )
    def test_indexed_by_obs_names(self, data, adata: AnnData, modality: Optional[None]):
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        index = data.draw(
            st.lists(
                st.sampled_from(adata.obs_names.to_list()),
                min_size=1,
                max_size=len(adata.obs_names),
                unique=True,
            )
        )
        keys = data.draw(
            st.lists(
                st.sampled_from(adata.var_names.to_list()),
                min_size=1,
                max_size=len(adata.var_names),
                unique=True,
            )
        )

        df = get_df(adata, keys=keys, index=index, layer=modality)

        assert isinstance(df, pd.DataFrame)
        assert (df.index == index).all()
        assert (df.columns == keys).all()
        np.testing.assert_equal(
            get_modality(adata=adata[index, keys], modality=modality), df.values
        )

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1", "layer_2"],
        ),
        modality=st.sampled_from([None, "X", "layer_1", "layer_2"]),
    )
    def test_indexed_by_var_names(self, data, adata: AnnData, modality: Optional[None]):
        adata.obs_names = [f"obs_{obs_id}" for obs_id in adata.obs_names]
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        index = data.draw(
            st.lists(
                st.sampled_from(adata.var_names.to_list()),
                min_size=1,
                max_size=len(adata.var_names),
                unique=True,
            )
        )
        keys = data.draw(
            st.lists(
                st.sampled_from(adata.obs_names.to_list()),
                min_size=1,
                max_size=len(adata.obs_names),
                unique=True,
            )
        )

        df = get_df(adata, keys=keys, index=index, layer=modality)

        assert isinstance(df, pd.DataFrame)
        assert (df.index == index).all()
        assert (df.columns == keys).all()
        np.testing.assert_equal(
            get_modality(adata=adata[keys, index], modality=modality).T, df.values
        )

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
        ),
    )
    def test_sorted_values(self, data, adata: AnnData):
        adata.obs_names = [f"obs_{obs_id}" for obs_id in adata.obs_names]
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        sort_values = data.draw(st.sampled_from([True, False] + [*adata.var_names]))

        df = get_df(
            adata, index=adata.obs_names, keys=adata.var_names, sort_values=sort_values
        )

        assert isinstance(df, pd.DataFrame)
        if isinstance(sort_values, str):
            assert (np.diff(df[sort_values]) <= 0).all()
        elif sort_values:
            assert (np.diff(df.values[:, 0]) <= 0).all()

    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1"],
            obsm_keys=["obsm_1"],
            obs_col_names=["obs_col_1"],
            var_col_names=["var_col_1"],
        ),
    )
    def test_keys_not_present(self, adata: AnnData):
        with pytest.raises(
            ValueError,
            match=(
                "'not_existing_column_name' not found in any of obs, var, obsm, varm, "
                "uns, layers, obsp, varp."
            ),
        ):
            _ = get_df(adata, keys="not_existing_column_name")

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1"],
            obsm_keys=["obsm_1"],
            varm_keys=["varm_1"],
            var_col_names=["var_col_1"],
            obs_col_names=["obs_col_1", "obs_col_2"],
        ),
        keys_as_string=st.booleans(),
    )
    def test_from_obs(self, data, adata: AnnData, keys_as_string: bool):
        adata.obs_names = [f"obs_{obs_id}" for obs_id in adata.obs_names]
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        if keys_as_string:
            keys = data.draw(st.sampled_from([*adata.obs.columns]))
        else:
            keys = data.draw(
                st.lists(
                    st.sampled_from(adata.obs.columns.to_list()),
                    min_size=1,
                    max_size=len(adata.obs.columns),
                    unique=True,
                )
            )

        df = get_df(adata, keys=keys)

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == keys).all()
        if isinstance(keys, str):
            assert (df == adata.obs[[keys]]).values.all()
        else:
            assert (df == adata.obs[keys]).values.all()

    @pytest.mark.parametrize("uns_name", ("neighbors", "random_name"))
    def test_from_uns(self, uns_name: str):
        adata = AnnData(
            np.eye(2),
            uns={uns_name: 5 * np.eye(2)},
            obs=pd.DataFrame({"obs_col_1": ["a", "b"], "obs_col_2": [0, 0]}),
        )

        df = get_df(adata, keys=uns_name)

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(df.values, 5 * np.eye(2))

    @pytest.mark.parametrize("categorical", (True, False))
    @pytest.mark.parametrize("shape", ((2, 2), (2, 1)))
    def test_from_uns_with_categorical_column(self, categorical: bool, shape: Tuple):
        adata = AnnData(
            np.eye(*shape),
            uns={"random_name": 5 * np.eye(*shape)},
            obs=pd.DataFrame({"obs_col_1": ["a", "b"], "obs_col_2": [0, 0]}),
        )

        if categorical:
            adata.obs["obs_col_1"] = adata.obs["obs_col_1"].astype("category")

        df = get_df(adata, keys="random_name")

        assert isinstance(df, pd.DataFrame)
        if categorical:
            assert (df.index == adata.obs["obs_col_1"].values).all()
            if shape[0] == shape[1]:
                assert (df.columns == adata.obs["obs_col_1"].values).all()

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1"],
            obsm_keys=["obsm_1"],
            varm_keys=["varm_1"],
            obs_col_names=["obs_col_1"],
            var_col_names=["var_col_1", "var_col_2"],
        ),
        keys_as_string=st.booleans(),
    )
    def test_from_var(self, data, adata: AnnData, keys_as_string: bool):
        adata.obs_names = [f"obs_{obs_id}" for obs_id in adata.obs_names]
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        if keys_as_string:
            keys = data.draw(st.sampled_from([*adata.var.columns]))
        else:
            keys = data.draw(
                st.lists(
                    st.sampled_from(adata.var.columns.to_list()),
                    min_size=1,
                    max_size=len(adata.var.columns),
                    unique=True,
                )
            )

        df = get_df(adata, keys=keys)

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == keys).all()
        if isinstance(keys, str):
            assert (df == adata.var[[keys]]).values.all()
        else:
            assert (df == adata.var[keys]).values.all()

    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["layer_1", "layer_2"],
            obsm_keys=["obsm_1"],
            varm_keys=["varm_1"],
            var_col_names=["var_col_1", "var_col_2"],
        ),
        layer_name=st.sampled_from(["layer_1", "layer_2"]),
    )
    def test_keys_as_layer(self, adata: AnnData, layer_name: str):
        adata.var_names = [f"var_{var_id}" for var_id in adata.var_names]

        df = get_df(adata, keys=layer_name)

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == adata.var_names).all()
        assert (df.values == adata.layers[layer_name]).all()

    @pytest.mark.parametrize("keys", (None, "col_1", "col_2", ["col_1", "col_2"]))
    def test_data_as_data_frame(self, keys):
        df = get_df(data=pd.DataFrame(np.eye(2), columns=["col_1", "col_2"]), keys=keys)

        assert isinstance(df, pd.DataFrame)
        if isinstance(keys, str):
            assert df.shape == (2, 1)
        else:
            assert df.shape == (2, 2)

    @pytest.mark.parametrize("dropna", (True, False, "any", "all"))
    def test_dropna(self, dropna: Union[bool, str]):
        df = get_df(
            data=pd.DataFrame(
                np.array([[1, np.nan, 0], [np.nan, 1, 0], [np.nan, np.nan, np.nan]]),
                columns=["col_1", "col_2", "col_3"],
            ),
            dropna=dropna,
        )

        assert isinstance(df, pd.DataFrame)
        assert (df.columns == ["col_1", "col_2", "col_3"]).all()
        if dropna == "all":
            np.testing.assert_equal(
                df.values, np.array([[1, np.nan, 0], [np.nan, 1, 0]])
            )
        elif dropna or dropna == "any":
            assert (
                df == pd.DataFrame(columns=["col_1", "col_2", "col_3"])
            ).values.all()
        else:
            np.testing.assert_equal(
                df.values,
                np.array([[1, np.nan, 0], [np.nan, 1, 0], [np.nan, np.nan, np.nan]]),
            )

    def test_index_from_obs_col(self):
        adata = AnnData(
            X=np.eye(2),
            layers={"layer_1": 2 * np.eye(2), "layer_2": 3 * np.eye(2)},
            obs=pd.DataFrame({"obs_col_1": ["a", "b"]}),
        )
        adata.var_names = ["var_name_1", "var_name_2"]

        df = get_df(adata, keys="layer_1", index="obs_col_1")

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(df.values, 2 * np.eye(2))
        assert (df.index == adata.obs["obs_col_1"]).all()
        assert (df.columns == ["var_name_1", "var_name_2"]).all()

    def test_columns_from_obs(self):
        adata = AnnData(
            X=np.eye(2),
            varm={"varm_1": 2 * np.eye(2)},
            obs=pd.DataFrame({"obs_col_1": ["a", "b"]}),
        )
        adata.var_names = ["var_name_1", "var_name_2"]

        df = get_df(adata, keys="varm_1", columns="obs_col_1")

        assert isinstance(df, pd.DataFrame)
        np.testing.assert_array_equal(df.values, 2 * np.eye(2))
        assert (df.index == ["var_name_1", "var_name_2"]).all()
        assert (df.columns == adata.obs["obs_col_1"]).all()

    @pytest.mark.parametrize("sparse", (True, False))
    @pytest.mark.parametrize("index", (None, ["index_1", "index_2"]))
    @pytest.mark.parametrize("columns", (None, ["col_1", "col_2"]))
    def test_data_as_array(
        self, index: Optional[List[str]], columns: Optional[List[str]], sparse: bool
    ):
        if sparse:
            data = csr_matrix(np.eye(2))
        else:
            data = np.eye(2)

        df = get_df(data, index=index, columns=columns)

        assert isinstance(df, pd.DataFrame)
        if index is None:
            assert (df.index == [0, 1]).all()
        else:
            assert (df.index == ["index_1", "index_2"]).all()
        if columns is None:
            assert (df.columns == [0, 1]).all()
        else:
            assert (df.columns == ["col_1", "col_2"]).all()


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


class TestMerge:
    def _assert_all_entries_present(
        self, returned_adata: AnnData, adata: AnnData, ldata: AnnData
    ):
        assert set(returned_adata.layers) == set(adata.layers).union(ldata.layers)
        assert set(returned_adata.obsm) == set(adata.obsm).union(ldata.obsm)
        assert set(returned_adata.varm) == set(adata.varm).union(ldata.varm)
        assert set(returned_adata.uns) == set(adata.uns).union(ldata.uns)
        assert set(returned_adata.obs.columns) == set(adata.obs.columns).union(
            ldata.obs.columns
        )
        assert set(returned_adata.var.columns) == set(adata.var.columns).union(
            ldata.var.columns
        )

    def _assert_copy_worked(self, returned_adata: AnnData, adata: AnnData, copy: bool):
        if copy:
            assert isinstance(returned_adata, AnnData)
        else:
            assert returned_adata is None
            returned_adata = adata

        return returned_adata

    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        ldata=get_adata(max_obs=5, max_vars=5),
        copy=st.booleans(),
    )
    def test_common_var_names(self, adata: AnnData, ldata: AnnData, copy: bool):
        adata.uns["a"] = ["a", 0, 3]
        ldata.uns["cluster_colors"] = {"cluster_1": "blue", "cluster_2": "red"}
        returned_adata = merge(adata=adata, ldata=ldata, copy=copy)

        returned_adata = self._assert_copy_worked(returned_adata, adata, copy)
        self._assert_all_entries_present(returned_adata, adata, ldata)

    @given(
        adata=get_adata(max_obs=5, max_vars=5, layer_keys=["spliced"]),
        ldata=get_adata(max_obs=5, max_vars=5),
        copy=st.booleans(),
    )
    def test_spliced_in_adata(self, adata: AnnData, ldata: AnnData, copy: bool):
        adata.uns["a"] = ["a", 0, 3]
        ldata.uns["cluster_colors"] = {"cluster_1": "blue", "cluster_2": "red"}
        returned_adata = merge(adata=adata, ldata=ldata, copy=copy)

        returned_adata = self._assert_copy_worked(returned_adata, adata, copy)
        self._assert_all_entries_present(returned_adata, adata, ldata)
        assert "initial_size" in returned_adata.obs.columns
        assert "initial_size_spliced" in returned_adata.obs.columns

    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        ldata=get_adata(max_obs=5, max_vars=5, layer_keys=["spliced"]),
        copy=st.booleans(),
    )
    def test_spliced_in_ldata(self, adata: AnnData, ldata: AnnData, copy: bool):
        adata.uns["a"] = ["a", 0, 3]
        ldata.uns["cluster_colors"] = {"cluster_1": "blue", "cluster_2": "red"}
        returned_adata = merge(adata=adata, ldata=ldata, copy=copy)

        returned_adata = self._assert_copy_worked(returned_adata, adata, copy)
        self._assert_all_entries_present(returned_adata, adata, ldata)
        assert "initial_size" in returned_adata.obs.columns
        assert "initial_size_spliced" in returned_adata.obs.columns

    @given(
        adata=get_adata(min_obs=3, max_obs=3, max_vars=5),
        ldata=get_adata(min_obs=3, max_obs=3, max_vars=5),
        copy=st.booleans(),
    )
    def test_no_common_obs_names(self, adata: AnnData, ldata: AnnData, copy: bool):
        adata.uns["a"] = ["a", 0, 3]
        ldata.uns["cluster_colors"] = {"cluster_1": "blue", "cluster_2": "red"}

        adata.obs_names = ["sample1_ABCD", "sample2_ABCD", "sample3_DCBA"]
        ldata.obs_names = ["_sample1_ABCD", "_sample2_ABCD", "_sample3_DCBA"]

        returned_adata = merge(adata=adata, ldata=ldata, copy=copy, id_length=4)

        returned_adata = self._assert_copy_worked(returned_adata, adata, copy)
        self._assert_all_entries_present(returned_adata, adata, ldata)
        assert returned_adata.obs_names.isin(["ABCD", "ABCD-1", "DCBA"]).all()


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

        if "X" in layers:
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
    @given(data=st.data(), adata=get_adata(max_obs=5, max_vars=5))
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
