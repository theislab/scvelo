import hypothesis.strategies as st
from hypothesis import given

import numpy as np

from anndata import AnnData

from scvelo.core import get_modality
from scvelo.core.tests import get_adata, TestBase
from scvelo.preprocessing._filter import CellCountFilter, GeneCountFilter


# TODO: Generalize test to test keyword `shared_counts`
# TODO: Decouple function `_get_n_vars_final` as it is the same for
# `TestGeneCountFilter`
class TestCellCountFilter(TestBase):
    def _get_n_vars_final(self, adata, modalities_to_filter, min_cells, max_cells):
        n_vars_final = {
            key: ((val.sum(axis=0) >= min_cells) & (val.sum(axis=0) <= max_cells)).sum()
            for key, val in modalities_to_filter.items()
            if key in adata.obsm
        }

        x_or_layers = {
            key: val
            for key, val in modalities_to_filter.items()
            if (key == "X" or key in adata.layers)
        }
        final_counts = np.sum(
            list(
                map(
                    np.all,
                    zip(
                        *[
                            (val.sum(axis=0) >= min_cells)
                            & (val.sum(axis=0) <= max_cells)
                            for key, val in x_or_layers.items()
                        ]
                    ),
                )
            )
        )
        n_vars_final.update({key: final_counts for key in x_or_layers})

        return n_vars_final

    @given(
        adata=get_adata(),
        inplace=st.booleans(),
        enforce=st.booleans(),
        min_cells=st.integers(min_value=-1000, max_value=1000),
        max_cells=st.integers(min_value=-1000, max_value=1000),
        n_modalities=st.integers(min_value=0),
    )
    def test_cell_count_filter(
        self,
        adata: AnnData,
        inplace: bool,
        enforce: bool,
        min_cells: int,
        max_cells: int,
        n_modalities: int,
    ):
        n_obs, n_vars = adata.shape

        modalities_to_filter = self._subset_modalities(adata, n_modalities)
        modalities_to_filter = {
            modality: get_modality(adata=adata, modality=modality) > 0
            for modality in modalities_to_filter
        }

        n_vars_final = self._get_n_vars_final(
            adata, modalities_to_filter, min_cells, max_cells
        )

        cell_count_filter = CellCountFilter(
            inplace=inplace, enforce=enforce, min_cells=min_cells, max_cells=max_cells
        )

        if inplace:
            cell_count_filter.transform(adata, modalities=modalities_to_filter)
        else:
            adata = cell_count_filter.transform(adata, modalities=modalities_to_filter)

        assert adata.n_obs == n_obs
        assert set(
            [
                "initial_size" if modality == "X" else f"initial_size_{modality}"
                for modality in modalities_to_filter
            ]
        ).issubset(adata.obs.columns)
        for modality in modalities_to_filter:
            assert get_modality(adata, modality).shape[1] == n_vars_final[modality]


# TODO: Generalize test to test keyword `shared_counts`
# TODO: Decouple function `_get_n_vars_final` as it is the same for
# `TestCellCountFilter`
class TestGeneCountFilter(TestBase):
    def _get_n_vars_final(self, adata, modalities_to_filter, min_counts, max_counts):
        n_vars_final = {
            key: (
                (val.sum(axis=0) >= min_counts) & (val.sum(axis=0) <= max_counts)
            ).sum()
            for key, val in modalities_to_filter.items()
            if key in adata.obsm
        }

        x_or_layers = {
            key: val
            for key, val in modalities_to_filter.items()
            if (key == "X" or key in adata.layers)
        }
        final_counts = np.sum(
            list(
                map(
                    np.all,
                    zip(
                        *[
                            (val.sum(axis=0) >= min_counts)
                            & (val.sum(axis=0) <= max_counts)
                            for key, val in x_or_layers.items()
                        ]
                    ),
                )
            )
        )
        n_vars_final.update({key: final_counts for key in x_or_layers})

        return n_vars_final

    @given(
        adata=get_adata(),
        inplace=st.booleans(),
        enforce=st.booleans(),
        min_counts=st.integers(min_value=-1000, max_value=1000),
        max_counts=st.integers(min_value=-1000, max_value=1000),
        n_modalities=st.integers(min_value=0),
    )
    def test_gene_count_filter(
        self,
        adata: AnnData,
        inplace: bool,
        enforce: bool,
        min_counts: int,
        max_counts: int,
        n_modalities: int,
    ):
        n_obs, n_vars = adata.shape

        modalities_to_filter = self._subset_modalities(adata, n_modalities)
        modalities_to_filter = {
            modality: get_modality(adata=adata, modality=modality)
            for modality in modalities_to_filter
        }

        n_vars_final = self._get_n_vars_final(
            adata, modalities_to_filter, min_counts, max_counts
        )

        gene_count_filter = GeneCountFilter(
            inplace=inplace,
            enforce=enforce,
            min_counts=min_counts,
            max_counts=max_counts,
        )
        if inplace:
            gene_count_filter.transform(adata, modalities=modalities_to_filter)
        else:
            adata = gene_count_filter.transform(adata, modalities=modalities_to_filter)

        assert adata.n_obs == n_obs
        assert set(
            [
                "initial_size" if modality == "X" else f"initial_size_{modality}"
                for modality in modalities_to_filter
            ]
        ).issubset(adata.obs.columns)
        for modality in modalities_to_filter:
            assert get_modality(adata, modality).shape[1] == n_vars_final[modality]
