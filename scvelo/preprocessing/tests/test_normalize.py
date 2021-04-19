from hypothesis import given
from hypothesis import strategies as st

from anndata import AnnData

from scvelo.core.tests import get_adata, TestBase
from scvelo.preprocessing._normalize import CountNormalizer


# TODO: Make test more sophisticated
class TestCountNormalizer(TestBase):
    @given(
        adata=get_adata(),
        inplace=st.booleans(),
        enforce=st.booleans(),
        n_modalities=st.integers(min_value=0),
    )
    def test_count_normalizer(
        self, adata: AnnData, inplace: bool, enforce: bool, n_modalities: int
    ):
        self._convert_to_float(adata=adata)
        n_obs, n_vars = adata.shape

        modalities_to_filter = self._subset_modalities(adata, n_modalities)

        count_normalizer = CountNormalizer(inplace=inplace, enforce=enforce)

        assert not any(
            [
                count_normalizer._is_normalized(adata, modality)
                for modality in modalities_to_filter
            ]
        )

        if inplace:
            count_normalizer.transform(adata, modalities=modalities_to_filter)
        else:
            adata = count_normalizer.transform(adata, modalities=modalities_to_filter)

        assert adata.shape == (n_obs, n_vars)
