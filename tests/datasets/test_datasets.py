import os
import sys
from typing import Optional

import pytest

from anndata import AnnData

import scvelo as scv


@pytest.fixture(scope="session")
def dentategyrus_adata(tmpdir_factory):
    path_to_file = tmpdir_factory.mktemp("dentategyrus").join("adata.h5ad")
    _ = scv.datasets.dentategyrus(file_path=path_to_file)

    return path_to_file


@pytest.mark.skipif(
    sys.version_info[:2] != (3, 10) or sys.platform != "linux",
    reason="Limit number of downloads to speed up testing.",
)
class TestDataSets:
    def test_dentategyrus_adjusted(self, dentategyrus_adata):
        adata = scv.datasets.dentategyrus(file_path=dentategyrus_adata, adjusted=True)

        assert isinstance(adata, AnnData)
        assert adata.shape == (2930, 13913)

    def test_dentategyrus_not_adjusted(self, tmpdir_factory):
        adata = scv.datasets.dentategyrus(
            file_path=tmpdir_factory.mktemp("dentategyrus").join("loomfile.loom"),
            adjusted=False,
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (3396, 25919)

    @pytest.mark.parametrize("n_obs", (None, 0, 1, 10))
    def test_toy_data(self, dentategyrus_adata, n_obs: Optional[int]):
        assert os.path.isfile(dentategyrus_adata)
        adata = scv.datasets.toy_data(file_path=dentategyrus_adata, n_obs=n_obs)

        assert isinstance(adata, AnnData)
        if n_obs is None:
            assert adata.shape == (2930, 13913)
        else:
            assert adata.shape == (n_obs, 13913)

    def test_forebrain(self, tmpdir_factory):
        adata = scv.datasets.forebrain(
            file_path=tmpdir_factory.mktemp("forebrain").join("loomfile.loom")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (1720, 32738)

    def test_pancreas(self, tmpdir_factory):
        adata = scv.datasets.pancreas(
            file_path=tmpdir_factory.mktemp("pancreas").join("adata.h5ad")
        )

        assert isinstance(adata, AnnData)
        assert adata.shape == (3696, 27998)
