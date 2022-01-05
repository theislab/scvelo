from datetime import timedelta

import pytest
from hypothesis import settings

import scanpy as sc
from anndata import AnnData

settings.register_profile("ci", deadline=timedelta(milliseconds=500))


# TODO: Make datasets smaller (less variables)
_dentategyrus_50obs = sc.read("tests/_data/dentategyrus_50obs.h5ad")
_dentategyrus_100obs = sc.read("tests/_data/dentategyrus_100obs.h5ad")
_pancreas_50obs = sc.read("tests/_data/pancreas_50obs.h5ad")
_pancreas_100obs = sc.read("tests/_data/pancreas_100obs.h5ad")


@pytest.fixture
def dentategyrus_50obs() -> AnnData:
    return _dentategyrus_50obs.copy()


@pytest.fixture
def dentategyrus_100obs() -> AnnData:
    return _dentategyrus_100obs.copy()


@pytest.fixture
def pancreas_50obs() -> AnnData:
    return _pancreas_50obs.copy()


@pytest.fixture
def pancreas_100obs() -> AnnData:
    return _pancreas_100obs.copy()
