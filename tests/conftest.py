from datetime import timedelta

import pytest
from hypothesis import settings

import scanpy as sc
from anndata import AnnData

settings.register_profile("ci", deadline=timedelta(milliseconds=500))


# TODO: Make datasets smaller (less variables)
_dentategyrus_50obs = sc.read("tests/_data/dentategyrus_50obs.h5ad")
_dentategyrus_50obs_preprocessed = sc.read(
    "tests/_data/dentategyrus_50obs_preprocessed.h5ad"
)
_dentategyrus_100obs = sc.read("tests/_data/dentategyrus_100obs.h5ad")
_dentategyrus_100obs_preprocessed = sc.read(
    "tests/_data/dentategyrus_100obs_preprocessed.h5ad"
)
_pancreas_50obs = sc.read("tests/_data/pancreas_50obs.h5ad")
_pancreas_50obs_preprocessed = sc.read("tests/_data/pancreas_50obs_preprocessed.h5ad")
_pancreas_100obs = sc.read("tests/_data/pancreas_100obs.h5ad")
_pancreas_100obs_preprocessed = sc.read("tests/_data/pancreas_100obs_preprocessed.h5ad")


@pytest.fixture
def dentategyrus_50obs() -> AnnData:
    return _dentategyrus_50obs.copy()


@pytest.fixture
def dentategyrus_50obs_preprocessed() -> AnnData:
    """Preprocessed dentategyrus dataset with 50 observations.

    The data has been preprocessed using

    .. code:: python
        import scanpy as sc
        import scvelo as scv

        adata = sc.read(f"tests/_data/dentategyrus_50obs.h5ad")

        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=200,
            retain_genes=None,
            subset_highly_variable=True,
            flavor="seurat",
            log=True,
            layers_normalize=None,
            copy=False,
        )
        scv.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=None,
            use_rep=None,
            use_highly_variable=True,
            knn=True,
            random_state=0,
            method="umap",
            metric="euclidean",
            metric_kwds=None,
            num_threads=-1,
            copy=False,
        )
        scv.pp.moments(
            adata,
            n_neighbors=30,
            n_pcs=None,
            mode="connectivities",
            method="umap",
            use_rep=None,
            use_highly_variable=True,
            copy=False,
        )
    """

    return _dentategyrus_50obs_preprocessed.copy()


@pytest.fixture
def dentategyrus_100obs() -> AnnData:
    return _dentategyrus_100obs.copy()


@pytest.fixture
def dentategyrus_100obs_preprocessed() -> AnnData:
    """Preprocessed dentategyrus dataset with 100 observations.

    The data has been preprocessed using

    .. code:: python
        import scanpy as sc
        import scvelo as scv

        adata = sc.read(f"tests/_data/dentategyrus_100obs.h5ad")

        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=200,
            retain_genes=None,
            subset_highly_variable=True,
            flavor="seurat",
            log=True,
            layers_normalize=None,
            copy=False,
        )
        scv.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=None,
            use_rep=None,
            use_highly_variable=True,
            knn=True,
            random_state=0,
            method="umap",
            metric="euclidean",
            metric_kwds=None,
            num_threads=-1,
            copy=False,
        )
        scv.pp.moments(
            adata,
            n_neighbors=30,
            n_pcs=None,
            mode="connectivities",
            method="umap",
            use_rep=None,
            use_highly_variable=True,
            copy=False,
        )
    """

    return _dentategyrus_100obs_preprocessed.copy()


@pytest.fixture
def pancreas_50obs() -> AnnData:
    return _pancreas_50obs.copy()


@pytest.fixture
def pancreas_50obs_preprocessed() -> AnnData:
    """Preprocessed pancreas dataset with 50 observations.

    The data has been preprocessed using

    .. code:: python
        import scanpy as sc
        import scvelo as scv

        adata = sc.read(f"tests/_data/pancreas_50obs.h5ad")

        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=200,
            retain_genes=None,
            subset_highly_variable=True,
            flavor="seurat",
            log=True,
            layers_normalize=None,
            copy=False,
        )
        scv.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=None,
            use_rep=None,
            use_highly_variable=True,
            knn=True,
            random_state=0,
            method="umap",
            metric="euclidean",
            metric_kwds=None,
            num_threads=-1,
            copy=False,
        )
        scv.pp.moments(
            adata,
            n_neighbors=30,
            n_pcs=None,
            mode="connectivities",
            method="umap",
            use_rep=None,
            use_highly_variable=True,
            copy=False,
        )
    """

    return _pancreas_50obs_preprocessed.copy()


@pytest.fixture
def pancreas_100obs() -> AnnData:
    return _pancreas_100obs.copy()


@pytest.fixture
def pancreas_100obs_preprocessed() -> AnnData:
    """Preprocessed dentategyrus dataset with 100 observations.

    The data has been preprocessed using

    .. code:: python
        import scanpy as sc
        import scvelo as scv

        adata = sc.read(f"tests/_data/pancreas_100obs.h5ad")

        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=200,
            retain_genes=None,
            subset_highly_variable=True,
            flavor="seurat",
            log=True,
            layers_normalize=None,
            copy=False,
        )
        scv.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=None,
            use_rep=None,
            use_highly_variable=True,
            knn=True,
            random_state=0,
            method="umap",
            metric="euclidean",
            metric_kwds=None,
            num_threads=-1,
            copy=False,
        )
        scv.pp.moments(
            adata,
            n_neighbors=30,
            n_pcs=None,
            mode="connectivities",
            method="umap",
            use_rep=None,
            use_highly_variable=True,
            copy=False,
        )
    """

    return _pancreas_100obs_preprocessed.copy()
