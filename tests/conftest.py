from datetime import timedelta
from typing import Tuple, Union

import pytest
from hypothesis import settings

import matplotlib as mpl

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


def pytest_sessionstart(session: pytest.Session) -> None:
    del session
    mpl.use("Agg")


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


def _get_dentategyrus_50obs(
    raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of dentategyrus dataset with 50 observations.

    Parameters
    ----------
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if raw and not preprocessed:
        return _dentategyrus_50obs.copy()
    elif not raw and preprocessed:
        return _dentategyrus_50obs_preprocessed.copy()
    elif raw and preprocessed:
        return _dentategyrus_50obs.copy(), _dentategyrus_50obs_preprocessed.copy()


def _get_dentategyrus_100obs(
    raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of dentategyrus dataset with 100 observations.

    Parameters
    ----------
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if raw and not preprocessed:
        return _dentategyrus_100obs.copy()
    elif not raw and preprocessed:
        return _dentategyrus_100obs_preprocessed.copy()
    elif raw and preprocessed:
        return _dentategyrus_100obs.copy(), _dentategyrus_100obs_preprocessed.copy()


def _get_dentategyrus_adata(
    n_obs: int, raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of raw or preprocessed dentategyrus dataset.

    Parameters
    ----------
    n_obs
        Number of observations of dataset to return.
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if n_obs == 50:
        return _get_dentategyrus_50obs(raw=raw, preprocessed=preprocessed)
    elif n_obs == 100:
        return _get_dentategyrus_100obs(raw=raw, preprocessed=preprocessed)


def _get_pancreas_50obs(
    raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of pancreas dataset with 50 observations.

    Parameters
    ----------
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if raw and not preprocessed:
        return _pancreas_50obs.copy()
    elif not raw and preprocessed:
        return _pancreas_50obs_preprocessed.copy()
    elif raw and preprocessed:
        return _pancreas_50obs.copy(), _pancreas_50obs_preprocessed.copy()


def _get_pancreas_100obs(
    raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of raw or preprocessed pancreas dataset with 100 observations.

    Parameters
    ----------
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if raw and not preprocessed:
        return _pancreas_100obs.copy()
    elif not raw and preprocessed:
        return _pancreas_100obs_preprocessed.copy()
    elif raw and preprocessed:
        return _pancreas_100obs.copy(), _pancreas_100obs_preprocessed.copy()


def _get_pancreas_adata(
    n_obs: int, raw: bool, preprocessed: bool
) -> Union[AnnData, Tuple[AnnData, AnnData]]:
    """Get AnnData of raw or preprocessed pancreas dataset.

    Parameters
    ----------
    n_obs
        Number of observations of dataset to return.
    raw
        Boolean identifier whether or not to return raw dataset.
    preprocessed
        Boolean identifier whether or not to return preprocessed dataset.

    Returns:
    -------
    Union[AnnData, Tuple[AnnData, AnnData]]
        Specified version of dataset.
    """

    if n_obs == 50:
        return _get_pancreas_50obs(raw=raw, preprocessed=preprocessed)
    elif n_obs == 100:
        return _get_pancreas_100obs(raw=raw, preprocessed=preprocessed)


@pytest.fixture
def adata():
    """Fixture to easily use available datasets in unit tests.

    The fixture returns a function to load the AnnData objects of a specified dataset
    (`"pancreas"` or `"dentategyrus"`). The function is then used in the unit test to
    load the needed version(s) (raw or preprocessed) of the dataset.
    """

    def _get_adata(dataset: str, n_obs: int, raw: bool, preprocessed: bool):
        if dataset == "pancreas":
            return _get_pancreas_adata(n_obs=n_obs, raw=raw, preprocessed=preprocessed)
        elif dataset == "dentategyrus":
            return _get_dentategyrus_adata(
                n_obs=n_obs, raw=raw, preprocessed=preprocessed
            )

    return _get_adata
