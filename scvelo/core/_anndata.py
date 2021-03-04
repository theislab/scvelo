from typing import List, Optional, Union

from numpy import ndarray
from pandas import DataFrame
from scipy.sparse import csr_matrix, issparse, spmatrix

from anndata import AnnData

import scvelo.logging as logg


def get_modality(adata: AnnData, modality: str) -> Union[ndarray, spmatrix]:
    """Extract data of one modality.

    Arguments
    ---------
    adata:
        Annotated data to extract modality from.
    modality:
        Modality for which data is needed.

    Returns
    -------
    Union[ndarray, spmatrix]
        Retrieved modality from :class:`~anndata.AnnData` object.

    """

    if modality == "X":
        return adata.X
    elif modality in adata.layers.keys():
        return adata.layers[modality]
    elif modality in adata.obsm.keys():
        if isinstance(adata.obsm[modality], DataFrame):
            return adata.obsm[modality].values
        else:
            return adata.obsm[modality]


def make_dense(
    adata: AnnData, modalities: Union[List[str], str], inplace: bool = True
) -> Optional[AnnData]:
    """Densify sparse AnnData entry.

    Arguments
    ---------
    adata:
        Annotated data object.
    modality:
        Modality to make dense.
    inplace:
        Boolean flag to perform operations inplace or not. Defaults to `True`.

    Returns
    -------
    Optional[AnnData]
        Copy of annotated data `adata` if `inplace=True` with dense modalities.

    """

    if not inplace:
        adata = adata.copy()

    if isinstance(modalities, str):
        modalities = [modalities]

    # Densify modalities
    for modality in modalities:
        count_data = get_modality(adata=adata, modality=modality)
        if issparse(count_data):
            set_modality(adata=adata, modality=modality, new_value=count_data.A)

    return adata if not inplace else None


# TODO: Allow choosing format of sparse matrix i.e., csr, csc, ...
def make_sparse(
    adata: AnnData, modalities: Union[List[str], str], inplace: bool = True
) -> Optional[AnnData]:
    """Make AnnData entry sparse.

    Arguments
    ---------
    adata:
        Annotated data object.
    modality:
        Modality to make sparse.
    inplace:
        Boolean flag to perform operations inplace or not. Defaults to `True`.

    Returns
    -------
    Optional[AnnData]
        Copy of annotated data `adata` with sparse modalities if `inplace=True`.

    """

    if not inplace:
        adata = adata.copy()

    if isinstance(modalities, str):
        modalities = [modalities]

    # Make modalities sparse
    for modality in modalities:
        count_data = get_modality(adata=adata, modality=modality)
        if modality == "X":
            logg.warn("Making `X` sparse is not supported.")
        elif not issparse(count_data):
            set_modality(
                adata=adata, modality=modality, new_value=csr_matrix(count_data)
            )

    return adata if not inplace else None


def set_modality(
    adata: AnnData,
    new_value: Union[ndarray, spmatrix, DataFrame],
    modality: Optional[str] = None,
    inplace: bool = True,
) -> Optional[AnnData]:
    """Set modality of annotated data object to new value.

    Arguments
    ---------
    adata:
        Annotated data object.
    new_value:
        New value of modality.
    modality:
        Modality to overwrite with new value. Defaults to `None`.
    inplace:
        Boolean flag to indicate whether setting of modality should be inplace or
            not. Defaults to `True`.

    Returns
    -------
    Optional[AnnData]
        Copy of annotated data `adata` with updated modality if `inplace=True`.

    """

    if not inplace:
        adata = adata.copy()

    if (modality == "X") or (modality is None):
        adata.X = new_value
    elif modality in adata.layers.keys():
        adata.layers[modality] = new_value
    elif modality in adata.obsm.keys():
        adata.obsm[modality] = new_value

    if not inplace:
        return adata
