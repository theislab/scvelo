from typing import List, Optional, Union

from numpy import ndarray
from scipy.sparse import csr_matrix, issparse, spmatrix

from pandas import DataFrame

from anndata import AnnData
import scvelo.logging as logg


def get_modality(adata: AnnData, modality: str) -> Union[ndarray, spmatrix, DataFrame]:
    """Extract data of one modality.

    Args:
        adata (AnnData): Annotated data to extract modality from.
        modality (str): Modality for which data is needed.
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


# TODO: Add unit test
def make_dense(
    adata: AnnData, modalities: Union[List[str], str], inplace: bool = True
) -> Optional[AnnData]:
    """Densify sparse AnnData entry

    Args:
        adata (AnnData): Annotated data object.
        modality (str): Modality to make dense.
        inplace (bool): Boolean flag to perform operations inplace or not. Defaults to
            `True`.

    Returns:
        Copy of annotated data `adata` if `inplace=True`.
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
    """Make AnnData entry sparse

    Args:
        adata (AnnData): Annotated data object.
        modality (str): Modality to make sparse.
        inplace (bool): Boolean flag to perform operations inplace or not. Defaults to
            `True`.

    Returns:
        Copy of annotated data `adata` if `inplace=True`.
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

    Args:
        adata (AnnData): Annotated data object.
        new_value (ndarray, spmatrix or DataFrame): New value of modality.
        modality (str or None): Modality to overwrite with new value. Defaults to
            `None`.
        inplace (bool): Boolean flag to indicate whether setting of modality should be
            inplace or not. Defaults to `True`.

    Returns:

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
