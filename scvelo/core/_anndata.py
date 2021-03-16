import re
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.sparse import csr_matrix, issparse, spmatrix

from anndata import AnnData

import scvelo.logging as logg
from ._arithmetic import sum


def clean_obs_names(data, base="[AGTCBDHKMNRSVWY]", ID_length=12, copy=False):
    """Clean up the obs_names.

    For example an obs_name 'sample1_AGTCdate' is changed to 'AGTC' of the sample
    'sample1_date'. The sample name is then saved in obs['sample_batch'].
    The genetic codes are identified according to according to
    https://www.neb.com/tools-and-resources/usage-guidelines/the-genetic-code.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    base: `str` (default: `[AGTCBDHKMNRSVWY]`)
        Genetic code letters to be identified.
    ID_length: `int` (default: 12)
        Length of the Genetic Codes in the samples.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    obs_names: list
        updated names of the observations
    sample_batch: `.obs`
        names of the identified sample batches
    """

    def get_base_list(name, base):
        base_list = base
        while re.search(base_list + base, name) is not None:
            base_list += base
        if len(base_list) == 0:
            raise ValueError("Encountered an invalid ID in obs_names: ", name)
        return base_list

    adata = data.copy() if copy else data

    names = adata.obs_names
    base_list = get_base_list(names[0], base)

    if len(np.unique([len(name) for name in adata.obs_names])) == 1:
        start, end = re.search(base_list, names[0]).span()
        newIDs = [name[start:end] for name in names]
        start, end = 0, len(newIDs[0])
        for i in range(end - ID_length):
            if np.any([ID[i] not in base for ID in newIDs]):
                start += 1
            if np.any([ID[::-1][i] not in base for ID in newIDs]):
                end -= 1

        newIDs = [ID[start:end] for ID in newIDs]
        prefixes = [names[i].replace(newIDs[i], "") for i in range(len(names))]
    else:
        prefixes, newIDs = [], []
        for name in names:
            match = re.search(base_list, name)
            newID = (
                re.search(get_base_list(name, base), name).group()
                if match is None
                else match.group()
            )
            newIDs.append(newID)
            prefixes.append(name.replace(newID, ""))

    adata.obs_names = newIDs
    if len(prefixes[0]) > 0 and len(np.unique(prefixes)) > 1:
        adata.obs["sample_batch"] = (
            pd.Categorical(prefixes)
            if len(np.unique(prefixes)) < adata.n_obs
            else prefixes
        )

    adata.obs_names_make_unique()
    return adata if copy else None


def get_initial_size(adata, layer=None, by_total_size=None):
    if by_total_size:
        sizes = [
            adata.obs[f"initial_size_{layer}"]
            for layer in {"spliced", "unspliced"}
            if f"initial_size_{layer}" in adata.obs.keys()
        ]
        return np.sum(sizes, axis=0)
    elif layer in adata.layers.keys():
        return (
            np.array(adata.obs[f"initial_size_{layer}"])
            if f"initial_size_{layer}" in adata.obs.keys()
            else get_size(adata, layer)
        )
    elif layer is None or layer == "X":
        return (
            np.array(adata.obs["initial_size"])
            if "initial_size" in adata.obs.keys()
            else get_size(adata)
        )
    else:
        return None


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


def get_size(adata, layer=None):
    X = adata.X if layer is None else adata.layers[layer]
    return sum(X, axis=1)


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


def set_initial_size(adata, layers=None):
    if layers is None:
        layers = ["spliced", "unspliced"]
    verify_dtypes(adata)
    layers = [
        layer
        for layer in layers
        if layer in adata.layers.keys()
        and f"initial_size_{layer}" not in adata.obs.keys()
    ]
    for layer in layers:
        adata.obs[f"initial_size_{layer}"] = get_size(adata, layer)
    if "initial_size" not in adata.obs.keys():
        adata.obs["initial_size"] = get_size(adata)


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


def verify_dtypes(adata):
    try:
        _ = adata[:, 0]
    except Exception:
        uns = adata.uns
        adata.uns = {}
        try:
            _ = adata[:, 0]
            logg.warn(
                "Safely deleted unstructured annotations (adata.uns), \n"
                "as these do not comply with permissible anndata datatypes."
            )
        except Exception:
            logg.warn(
                "The data might be corrupted. Please verify all annotation datatypes."
            )
            adata.uns = uns
