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


def cleanup(data, clean="layers", keep=None, copy=False):
    """Delete not needed attributes.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    clean: `str` or list of `str` (default: `layers`)
        Which attributes to consider for freeing memory.
    keep: `str` or list of `str` (default: None)
        Which attributes to keep.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with selection of attributes kept.
    """
    adata = data.copy() if copy else data
    verify_dtypes(adata)

    keep = list([keep] if isinstance(keep, str) else {} if keep is None else keep)
    keep.extend(["spliced", "unspliced", "Ms", "Mu", "clusters", "neighbors"])

    ann_dict = {
        "obs": adata.obs_keys(),
        "var": adata.var_keys(),
        "uns": adata.uns_keys(),
        "layers": list(adata.layers.keys()),
    }

    if "all" not in clean:
        ann_dict = {ann: values for (ann, values) in ann_dict.items() if ann in clean}

    for (ann, values) in ann_dict.items():
        for value in values:
            if value not in keep:
                del getattr(adata, ann)[value]

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


def merge(adata, ldata, copy=True):
    """Merge two annotated data matrices.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix (reference data set).
    ldata: :class:`~anndata.AnnData`
        Annotated data matrix (to be merged into adata).

    Returns
    -------
    Returns a :class:`~anndata.AnnData` object
    """
    adata.var_names_make_unique()
    ldata.var_names_make_unique()

    if (
        "spliced" in ldata.layers.keys()
        and "initial_size_spliced" not in ldata.obs.keys()
    ):
        set_initial_size(ldata)
    elif (
        "spliced" in adata.layers.keys()
        and "initial_size_spliced" not in adata.obs.keys()
    ):
        set_initial_size(adata)

    common_obs = pd.unique(adata.obs_names.intersection(ldata.obs_names))
    common_vars = pd.unique(adata.var_names.intersection(ldata.var_names))

    if len(common_obs) == 0:
        clean_obs_names(adata)
        clean_obs_names(ldata)
        common_obs = adata.obs_names.intersection(ldata.obs_names)

    if copy:
        _adata = adata[common_obs].copy()
        _ldata = ldata[common_obs].copy()
    else:
        adata._inplace_subset_obs(common_obs)
        _adata, _ldata = adata, ldata[common_obs].copy()

    _adata.var_names_make_unique()
    _ldata.var_names_make_unique()

    same_vars = len(_adata.var_names) == len(_ldata.var_names) and np.all(
        _adata.var_names == _ldata.var_names
    )
    join_vars = len(common_vars) > 0

    if join_vars and not same_vars:
        _adata._inplace_subset_var(common_vars)
        _ldata._inplace_subset_var(common_vars)

    for attr in _ldata.obs.keys():
        if attr not in _adata.obs.keys():
            _adata.obs[attr] = _ldata.obs[attr]
    for attr in _ldata.obsm.keys():
        if attr not in _adata.obsm.keys():
            _adata.obsm[attr] = _ldata.obsm[attr]
    for attr in _ldata.uns.keys():
        if attr not in _adata.uns.keys():
            _adata.uns[attr] = _ldata.uns[attr]
    if join_vars:
        for attr in _ldata.layers.keys():
            if attr not in _adata.layers.keys():
                _adata.layers[attr] = _ldata.layers[attr]

        if _adata.shape[1] == _ldata.shape[1]:
            same_vars = len(_adata.var_names) == len(_ldata.var_names) and np.all(
                _adata.var_names == _ldata.var_names
            )
            if same_vars:
                for attr in _ldata.var.keys():
                    if attr not in _adata.var.keys():
                        _adata.var[attr] = _ldata.var[attr]
                for attr in _ldata.varm.keys():
                    if attr not in _adata.varm.keys():
                        _adata.varm[attr] = _ldata.varm[attr]
            else:
                raise ValueError("Variable names are not identical.")

    return _adata if copy else None


def obs_df(adata, keys, layer=None):
    lookup_keys = [k for k in keys if k in adata.var_names]
    if len(lookup_keys) < len(keys):
        logg.warn(
            f"Keys {[k for k in keys if k not in adata.var_names]} "
            f"were not found in `adata.var_names`."
        )

    df = pd.DataFrame(index=adata.obs_names)
    for lookup_key in lookup_keys:
        df[lookup_key] = adata.obs_vector(lookup_key, layer=layer)
    return df


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


def var_df(adata, keys, layer=None):
    lookup_keys = [k for k in keys if k in adata.obs_names]
    if len(lookup_keys) < len(keys):
        logg.warn(
            f"Keys {[k for k in keys if k not in adata.obs_names]} "
            f"were not found in `adata.obs_names`."
        )

    df = pd.DataFrame(index=adata.var_names)
    for lookup_key in lookup_keys:
        df[lookup_key] = adata.var_vector(lookup_key, layer=layer)
    return df


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
