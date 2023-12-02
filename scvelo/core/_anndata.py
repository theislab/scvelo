import re
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from pandas.api.types import is_categorical_dtype
from scipy.sparse import csr_matrix, issparse, spmatrix

from anndata import AnnData

from scvelo import logging as logg
from ._arithmetic import sum
from ._utils import deprecated_arg_names


@deprecated_arg_names(
    {"data": "adata", "copy": "inplace", "ID_length": "id_length", "base": "alphabet"}
)
def clean_obs_names(
    adata: AnnData,
    alphabet: str = "[AGTCBDHKMNRSVWY]",
    id_length: int = 12,
    inplace: bool = True,
) -> Optional[AnnData]:
    """Clean up the obs_names.

    For example an obs_name 'sample1_AGTCdate' is changed to 'AGTC' of the sample 'sample1_date'. The sample name is
    then saved in obs['sample_batch']. The genetic codes are identified according to according to
    https://www.neb.com/tools-and-resources/usage-guidelines/the-genetic-code.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    alphabet
        Genetic code letters to be identified.
    id_length
        Length of the Genetic Codes in the samples.
    inplace
        Whether to update `adata` inplace or not.

    Returns
    -------
    Optional[AnnData]
        Returns or updates `adata` with updated names of observations, and names of identified sample batches added
        as column `"sample_batch"` to `.obs`.
    """
    if not inplace:
        adata = adata.copy()

    if adata.obs_names.map(len).unique().size == 1:
        start, end = re.search(alphabet * id_length, adata.obs_names[0]).span()
        new_obs_names = [obs_name[start:end] for obs_name in adata.obs_names]

        prefixes = [
            obs_name.replace(new_obs_name, "")
            for obs_name, new_obs_name in zip(adata.obs_names, new_obs_names)
        ]
    else:
        prefixes, new_obs_names = [], []
        for obs_name in adata.obs_names:
            start, end = re.search(alphabet * id_length, adata.obs_names[0]).span()
            new_obs_names.append(obs_name[start:end])
            prefixes.append(obs_name.replace(obs_name[start:end], ""))

    adata.obs_names = new_obs_names
    adata.obs_names_make_unique()

    if len(prefixes[0]) > 0 and len(np.unique(prefixes)) > 1:
        adata.obs["sample_batch"] = (
            pd.Categorical(prefixes)
            if len(np.unique(prefixes)) < adata.n_obs
            else prefixes
        )

    if not inplace:
        return adata


@deprecated_arg_names({"data": "adata", "copy": "inplace"})
def cleanup(
    adata: AnnData,
    clean: Union[
        Literal["layers", "obs", "var", "uns", "all"],
        List[Literal["layers", "obs", "var", "uns"]],
    ] = "layers",
    keep: Optional[Union[str, List[str]]] = None,
    inplace: bool = True,
) -> Optional[AnnData]:
    """Delete not needed attributes.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    clean
        Which attributes to consider for freeing memory.
    keep
        Which attributes to keep.
    inplace
        Whether to update `adata` inplace or not.

    Returns
    -------
    Optional[AnnData]
        Returns or updates `adata` with selection of attributes kept.
    """
    if not inplace:
        adata = adata.copy()
    verify_dtypes(adata)

    keep = list([keep] if isinstance(keep, str) else {} if keep is None else keep)
    keep.extend(["unspliced", "spliced", "Mu", "Ms", "clusters", "neighbors"])

    attributes_to_remove = {
        "obs": adata.obs_keys(),
        "var": adata.var_keys(),
        "uns": adata.uns_keys(),
        "layers": list(adata.layers.keys()),
    }

    if "all" not in clean:
        attributes_to_remove = {
            attr: attr_keys
            for (attr, attr_keys) in attributes_to_remove.items()
            if attr in clean
        }

    for attr, attr_keys in attributes_to_remove.items():
        for key in attr_keys:
            if key not in keep:
                del getattr(adata, attr)[key]

    if not inplace:
        return adata


# TODO: Add unit test for `precision` argument
def get_df(
    data: AnnData,
    keys: Optional[Union[str, List[str]]] = None,
    layer: Optional[str] = None,
    index: List = None,
    columns: List = None,
    sort_values: bool = None,
    dropna: Literal["all", "any"] = "all",
    precision: int = None,
) -> DataFrame:
    """Get dataframe for a specified adata key.

    Return values for specified key
    (in obs, var, obsm, varm, obsp, varp, uns, or layers) as a dataframe.

    Arguments:
    ---------
    data
        AnnData object or a numpy array to get values from.
    keys
        Keys from `.var_names`, `.obs_names`, `.var`, `.obs`,
        `.obsm`, `.varm`, `.obsp`, `.varp`, `.uns`, or `.layers`.
    layer
        Layer of `adata` to use as expression values.
    index
        List to set as index.
    columns
        List to set as columns names.
    sort_values
        Wether to sort values by first column (sort_values=True) or a specified column.
    dropna
        Drop columns/rows that contain NaNs in all ('all') or in any entry ('any').
    precision
        Set precision for pandas dataframe.

    Returns
    -------
    :class:`pd.DataFrame`
        A dataframe.
    """
    if precision is not None:
        pd.set_option("display.precision", precision)

    if isinstance(data, AnnData):
        keys, keys_split = (
            keys.split("*") if isinstance(keys, str) and "*" in keys else (keys, None)
        )
        keys, key_add = (
            keys.split("/") if isinstance(keys, str) and "/" in keys else (keys, None)
        )
        keys = [keys] if isinstance(keys, str) else keys
        key = keys[0]

        s_keys = ["obs", "var", "obsm", "varm", "uns", "layers"]
        d_keys = [
            data.obs.keys(),
            data.var.keys(),
            data.obsm.keys(),
            data.varm.keys(),
            data.uns.keys(),
            data.layers.keys(),
        ]

        if hasattr(data, "obsp") and hasattr(data, "varp"):
            s_keys.extend(["obsp", "varp"])
            d_keys.extend([data.obsp.keys(), data.varp.keys()])

        if keys is None:
            df = data.to_df()
        elif key in data.var_names:
            df = obs_df(data, keys, layer=layer)
        elif key in data.obs_names:
            df = var_df(data, keys, layer=layer)
        else:
            if keys_split is not None:
                keys = [
                    k
                    for k in list(data.obs.keys()) + list(data.var.keys())
                    if key in k and keys_split in k
                ]
                key = keys[0]
            s_key = [s for (s, d_key) in zip(s_keys, d_keys) if key in d_key]
            if len(s_key) == 0:
                raise ValueError(f"'{key}' not found in any of {', '.join(s_keys)}.")
            if len(s_key) > 1:
                logg.warn(f"'{key}' found multiple times in {', '.join(s_key)}.")

            s_key = s_key[-1]
            df = getattr(data, s_key)[keys if len(keys) > 1 else key]
            if key_add is not None:
                df = df[key_add]
            if index is None:
                index = (
                    data.var_names
                    if s_key == "varm"
                    else data.obs_names
                    if s_key in {"obsm", "layers"}
                    else None
                )
                if index is None and s_key == "uns" and hasattr(df, "shape"):
                    key_cats = np.array(
                        [
                            key
                            for key in data.obs.keys()
                            if is_categorical_dtype(data.obs[key])
                        ]
                    )
                    num_cats = [
                        len(data.obs[key].cat.categories) == df.shape[0]
                        for key in key_cats
                    ]
                    if np.sum(num_cats) == 1:
                        index = data.obs[key_cats[num_cats][0]].cat.categories
                        if (
                            columns is None
                            and len(df.shape) > 1
                            and df.shape[0] == df.shape[1]
                        ):
                            columns = index
            elif isinstance(index, str) and index in data.obs.keys():
                index = pd.Categorical(data.obs[index]).categories
            if columns is None and s_key == "layers":
                columns = data.var_names
            elif isinstance(columns, str) and columns in data.obs.keys():
                columns = pd.Categorical(data.obs[columns]).categories
    elif isinstance(data, pd.DataFrame):
        if isinstance(keys, str) and "*" in keys:
            keys, keys_split = keys.split("*")
            keys = [k for k in data.columns if keys in k and keys_split in k]
        df = data[keys] if keys is not None else data
    else:
        df = data

    if issparse(df):
        df = np.array(df.A)
    if columns is None and hasattr(df, "names"):
        columns = df.names

    df = pd.DataFrame(df, index=index, columns=columns)

    if dropna:
        df.replace("", np.nan, inplace=True)
        how = dropna if isinstance(dropna, str) else "any" if dropna is True else "all"
        df.dropna(how=how, axis=0, inplace=True)
        df.dropna(how=how, axis=1, inplace=True)

    if sort_values:
        sort_by = (
            sort_values
            if isinstance(sort_values, str) and sort_values in df.columns
            else df.columns[0]
        )
        df = df.sort_values(by=sort_by, ascending=False)

    return df


# TODO: Generalize to arbitrary modality
def get_initial_size(
    adata: AnnData, layer: Optional[str] = None, by_total_size: bool = False
) -> Optional[ndarray]:
    """Get initial counts per observation of a layer.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    layer
        Name of layer for which to retrieve initial size.
    by_total_size
        Whether or not to return the combined initial size of the spliced and unspliced
        layers.

    Returns
    -------
    np.ndarray
        Initial counts per observation in the specified layer.
    """
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


def get_modality(adata: AnnData, modality: Optional[str]) -> Union[ndarray, spmatrix]:
    """Extract data of one modality.

    Arguments:
    ---------
    adata
        Annotated data to extract modality from.
    modality
        Modality for which data is needed.

    Returns
    -------
    Union[ndarray, spmatrix]
        Retrieved modality from :class:`~anndata.AnnData` object.
    """
    if modality in ["X", None]:
        return adata.X
    elif modality in adata.layers.keys():
        return adata.layers[modality]
    elif modality in adata.obsm.keys():
        if isinstance(adata.obsm[modality], DataFrame):
            return adata.obsm[modality].values
        else:
            return adata.obsm[modality]


@deprecated_arg_names({"layer": "modality"})
def get_size(adata: AnnData, modality: Optional[str] = None) -> ndarray:
    """Get counts per observation in a modality.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    modality
        Name of modality for which to retrieve size.

    Returns
    -------
    np.ndarray
        Counts per observation in the specified modality.
    """
    X = get_modality(adata=adata, modality=modality)
    return sum(X, axis=1)


def make_dense(
    adata: AnnData, modalities: Union[List[str], str], inplace: bool = True
) -> Optional[AnnData]:
    """Densify sparse AnnData entry.

    Arguments:
    ---------
    adata
        Annotated data object.
    modality
        Modality to make dense.
    inplace
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

    Arguments:
    ---------
    adata
        Annotated data object.
    modality
        Modality to make sparse.
    inplace
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


def merge(
    adata: AnnData, ldata: AnnData, copy: bool = True, **kwargs
) -> Optional[AnnData]:
    """Merge two annotated data matrices.

    Arguments:
    ---------
    adata
        Annotated data matrix (reference data set).
    ldata
        Annotated data matrix (to be merged into adata).
    copy
        Boolean flag to manipulate original AnnData or a copy of it.

    Returns
    -------
    Optional[:class:`anndata.AnnData`]
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
        clean_obs_names(adata, **kwargs)
        clean_obs_names(ldata, **kwargs)
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


def obs_df(adata: AnnData, keys: List[str], layer: Optional[str] = None) -> DataFrame:
    """Extract layer as Pandas DataFrame indexed by observation.

    Arguments:
    ---------
    adata
        Annotated data matrix (reference data set).
    keys
        Variables for which to extract data.
    layer
        Name of layer to turn into a Pandas DataFrame.

    Returns
    -------
    DataFrame
        DataFrame indexed by observations. Columns correspond to variables of specified
        layer.
    """
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


# TODO: Generalize to arbitrary modality
def set_initial_size(adata: AnnData, layers: Optional[str] = None) -> None:
    """Set current counts per observation of a layer as its initial size.

    The initial size is only set if it does not already exist.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    layers
        Name of layers for which to calculate initial size.

    Returns
    -------
    None
    """
    if layers is None:
        layers = ["unspliced", "spliced"]
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

    Arguments:
    ---------
    adata
        Annotated data object.
    new_value
        New value of modality.
    modality
        Modality to overwrite with new value. Defaults to `None`.
    inplace
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


def show_proportions(
    adata: AnnData, layers: Optional[str] = None, use_raw: bool = True
) -> None:
    """Proportions of abundances of modalities in layers.

    The proportions are printed.

    Arguments:
    ---------
    adata
        Annotated data matrix.
    layers
        Layers to consider.
    use_raw
        Use initial sizes, i.e., raw data, to determine proportions.

    Returns
    -------
    None
    """
    if layers is None:
        layers = ["unspliced", "spliced", "ambiguous"]
    layers_keys = [key for key in layers if key in adata.layers.keys()]
    counts_layers = [sum(adata.layers[key], axis=1) for key in layers_keys]
    if use_raw:
        size_key, obs = "initial_size_", adata.obs
        counts_layers = [
            obs[size_key + layer] if size_key + layer in obs.keys() else c
            for layer, c in zip(layers_keys, counts_layers)
        ]

    counts_per_cell_sum = np.sum(counts_layers, 0)
    counts_per_cell_sum += counts_per_cell_sum == 0

    mean_abundances = [
        np.mean(counts_per_cell / counts_per_cell_sum)
        for counts_per_cell in counts_layers
    ]

    print(f"Abundance of {layers_keys}: {np.round(mean_abundances, 2)}")


def var_df(adata: AnnData, keys: List[str], layer: Optional[str] = None):
    """Extract layer as Pandas DataFrame indexed by features.

    Arguments:
    ---------
    adata
        Annotated data matrix (reference data set).
    keys
        Observations for which to extract data.
    layer
        Name of layer to turn into a Pandas DataFrame.

    Returns
    -------
    DataFrame
        DataFrame indexed by features. Columns correspond to observations of specified
        layer.
    """
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


# TODO: Deprecate function
# TODO: Find better function name
def verify_dtypes(adata: AnnData) -> None:
    """Verify that AnnData object is not corrupted.

    Arguments:
    ---------
    adata
        Annotated data matrix to check.

    Returns
    -------
    None
    """
    try:
        _ = adata[:, 0]
    except IndexError:
        uns = adata.uns
        adata.uns = {}
        try:
            _ = adata[:, 0]
            logg.warn(
                "Safely deleted unstructured annotations (adata.uns), \n"
                "as these do not comply with permissible anndata datatypes."
            )
        except IndexError:
            logg.warn(
                "The data might be corrupted. Please verify all annotation datatypes."
            )
            adata.uns = uns
