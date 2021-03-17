import os
import warnings
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse

from anndata import AnnData

from scvelo.core import clean_obs_names as _clean_obs_names
from scvelo.core import merge as _merge
from scvelo.core._anndata import obs_df as _obs_df
from . import logging as logg


def load(filename, backup_url=None, header="infer", index_col="infer", **kwargs):
    """Load a csv, txt, tsv or npy file."""
    numpy_ext = {"npy", "npz"}
    pandas_ext = {"csv", "txt", "tsv"}

    if not os.path.exists(filename) and backup_url is None:
        raise FileNotFoundError(f"Did not find file {filename}.")

    elif not os.path.exists(filename):
        d = os.path.dirname(filename)
        if not os.path.exists(d):
            os.makedirs(d)
        urlretrieve(backup_url, filename)

    ext = Path(filename).suffixes[-1][1:]

    if ext in numpy_ext:
        return np.load(filename, **kwargs)

    elif ext in pandas_ext:
        df = pd.read_csv(
            filename,
            header=header,
            index_col=None if index_col == "infer" else index_col,
            **kwargs,
        )
        if index_col == "infer" and len(df.columns) > 1:
            is_int_index = all(np.arange(0, len(df)) == df.iloc[:, 0])
            is_str_index = isinstance(df.iloc[0, 0], str) and all(
                [not isinstance(d, str) for d in df.iloc[0, 1:]]
            )
            if is_int_index or is_str_index:
                df.set_index(df.columns[0], inplace=True)
        return df

    else:
        raise ValueError(
            f"'{filename}' does not end on a valid extension.\n"
            "Please, provide one of the available extensions.\n"
            f"{numpy_ext | pandas_ext}\n"
        )


read_csv = load


def clean_obs_names(data, base="[AGTCBDHKMNRSVWY]", ID_length=12, copy=False):
    warnings.warn(
        "`scvelo.read_load.clean_obs_names` is deprecated since scVelo v0.2.4 and will "
        "be removed in a future version. Please use `scvelo.core.clean_obs_names` "
        "instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return _clean_obs_names(data=data, base=base, ID_length=ID_length, copy=copy)


def merge(adata, ldata, copy=True):
    warnings.warn(
        "`scvelo.read_load.merge` is deprecated since scVelo v0.2.4 and will be "
        "removed in a future version. Please use `scvelo.core.merge` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return _merge(adata=adata, ldata=ldata, copy=True)


def obs_df(adata, keys, layer=None):
    warnings.warn(
        "`scvelo.read_load.obs_df` is deprecated since scVelo v0.2.4 and will be "
        "removed in a future version. Please use `scvelo.core._anndata.obs_df` "
        "instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return _obs_df(adata=adata, keys=keys, layer=layer)


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


def get_df(
    data,
    keys=None,
    layer=None,
    index=None,
    columns=None,
    sort_values=None,
    dropna="all",
    precision=None,
):
    """Get dataframe for a specified adata key.

    Return values for specified key
    (in obs, var, obsm, varm, obsp, varp, uns, or layers) as a dataframe.

    Arguments
    ------
    adata
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
    A dataframe.
    """
    if precision is not None:
        pd.set_option("precision", precision)

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
            df = _obs_df(data, keys, layer=layer)
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

    if hasattr(data, "var_names"):
        if df.index[0] in data.var_names:
            df.var_names = df.index
        elif df.columns[0] in data.var_names:
            df.var_names = df.columns
    if hasattr(data, "obs_names"):
        if df.index[0] in data.obs_names:
            df.obs_names = df.index
        elif df.columns[0] in data.obs_names:
            df.obs_names = df.columns

    return df


DataFrame = get_df


def load_biomart():
    # human genes from https://biomart.genenames.org/martform
    # mouse genes from http://www.ensembl.org/biomart/martview
    # antibodies from https://www.biolegend.com/en-us/totalseq
    nb_url = "https://github.com/theislab/scvelo_notebooks/raw/master/"

    filename = "data/biomart/mart_export_human.txt"
    df = load(filename, sep="\t", backup_url=f"{nb_url}{filename}")
    df.columns = ["ensembl", "gene name"]
    df.index = df.pop("ensembl")

    filename = "data/biomart/mart_export_mouse.txt"
    df2 = load(filename, sep="\t", backup_url=f"{nb_url}{filename}")
    df2.columns = ["ensembl", "gene name"]
    df2.index = df2.pop("ensembl")

    df = pd.concat([df, df2])
    return df


def convert_to_gene_names(ensembl_names=None):
    """Retrieve gene names from ensembl IDs."""
    df = load_biomart()
    if ensembl_names is not None:
        if isinstance(ensembl_names, str):
            ensembl_names = ensembl_names
        valid_names = [name for name in ensembl_names if name in df.index]
        if len(valid_names) > 0:
            df = df.loc[valid_names]

        gene_names = np.array(ensembl_names)
        idx = pd.DataFrame(ensembl_names).isin(df.index).values.flatten()
        gene_names[idx] = df["gene name"].values

        df = pd.DataFrame([ensembl_names, gene_names]).T
        df.columns = ["ensembl", "gene name"]
        df.index = df.pop("ensembl")
    return df


def convert_to_ensembl(gene_names=None):
    """Retrieve ensembl IDs from a list of gene names."""
    df = load_biomart()
    if gene_names is not None:
        if isinstance(gene_names, str):
            gene_names = [gene_names]
        valid_names = [name for name in gene_names if name in df["gene name"].tolist()]
        if len(valid_names) > 0:
            index = [i in valid_names for i in df["gene name"].tolist()]
            df = df[index]

        df["ensembl"] = df.index
        df = df.set_index("gene name")
    return df


def gene_info(name, fields="name,symbol,refseq,generif,ensembl"):
    """Retrieve gene information from biothings client."""
    try:
        from biothings_client import get_client
    except ImportError:
        raise ImportError(
            "Please install Biothings first via `pip install biothings_client`."
        )

    class MyGeneInfo(get_client("gene", instance=False)):
        def __init__(self):
            super(MyGeneInfo, self).__init__()

    if not name.startswith("ENS"):
        df = convert_to_gene_names()
        df.reset_index(inplace=True)
        df.set_index("gene name", inplace=True)
        if name in df.index:
            name = df.loc[name][0]

    info = MyGeneInfo().getgene(name, fields)
    return info
