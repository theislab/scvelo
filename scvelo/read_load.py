import os
import warnings
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd

from scvelo.core import clean_obs_names as _clean_obs_names
from scvelo.core import get_df as _get_df
from scvelo.core import merge as _merge
from scvelo.core._anndata import obs_df as _obs_df
from scvelo.core._anndata import var_df as _var_df


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
    warnings.warn(
        "`scvelo.read_load.var_df` is deprecated since scVelo v0.2.4 and will be "
        "removed in a future version. Please use `scvelo.core._anndata.var_df` "
        "instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return _var_df(adata=adata, keys=keys, layer=layer)


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
    warnings.warn(
        "`scvelo.read_load.get_df` is deprecated since scVelo v0.2.4 and will be "
        "removed in a future version. Please use `scvelo.core.get_df` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return _get_df(
        data=data,
        keys=keys,
        layer=layer,
        index=index,
        columns=columns,
        sort_values=sort_values,
        dropna=dropna,
        precision=precision,
    )


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
    df = df.drop_duplicates()
    return df


def convert_to_gene_names(ensembl_names=None):
    """Retrieve gene names from ensembl IDs."""
    df = load_biomart()
    if ensembl_names is not None:
        if isinstance(ensembl_names, str):
            ensembl_names = [ensembl_names]
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
