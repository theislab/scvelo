import os
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd


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
