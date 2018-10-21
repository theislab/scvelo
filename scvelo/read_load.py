import os, re
import numpy as np
import pandas as pd
from urllib.request import urlretrieve
from pathlib import Path
from scanpy.api import AnnData, read, read_loom


def load(filename, backup_url=None, **kwargs):
    numpy_ext = {'npy', 'npz'}
    pandas_ext = {'csv', 'txt'}

    if not os.path.exists(filename) and backup_url is None:
        raise FileNotFoundError('Did not find file {}.'.format(filename))

    elif not os.path.exists(filename):
        d = os.path.dirname(filename)
        if not os.path.exists(d): os.makedirs(d)
        urlretrieve(backup_url, filename)

    ext = Path(filename).suffixes[-1][1:]

    if ext in numpy_ext: return np.load(filename, **kwargs)
    elif ext in pandas_ext: return pd.read_csv(filename, **kwargs)
    else: raise ValueError('"{}" does not end on a valid extension.\n'
                           'Please, provide one of the available extensions.\n{}\n'
                           .format(filename, numpy_ext|pandas_ext))


def clean_obs_names(adata, base='[AGTCBDHKMNRSVWY]', ID_length=12, copy=False):
    """clean in obs_names and append the rest to obs.sample as follows:

        obs_name "SomePrefix_AGTGATCCATG_SomeSuffix" is changed into
        obs_name "AGTGATCCATG" and obs.sample "SomePrefix_SomeSuffix"

        according to https://www.neb.com/tools-and-resources/usage-guidelines/the-genetic-code
    """
    # change cell ID "SomePrefix_AGTGATCCATG_SomeSuffix" into:
    # cell ID "AGTGATCCATG" of sample "SomePrefix_SomeSuffix"
    # GENETIC CODE: according to https://www.neb.com/tools-and-resources/usage-guidelines/the-genetic-code
    def get_base_list(name, base):
        base_list = base
        while re.search(base_list + base, name) is not None:
            base_list += base
        if len(base_list) == 0:
            raise ValueError('Encountered an invalid ID in obs_names: ', name)
        return base_list

    names = adata.obs_names
    base_list = get_base_list(names[0], base)

    if len(np.unique([len(name) for name in adata.obs_names])) == 1:
        start, end = re.search(base_list, names[0]).span()
        newIDs = [name[start:end] for name in names]
        start, end = 0, len(newIDs[0])
        for i in range(end - ID_length):
            if np.any([ID[i] not in base for ID in newIDs]): start += 1
            if np.any([ID[::-1][i] not in base for ID in newIDs]): end -= 1

        newIDs = [ID[start:end] for ID in newIDs]
        prefixes = [names[i].replace(newIDs[i], '') for i in range(len(names))]
    else:
        prefixes, newIDs = [], []
        for name in names:
            match = re.search(base_list, name)
            newID = re.search(get_base_list(name, base), name).group() if match is None else match.group()
            newIDs.append(newID)
            prefixes.append(name.replace(newID, ''))

    adata.obs_names = newIDs
    if len(prefixes[0]) > 0 and len(np.unique(prefixes)) > 1:
        #idx_names = np.random.choice(len(names), size=20, replace=False)
        #for i in range(len(names[0])):
        #    if np.all([re.search(names[0][:i], names[ix]) for ix in idx_names]) is not None: obs_key = names[0][:i]
        adata.obs['sample_batch'] = prefixes
        adata._sanitize()

    adata.obs_names_make_unique()
    return adata if copy else None


def merge(adata, ldata):
    common_obs = adata.obs_names.intersection(ldata.obs_names)

    if len(common_obs) == 0:
        clean_obs_names(adata)
        clean_obs_names(ldata)
        common_obs = adata.obs_names.intersection(ldata.obs_names)

    _adata = adata.copy() if adata.shape[1] >= ldata.shape[1] else ldata.copy()
    _ldata = ldata.copy() if adata.shape[1] >= ldata.shape[1] else adata.copy()

    _adata = _adata[common_obs]
    _ldata = _ldata[common_obs]

    for attr in _ldata.obs.keys():
        _adata.obs[attr] = _ldata.obs[attr]
    for attr in _ldata.obsm.keys():
        _adata.obsm[attr] = _ldata.obsm[attr]
    for attr in _ldata.uns.keys():
        _adata.uns[attr] = _ldata.uns[attr]
    for attr in _ldata.layers.keys():
        _adata.layers[attr] = _ldata.layers[attr]

    if _adata.shape[1] == _ldata.shape[1]:
        if np.all(adata.var_names == ldata.var_names):
            for attr in _ldata.var.keys():
                _adata.var[attr] = _ldata.var[attr]
            for attr in _ldata.varm.keys():
                _adata.varm[attr] = _ldata.varm[attr]
        else:
            raise ValueError('Variable names are not identical.')

    return _adata



