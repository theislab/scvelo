import os, numpy, pandas
from urllib.request import urlretrieve
from pathlib import Path
from scanpy.api import read


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

    if ext in numpy_ext: return numpy.load(filename, **kwargs)
    elif ext in pandas_ext: return pandas.read_csv(filename, **kwargs)
    else: raise ValueError('"{}" does not end on a valid extension.\n'
                           'Please, provide one of the available extensions.\n{}\n'
                           .format(filename, numpy_ext|pandas_ext))
