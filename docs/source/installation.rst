Installation
------------

scVelo requires Python 3.6 or later. We recommend to use Miniconda_.

PyPI
^^^^

Install scVelo from PyPI_ using::

    pip install -U scvelo


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    pip install git+https://github.com/theislab/scvelo

or::

    git clone git+https://github.com/theislab/scvelo
    pip install -e scvelo


Dependencies
^^^^^^^^^^^^

- `anndata <https://anndata.readthedocs.io/>`_ - object to keep track of data with annotations.
- `scanpy <https://scanpy.readthedocs.io/>`_ - toolkit for analyzing single-cell data.
- `numpy <https://docs.scipy.org/>`_, `scipy <https://docs.scipy.org/>`_, `pandas <https://pandas.pydata.org/>`_, `scikit-learn <https://scikit-learn.org/>`_ and `matplotlib <https://matplotlib.org/>`_.


Parts of scVelo require (optional)::

    conda install -c conda-forge numba pytables louvain


If you run into issues, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/scvelo
.. _Github: https://github.com/theislab/scvelo
.. _`Github issue`: https://github.com/theislab/scvelo/issues/new/choose