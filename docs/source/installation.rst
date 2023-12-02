Installation
------------

scVelo requires Python 3.6 or later. We recommend to use Miniconda_.

PyPI
^^^^

Install scVelo from PyPI_ using::

    pip install -U scvelo

``-U`` is short for ``--upgrade``.
If you get a ``Permission denied`` error, use ``pip install -U scvelo --user`` instead.


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    pip install git+https://github.com/theislab/scvelo@main

or::

    git clone https://github.com/theislab/scvelo && cd scvelo
    git checkout --track origin/main
    pip install -e .

``-e`` is short for ``--editable`` and links the package to the original cloned
location such that pulled changes are also reflected in the environment.

To contribute to scVelo, ``cd`` into the cloned directory and
install the latest packages required for development together with the pre-commit hooks::

    pip install -e ".[dev]"
    pre-commit install


Dependencies
^^^^^^^^^^^^

- `anndata <https://anndata.readthedocs.io/>`_ - annotated data object.
- `scanpy <https://scanpy.readthedocs.io/>`_ - toolkit for single-cell analysis.
- `numpy <https://docs.scipy.org/>`_, `scipy <https://docs.scipy.org/>`_, `pandas <https://pandas.pydata.org/>`_, `scikit-learn <https://scikit-learn.org/>`_, `matplotlib <https://matplotlib.org/>`_.


Parts of scVelo (directed PAGA and Louvain modularity) require (optional)::

    pip install igraph louvain


Using fast neighbor search via `hnswlib <https://github.com/nmslib/hnswlib>`_ further requires (optional)::

    pip install pybind11 hnswlib


Jupyter Notebook
^^^^^^^^^^^^^^^^

To run the tutorials in a notebook locally, please install::

   conda install notebook

and run ``jupyter notebook`` in the terminal. If you get the error ``Not a directory: 'xdg-settings'``,
use ``jupyter notebook --no-browser`` instead and open the url manually (or use this
`bugfix <https://github.com/jupyter/notebook/issues/3746#issuecomment-444957821>`_).


If you run into issues, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/scvelo
.. _Github: https://github.com/theislab/scvelo
.. _`Github issue`: https://github.com/theislab/scvelo/issues/new/choose
