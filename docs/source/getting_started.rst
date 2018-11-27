Getting Started
---------------

Welcome to scVelo!

scVelo is a scalable toolkit for estimating and analyzing stochastic RNA velocities in single cells.

If you don't have working Python 3.6 yet, consider installing Miniconda_.

Once you are set, install scVelo from PyPI_ with::

    pip install scvelo

If you want to work with the latest version on GitHub_, install scVelo from source::

    git clone https://github.com/theislab/scvelo.git
    cd scvelo
    pip install .

Basic Usage
-----------

Import scVelo as::

    import scvelo as scv

For beautiful visualization you can change the matplotlib settings to our defaults with::

    scv.settings.set_figure_params('scvelo')

Read your data into an object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read your data file (loom, h5ad, xlsx, csv, tab, txt ...) to an :class:`~anndata.AnnData` object::

    adata = scv.read(filename, cache=True)

which stores the data matrix (``adata.X``) with dimension :math:`n_{\mathrm{obs}} \times n_{\mathrm{vars}}`,
annotation of observations (``adata.obs``) and variables (``adata.var``), unstructured annotation (``adata.uns``) and
additional layers (``adata.layers``).

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

For instance, the data matrices relevant for velocity analysis can be retrieved via
``adata.layers['spliced']`` and ``adata.layers['unspliced']``.

he typical workflow consists of subsequent calls of preprocessing (``scv.pp.*``), analysis tools (``scv.tl.*``) and plotting (``scv.pl.*``).

Basic preprocessing
^^^^^^^^^^^^^^^^^^^

You are probably familiar with preprocessing. The very basic steps include gene selection by detection and variability, and normalization of each cell by total counts.
Simply run::

    scv.pp.filter_and_normalize(adata, **params)

I recommend using scanpy_ (which perfectly harmonizes with scVelo) to explore further preprocessing steps (such as correcting for batch effects).

For processing of spliced and unspliced counts it suffices to compute their moments (which automatically normalizes the counts)::

    scv.pp.moments(adata, **params)

That's all, no extensive preparation is needed.

Velocity Tools
^^^^^^^^^^^^^^

Now you are hitting the core of the package.

Estimating the velocities for each individual cell is done in a single line::

    scv.tl.velocity(adata, mode='stochastic', **params)

The velocities are vectors in gene expression space obtained by using a closed-form solution that
solves a stochastic model of transcriptional dynamics. The stochastic model incorporates intrinsic expression variability.
The solution to the deterministic model is obtained by setting mode to 'deterministic'.

The velocities are stored in ``adata.layers`` just like the count matrices.

Given these velocities we are interested in cell transitions that are likely. These are computed using cosine correlation
(i.e. find potential transitions that correlate with the velocity vector) and are stored in a matrix that we call velocity graph::

    scv.tl.velocity_graph(adata, **params)

Using the graph you can then project the velocities into any embedding (such as UMAP, e.g. obtained with scanpy_)::

    scv.tl.velocity_embedding(adata, basis='umap', **params)

Visualization
^^^^^^^^^^^^^
The velocities for all individual cells can be visualized using::

    scv.pl.velocity_embedding(adata, basis='umap', **params)

For big datasets it might be useful to visualize the velocities on a grid::

    scv.pl.velocity_embedding_grid(adata, basis='umap', **params)

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/scvelo
.. _GitHub: https://github.com/theislab/scvelo
.. _scanpy: https://scanpy.readthedocs.io/en/latest/api
