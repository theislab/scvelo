Getting Started
---------------

Welcome to scVelo!

scVelo is a scalable toolkit for estimating and analyzing RNA velocities in single cells.


Installation
^^^^^^^^^^^^

scVelo requires Python 3.6 or later. We recommend to use Miniconda_.

Install scVelo from PyPI_ using::

    pip install -U scvelo

or from source using::

    pip install git+https://github.com/theislab/scvelo


Parts of scVelo require (optional)::

    conda install -c conda-forge numba pytables louvain

The splicing data can be obtained using the `velocyto command line interface`_.

Basic Usage
^^^^^^^^^^^

Import scVelo as::

    import scvelo as scv

For beautiful visualization you can change the matplotlib settings to our defaults with::

    scv.settings.set_figure_params('scvelo')

Read your data
''''''''''''''
Read your data file (loom, h5ad, csv, ...) using::

    adata = scv.read(filename, cache=True)

which stores the data matrix (``adata.X``),
annotation of cells / observations (``adata.obs``) and genes / variables (``adata.var``), unstructured annotation such
as graphs (``adata.uns``) and additional data layers where spliced and unspliced counts are stored (``adata.layers``) .

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

If you already have an existing preprocessed adata object you can simply merge the spliced/unspliced counts via::

    ldata = scv.read(filename.loom, cache=True)
    adata = scv.utils.merge(adata, ldata)

If you do not have a datasets yet, you can still play around using one of the in-built datasets, e.g.::

    adata = scv.datasets.dentategyrus()

The typical workflow consists of subsequent calls of preprocessing (``scv.pp.*``), analysis tools (``scv.tl.*``) and plotting (``scv.pl.*``).

Basic preprocessing
'''''''''''''''''''

For velocity estimation basic preprocessing (i.e. gene selection and normalization) is sufficient, e.g. using::

    scv.pp.filter_and_normalize(adata, **params)

For velocity estimation we need the first- and second-order moments (basically means and variances), computed with::

    scv.pp.moments(adata, **params)

Velocity Tools
''''''''''''''

The core of the software is the efficient and robust estimation of velocities, obtained with::

    scv.tl.velocity(adata, mode='stochastic', **params)

The velocities are vectors in gene expression space obtained by solving a stochastic model of transcriptional dynamics.
The solution to the deterministic model is obtained by setting ``mode='deterministic'``.

The velocities are stored in ``adata.layers`` just like the count matrices.

Now we would like to predict cell transitions that are in accordance with the velocity directions. These are computed
using cosine correlation (i.e. find potential cell transitions that correlate with the velocity vector) and are stored
in a matrix called velocity graph::

    scv.tl.velocity_graph(adata, **params)

Using the graph you can then project the velocities into any embedding (such as UMAP, e.g. obtained with scanpy_)::

    scv.tl.velocity_embedding(adata, basis='umap', **params)

Note, that translation of velocities into a graph is only needed for non-linear embeddings.
In PCA space you can skip the velocity graph and directly project into the embedding using ``scv.tl.velocity_embedding(adata, basis='pca', direct_projection=True)``.


Visualization
'''''''''''''

Finally the velocities can be projected and visualized in any embedding (e.g. UMAP) on single cell level, grid level, or as streamplot::

    scv.pl.velocity_embedding(adata, basis='umap', **params)
    scv.pl.velocity_embedding_grid(adata, basis='umap', **params)
    scv.pl.velocity_embedding_stream(adata, basis='umap', **params)

For every tool module there is a plotting counterpart, which allows you to examine your results in detail, e.g.::

    scv.pl.velocity(adata, var_names=['gene_A', 'gene_B'], **params)
    scv.pl.velocity_graph(adata, **params)


.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/scvelo
.. _GitHub: https://github.com/theislab/scvelo
.. _scanpy: https://scanpy.readthedocs.io/en/latest/api
.. _`velocyto command line interface`: http://velocyto.org/velocyto.py/tutorial/cli.html
