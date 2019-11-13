Getting Started
---------------

Welcome to scVelo!

scVelo is a scalable toolkit for estimating and analyzing RNA velocities in single cells.


Installation
^^^^^^^^^^^^
scVelo requires Python 3.6 or later. We recommend to use Miniconda_.

Install scVelo from PyPI using::

    pip install -U scvelo


To work with the latest development version, install from source using::

    pip install git+https://github.com/theislab/scvelo

or::

    git clone git+https://github.com/theislab/scvelo
    pip install -e scvelo

Parts of scVelo require (optional)::

    conda install -c conda-forge numba pytables louvain



Alignment
^^^^^^^^^
The splicing data can be obtained using one of the following read counting pipelines:

- `velocyto pipeline`_
- `kallisto pipeline via loompy`_
- `kallisto pipeline via kb`_

scVelo in action
^^^^^^^^^^^^^^^^
Import scvelo as::

    import scvelo as scv

For beautified visualization you can change the matplotlib settings to our defaults with::

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
After basic preprocessing (gene selection and normalization is sufficient),
we compute the first- and second-order moments (basically means and variances) for velocity estimation::

    scv.pp.filter_and_normalize(adata, **params)
    scv.pp.moments(adata, **params)

Velocity Tools
''''''''''''''
The core of the software is the efficient and robust estimation of velocities, obtained with::

    scv.tl.velocity(adata, mode='stochastic', **params)

The velocities are vectors in gene expression space obtained by solving a stochastic model of transcriptional dynamics.
The solution to the deterministic model is obtained by setting ``mode='deterministic'``.

The solution to the dynamical model is obtained by setting ``mode='dynamical'``, which requires to run
``scv.tl.recover_dynamics(adata, **params)`` beforehand.

The velocities are stored in ``adata.layers`` just like the count matrices.

The velocities are projected into a lower-dimensional embedding by translating them into likely cell transitions.
That is, for each velocity vector we find the likely cell transitions that are accordance with that direction.
The probabilities of one cell transitioning into another cell are computed using cosine correlation
(btw. the potential cell transition and the velocity vector) and are stored in a matrix denoted as velocity graph::

    scv.tl.velocity_graph(adata, **params)

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
.. _`velocyto pipeline`: http://velocyto.org/velocyto.py/tutorial/cli.html
.. _`kallisto pipeline via loompy`: https://linnarssonlab.org/loompy/kallisto/index.html
.. _`kallisto pipeline via kb`: https://www.kallistobus.tools/kb_velocity_tutorial.html