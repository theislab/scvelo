Getting Started
---------------

Here, you will be briefly guided through the basics of how to use scVelo.
Once you are set, the following tutorials go straight into analysis of RNA velocity,
latent time, driver identification and many more.

First of all, the input data for scVelo are two count matrices of pre-mature (unspliced) and mature (spliced) abundances,
which can be obtained from standard sequencing protocols, using the `velocyto`_ or `kallisto`_
counting pipeline.

scVelo workflow at a glance
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Import scvelo as::

    import scvelo as scv

For beautified visualization you can change the matplotlib settings to our defaults with::

    scv.set_figure_params()

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

    adata = scv.datasets.pancreas()

The typical workflow consists of subsequent calls of preprocessing (``scv.pp.*``), analysis tools (``scv.tl.*``) and plotting (``scv.pl.*``).

Basic preprocessing
'''''''''''''''''''
After basic preprocessing (gene selection and normalization),
we compute the first- and second-order moments (means and uncentered variances) for velocity estimation::

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
That is, for each velocity vector we find the likely cell transitions that are in accordance with that direction.
The probabilities of one cell transitioning into another cell are computed using cosine correlation
(between the potential cell transition and the velocity vector) and are stored in a matrix denoted as velocity graph::

    scv.tl.velocity_graph(adata, **params)

Visualization
'''''''''''''

Finally, the velocities can be projected and visualized in any embedding (e.g. UMAP) on single cell level, as gridlines, or as streamlines::

    scv.pl.velocity_embedding(adata, basis='umap', **params)
    scv.pl.velocity_embedding_grid(adata, basis='umap', **params)
    scv.pl.velocity_embedding_stream(adata, basis='umap', **params)

For every tool module there is a plotting counterpart, which allows you to examine your results in detail, e.g.::

    scv.pl.velocity(adata, var_names=['gene_A', 'gene_B'], **params)
    scv.pl.velocity_graph(adata, **params)


.. _`velocyto`: http://velocyto.org/velocyto.py/tutorial/cli.html
.. _`kallisto`: https://www.kallistobus.tools/tutorials/kb_velocity/python/kb_velocity/#generate-rna-velocity-count-matrices
