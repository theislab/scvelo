Usage Principles
----------------

Import scvelo as::

    import scvelo as scv


Read the file into a data object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read the data file (loom, h5ad, xlsx, csv, tab, txt ...) to an :class:`~anndata.AnnData` object::

   adata = scv.read(filename, cache=True)

which stores the data matrix (``adata.X``) with dimension :math:`n_{\mathrm{obs}} \times n_{\mathrm{vars}}`,
annotation of observations (``adata.obs``) and variables (``adata.var``), unstructured annotation (``adata.uns``) and
additional layers (``adata.layers``).

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

For instance, the data matrices relevant for velocity analysis can be retrieved via ``adata.layers['spliced']`` and ``adata.layers['unspliced']``.

Basic preprocessing
^^^^^^^^^^^^^^^^^^^

If you have not preprocessed you data yet you may simply run the following line which selects genes (according to detection and
variability) and normalizes each cell by total counts::

    scv.pp.filter_and_normalize(adata, **params)

You may consider using scanpy_ for further preprocessing (such as correcting for batch effects).

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

Compute the velocity graph which finds the most likely cell transitions according to the velocity prediction
(based on cosine correlation) and serves as a Markov transition matrix::

   scv.tl.velocity_graph(adata, **params)


Using the graph you can easily project the velocities into any embedding (such as UMAP, e.g. obtained with scanpy_)::

   scv.tl.velocity_embedding(adata, basis='umap', **params)


Visualization
^^^^^^^^^^^^^
The velocities for all individual cells can be visualized using::

   scv.pl.velocity_embedding(adata, basis='umap', **params)

For big datasets it might be useful to visualize the velocities on a grid::

   scv.pl.velocity_embedding_grid(adata, basis='umap', **params)


.. _scanpy: https://scanpy.readthedocs.io/en/latest/api