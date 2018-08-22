Usage Principles
----------------

Import scvelo (velocity specific workflows) and scanpy (basic workflows) as::

    import scvelo as scv
    import scanpy.api as sc

Read the file into a data object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read the data file (loom, h5ad, xlsx, csv, tab, txt ...) to an :class:`~anndata.AnnData` object::

   adata = sc.read(filename, cache=True)


which stores the data matrix (``adata.X``) with dimension :math:`n_{\mathrm{obs}} \times n_{\mathrm{vars}}`,
annotation of observations (``adata.obs``) and variables (``adata.var``), unstructured annotation (``adata.uns``) and
additional layers (``adata.layers``).

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

For instance, the data matrices relevant for velocity analysis can be retrieved via ``adata.layers['spliced']`` and ``adata.layers['unspliced']``.

Basic preprocessing
^^^^^^^^^^^^^^^^^^^
Select genes (according to detection and variability) and normalize each cell by total counts::

    sc.pp.filter_genes(adata, **params)
    sc.pp.filter_genes_dispersion(adata, **params)

    sc.pp.normalize_per_cell(adata, layers='all', **params)

Then compute the neighbor graph in PCA space on the logarithmized data::

   sc.pp.log1p(adata, **params)
   sc.pp.pca(adata, **params)
   sc.pp.neighbors(adata, use_rep='X_pca', **params)

Finally, compute the moments of spliced and unspliced counts which will be used for velocity estimation::

   scv.pp.moments(adata, **params)

That's all, no extensive preparation needed and it takes a few seconds to run these lines.

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


Using the graph you can easily project the velocities into any embedding (such as UMAP)::

   scv.tl.velocity_embedding(adata, basis='umap', **params)


Visualization
^^^^^^^^^^^^^
The velocities for all individual cells can be visualized using::

   scv.pl.velocity_embedding(adata, basis='umap', **params)

For big datasets it might be useful to visualize the velocities on a grid::

   scv.pl.velocity_embedding_grid(adata, basis='umap', **params)

