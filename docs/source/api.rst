.. automodule:: scvelo

API
===

Import scVelo as::

   import scvelo as scv


After reading the data (``scv.read``) or loading an in-built dataset (``scv.datasets.*``),
the typical workflow consists of subsequent calls of
preprocessing (``scv.pp.*``), analysis tools (``scv.tl.*``) and plotting (``scv.pl.*``).
Further, several utilities (``scv.utils.*``) are provided to facilitate data analysis.


Read / Load
-----------

.. autosummary::
   :toctree: .

   read
   read_loom


Preprocessing (pp)
------------------

**Basic preprocessing** (gene selection and normalization)

.. autosummary::
   :toctree: .

   pp.filter_genes
   pp.filter_genes_dispersion
   pp.normalize_per_cell
   pp.log1p
   pp.filter_and_normalize

**Moments** (across nearest neighbors in PCA space)

.. autosummary::
   :toctree: .

   pp.pca
   pp.neighbors
   pp.moments


Tools (tl)
----------

**Clustering and embedding**
(more at `scanpy-docs <https://scanpy.readthedocs.io/en/stable/api/>`_)

.. autosummary::
   :toctree: .

   tl.louvain
   tl.umap

**Velocity estimation**

.. autosummary::
   :toctree: .

   tl.velocity
   tl.velocity_graph
   tl.velocity_embedding

**Dynamical modeling**

.. autosummary::
   :toctree: .

   tl.recover_dynamics
   tl.differential_kinetic_test

**Dynamical genes**

.. autosummary::
   :toctree: .

   tl.rank_velocity_genes
   tl.rank_dynamical_genes


**Pseudotime and trajectory inference**

.. autosummary::
   :toctree: .

   tl.terminal_states
   tl.velocity_pseudotime
   tl.latent_time
   tl.paga

**Further tools**

.. autosummary::
   :toctree: .

   tl.velocity_clusters
   tl.velocity_confidence
   tl.score_genes_cell_cycle


Plotting (pl)
-------------

**Base scatter plot**

.. autosummary::
   :toctree: .

   pl.scatter

**Velocity embeddings**

.. autosummary::
   :toctree: .

   pl.velocity_embedding
   pl.velocity_embedding_grid
   pl.velocity_embedding_stream

**Velocity graph**

.. autosummary::
   :toctree: .

   pl.velocity
   pl.velocity_graph
   pl.paga

**Further plotting**

.. autosummary::
   :toctree: .

   pl.proportions
   pl.heatmap
   pl.hist


Datasets
--------

.. autosummary::
   :toctree: .

   datasets.pancreas
   datasets.dentategyrus
   datasets.forebrain
   datasets.simulation


Utils
-----

**Get data by key**

.. autosummary::
   :toctree: .

   get_df

**Data preparation**

.. autosummary::
   :toctree: .

   utils.cleanup
   utils.clean_obs_names
   utils.merge
   utils.show_proportions

**Getters**

.. autosummary::
   :toctree: .

   utils.get_moments
   utils.get_transition_matrix
   utils.get_cell_transitions
   utils.get_extrapolated_state

**Least squares and correlation**

.. autosummary::
   :toctree: .

   utils.leastsq
   utils.vcorrcoef
   utils.test_bimodality


Settings
--------

.. autosummary::
   :toctree: .

   set_figure_params
