.. automodule:: scvelo

API
===

Import scVelo as::

   import scvelo as scv


Read / Load
~~~~~~~~~~~

.. autosummary::
   :toctree: .

   read
   read_loom


Preprocessing (pp)
~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.filter_genes
   pp.filter_genes_dispersion
   pp.normalize_per_cell
   pp.filter_and_normalize
   pp.moments


Tools (tl)
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.recover_dynamics
   tl.velocity
   tl.velocity_graph
   tl.transition_matrix
   tl.velocity_embedding

   tl.terminal_states
   tl.latent_time
   tl.velocity_pseudotime

   tl.velocity_clusters
   tl.rank_velocity_genes
   tl.velocity_confidence


Plotting (pl)
~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pl.scatter
   pl.velocity
   pl.velocity_graph
   pl.velocity_embedding
   pl.velocity_embedding_grid
   pl.velocity_embedding_stream

Datasets
~~~~~~~~

.. autosummary::
   :toctree: .

   datasets.toy_data
   datasets.dentategyrus
   datasets.pancreatic_endocrinogenesis
   datasets.forebrain


Utils
~~~~~

.. autosummary::
   :toctree: .

   utils.show_proportions
   utils.cleanup
   utils.clean_obs_names
   utils.merge


Settings
~~~~~~~~

.. autosummary::
   :toctree: .

   settings.set_figure_params
