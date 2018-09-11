.. automodule:: docs.source.api

API
===

Import scvelo and scanpy as::

   import scvelo as scv
   import scanpy.api as sc


Preprocessing (pp)
~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pp.show_proportions
   pp.filter_and_normalize
   pp.moments


Tools (tl)
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tl.velocity
   tl.velocity_graph
   tl.velocity_embedding

   tl.transition_matrix
   tl.terminal_states
   tl.score_smoothness
   tl.score_transition


Plotting (pl)
~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   pl.velocity
   pl.velocity_embedding
   pl.velocity_embedding_grid


Datasets
~~~~~~~~

.. autosummary::
   :toctree: .

   datasets.toy_data
   datasets.dentategyrus
   datasets.forebrain
