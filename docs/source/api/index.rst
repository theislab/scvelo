.. automodule:: docs.source.api

API
===

Import scveloas::

   import scvelo as scv


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

   pl.scatter
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
