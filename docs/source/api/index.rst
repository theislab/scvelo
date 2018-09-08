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

   datasets.dentategyrus
   datasets.toy_data