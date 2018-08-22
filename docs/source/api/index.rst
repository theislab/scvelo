.. automodule:: scvelo

API
===

Import scvelo and scanpy as::

   import scvelo as scv
   import scanpy.api as sc


Preprocessing (pp)
~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   preprocessing.show_proportions
   preprocessing.moments
   preprocessing.recipe_velocity


Tools (tl)
~~~~~~~~~~

.. autosummary::
   :toctree: .

   tools.velocity
   tools.velocity_graph
   tools.velocity_embedding


Plotting (pl)
~~~~~~~~~~~~~

.. autosummary::
   :toctree: .

   plotting.velocity
   plotting.velocity_embedding
   plotting.velocity_embedding_grid


Datasets
~~~~~~~~

.. autosummary::
   :toctree: .

   datasets.dentategyrus
   datasets.toy_data