.. role:: small
.. role:: smaller


Version 0.1.14 :small:`Dez 7, 2018`
------------------------------------
Plotting:

- pl.velocity_graph: scatterplot on the embedding
- added support for paga (from scanpy)
- velocity_embedding_stream: new function for visualizing cell differenciation paths

Preprocessing:

- more intelligent preprocessing now automatically detects which pp steps have already been done
- louvain renamed to clusters

Tools:

- cell_fate and terminal_states: compute final cell states
- velocity_pseudotime allows for pseudotemporal ordering of cells
- enable approximated group fitting


Version 0.1.11 :small:`Okt 27, 2018`
-----------------------------------
Plotting:

- adopting scanpy plot defaults: to use them call set_rcParams

Utils:

- merge function added for merging two datasets (e.g. loom into h5ad)
- improved performance by operations on sparse matrices

Tools:

- neighbors can now be computed on any embedding
- added option for weighted fitting in velocity estimation
- enable velocity estimation on subsets and with raw data
- score function added: velocity_confidence, -_transition and rank_velocity_genes, velocity_clusters
- improved the finding of roots
- only use informative genes (high abundance and variability) for velocity estimation
- transition matrix now includes negative cosines and self-transitions


Version 0.1.8 :small:`Sep 12, 2018`
-----------------------------------
Plotting:

- support saving plots as pdf, png etc.
- support multiple colors and layers
- quiver autoscaling for velocity plots
- attributes added: figsize and dpi

Preprocessing:

- filter_and_normalize() instead of recipe_velocity()
- normalization of layers is done automatically when computing moments

Tools:

- terminal_states: computes root and end points via eigenvalue decomposition :smaller:`thanks to M Lange`


Version 0.1.5 :small:`Sep 4, 2018`
----------------------------------
- Support writing loom files
- Support both dense and sparse layers
- Plotting bugfixes
- Added pp.recipe_velocity()

Version 0.1.2 :small:`Aug 21, 2018`
-----------------------------------
First alpha release of scvelo.