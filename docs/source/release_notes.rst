.. role:: small
.. role:: smaller

Release Notes
=============

Version 0.2.5 :small:`Oct 14, 2022`
-----------------------------------

Changes:

- Catch non-positive parameter values and raise a `ValueError` if necessary (`PR 614 <https://github.com/theislab/scvelo/pull/614>`_).
- `get_mean_var` uses the same size parameter for mean and variance (`PR 698 <https://github.com/theislab/scvelo/pull/698>`_).

Bugfixes:

- `filter_genes` now works with `adata.layers['unspliced']` being sparse and `adata.layers['spliced']` dense (`PR 537 <https://github.com/theislab/scvelo/pull/537>`_).
- `show_proportions` actually considers the layer `"ambiguous"` if present (`PR 587 <https://github.com/theislab/scvelo/pull/587>`_).
- Fix calculation of Pearson's correlation in `csr_vcorrcoef` (`PR 679 <https://github.com/theislab/scvelo/pull/679>`_).
- Fix `get_mean_var` to work with sparse input and `ignore_zeros=True` (`PR 698 <https://github.com/theislab/scvelo/pull/698>`_).
- Fix bug in neighbor calculation (`PR 797 <https://github.com/theislab/scvelo/pull/797>`_).
- Fix `optimization.py::get_weight` to work with numeric, non-integer values (`PR 839 <https://github.com/theislab/scvelo/pull/839>`_).
- Fix inference with `fit_scaling=False` (`PR 848 <https://github.com/theislab/scvelo/pull/848>`_).
- Fix saving of velocity embedding stream (`PR 900 <https://github.com/theislab/scvelo/pull/900>`_).
- Fix Pandas' display precison when passed to `get_df` (`PR 907 <https://github.com/theislab/scvelo/pull/907>`_).

Version 0.2.4 :small:`Aug 26, 2021`
-----------------------------------

Perspectives:

- Landing page and two notebooks accompanying the perspectives manuscript at MSB.
- New datasets: Gastrulation, bone marrow, and PBMCs.

New capabilities:

- Added vignettes accompanying the NBT manuscript.
- Kinetic simulations with time-dependent rates.
- New arguments for `tl.velocity_embedding_stream` (`PR 492 <https://github.com/theislab/scvelo/pull/492>`_).
- Introduced automated code formatting `flake8` and `isort` (`PR 360 <https://github.com/theislab/scvelo/pull/360>`_, `PR 374 <https://github.com/theislab/scvelo/pull/374>`_).
- `tl.velocity_graph` parallelized (`PR 392 <https://github.com/theislab/scvelo/pull/392>`_).
- `legend_align_text` parameter in `pl.scatter` for smart placing of labels without overlapping.
- Save option for `pl.proportions`.

Bugfixes:

- Pinned `sphinx<4.0` and `nbsphinx<0.8.7`.
- Fix IPython import at CLI.


Version 0.2.3 :small:`Feb 13, 2021`
-----------------------------------

- `tl.recover_dynamics`: Multicore implementation :smaller:`thanks to M Klein, Y Schaelte, P Weiler`
- CI now runs on GitHub Actions

New utility functions:

- `utils.gene_info`: Retrieve gene information from biothings client.
- `utils.convert_to_ensembl` and `utils.convert_to_gene_names`: Converting ensembl IDs into gene names and vice versa.


Version 0.2.2 :small:`July 22, 2020`
------------------------------------

- `tl.paga`: PAGA graph with velocity-directed edges.
- Black code style


Version 0.2.1 :small:`May 28, 2020`
--------------------------------------
Bugfixes:

- Correct identification of root cells in `tl.latent_time` :smaller:`thanks to M Lange`
- Correct usage of latent_time prior in `tl.paga` :smaller:`thanks to G Lubatti`


Version 0.2.0 :small:`May 12, 2020`
--------------------------------------
New vignettes:

- RNA velocity basics
- Dynamical Modeling
- Differential Kinetics

Tools:

- `tl.differential_kinetic_test`: introduced a statistical test to detect different kinetic regimes.
- `tl.rank_dynamical_genes`: introduced a gene ranking by cluster-wise likelihoods.
- `tl.paga`: introduced directed PAGA graph

Plotting:

- `pl.scatter` enhancements: linear and polynomical fits, gradient coloring
- `pl.proportions`: Pie and bar chart of spliced/unspliced proprtions.
- `GridSpec`: multiplot environment.


Version 0.1.20 :small:`Sep 5, 2019`
-----------------------------------
Tools:

- `tl.recover_dynamics`: introduced a dynamical model inferring the full splicing kinetics, thereby identifying all kinetic rates of transcription, splicing and degradation.
- `tl.recover_latent_time`: infers a shared latent time across all genes based on the learned splicing dynamics.

Plotting:

- `pl.scatter` enhancements: multiplots, rugplot, linear and polynomial fits, density plots, etc.
- `pl.heatmap`: heatmap / clustermap of genes along time coordinate sorted by expression along dynamics.

Preprocessing:

- New attributes in pp.filter_genes: `min_shared_counts` and `min_shared_genes`.
- Added fast neighbor search method: Hierarchical Navigable Small World graphs (HNSW)


Version 0.1.14 :small:`Dec 7, 2018`
-----------------------------------
Plotting:

- New attriutes `arrow_length` and `arrow_size` for flexible adjustment of embedded velocities.
- `pl.velocity_graph`: Scatter plot of embedding with cell-to-cell transition connectivities.
- `pl.velocity_embedding_stream`: Streamplot visualization of velocities.
- Improve visualization of embedded single cell velocities (autosize, colors etc.)

Tools:

- `tl.cell_fate`: compute cell-specific terminal state likelihood
- New attribute `approx=True` in `tl.velocity_graph` to enable approximate graph computation by performing cosine correlations on PCA space.

Preprocessing:

- Automatically detect whether data is already preprocessed.


Version 0.1.11 :small:`Oct 27, 2018`
------------------------------------
Plotting:

- `settings.set_figure_params()`: adjust matplotlib defaults for beautified plots
- improved default point and arrow sizes; improved quiver autoscale
- enable direct plotting of

Tools:

- `tl.velocity_confidence`: Added two confidence measures 'velocity_confidence' and 'velocity_confidence_transition'.
- `tl.rank_velocity_genes`: Added functionality to rank genes for velocity characterizing groups using a t-test.
- New attribute `perc` in `tl.velocity` enables extreme quantile fit, e.g. set `perc=95`.
- New attribute `groups` in `tl.velocity` enables velocity estimation only on a subset of the data.
- Improved `tl.transition_matrix` by incorporating self-loops via `self_transitions=True`
  and state changes that have negative correlation with velocity (opposite direction) via `use_negative_cosines=True`

Utils:

- `utils.merge` to merge to AnnData objects such as already existing AnnData and newly generated Loom File.



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
