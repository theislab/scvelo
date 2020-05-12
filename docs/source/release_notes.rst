.. role:: small
.. role:: smaller

Release Notes
=============


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

- enhancements in `pl.scatter`: linear and polynomical fits, gradient coloring
- `pl.proportions`: Pie and bar chart of spliced/unspliced proprtions.
- `GridSpec`: multiplot environment.


Version 0.1.20 :small:`Sep 5, 2019`
-----------------------------------
Tools:

- `tl.recover_dynamics`: introduced a dynamical model inferring the full splicing kinetics, thereby identifying all kinetic rates of transcription, splicing and degradation.
- `tl.recover_latent_time`: infers a shared latent time across all genes based on the learned splicing dynamics.

Plotting:

- enhancements in `pl.scatter`: multiplots, rugplot, linear and polynomial fits, density plots, etc.
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