.. role:: small


Version 0.1.7 :small:`Sep 11, 2018`
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

- terminal_states: computes root and end points via eigenvalue decomposition


Version 0.1.5 :small:`Sep 4, 2018`
----------------------------------
- Support writing loom files
- Support both dense and sparse layers
- Plotting bugfixes
- Added pp.recipe_velocity()

Version 0.1.2 :small:`Aug 21, 2018`
-----------------------------------
First alpha release of scvelo.