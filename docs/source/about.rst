About RNA velocity
------------------

RNA velocity enables you to infer directionality in your data by superimposing splicing information.

With the astounding finding that unspliced and spliced mRNA abundances can be distinguished in standard single cell
protocols, `La Manno et al., (2018) <https://doi.org/10.1038/s41586-018-0414-6>`_ introduced the concept of RNA velocity.
Inference of directional trajectories is explored by connecting measurements to the underlying mRNA splicing kinetics:
Transcriptional induction for a particular gene results in an increase of (newly transcribed) precursor unspliced mRNAs
while, conversely, repression or absence of transcription results in a decrease of unspliced mRNAs.
Hence, by distinguishing unspliced from mature spliced mRNA, the change of mRNA abundance, i.e. its time derivative,
denoted as RNA velocity, can be approximated. The combination of velocities across mRNAs can then be used to estimate
the future state of an individual cell.
The movement of the cells is visualized by projecting the velocities into a lower-dimensional embedding.

How scVelo estimates RNA velocities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RNA velocity estimation can currently be tackled with three existing approaches:

- steady-state / deterministic model (as being used in velocyto)
- stochastic model (using second-order moments),
- dynamical model (using a likelihood-based framework).

The **steady-state / deterministic model**, as being used in velocyto, estimates velocities as follows: Under the assumption
that transcriptional phases (induction and repression) endure long enough to reach a steady-state equilibrium
(active and inactive), velocities are quantified as how the actual observations deviate from the steady-state equilibrium.
The equilibrium mRNA levels are approximated with a linear regression on the presumed steady states in the lower and upper quantiles.
This simplification is obtained by assuming of a common splicing rate across genes and steady-state mRNA levels to be
reflected in the data. It can lead to errors in velocity estimates and cellular states as the assumptions are often
violated, in particular when a population comprises multiple heterogeneous subpopulation dynamics.

The **stochastic model** aims to better capture the steady states, yet relying on the same assumptions as the steady-state model.
It is obtained by treating transcription, splicing and degradation as probabilistic events,
thereby incorporating second-order moments. That is, steady state levels are
approximated not only from mRNA levels, but also from intrinsic expression variability.

The **dynamical model** (most powerful while computationally most expensive) solves the full dynamics of splicing kinetics
for each gene. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations,
as it does not rely on the restrictions of a common splicing rate or steady states to be sampled.
The splicing dynamics is solved in a likelihood-based expectation-maximization framework, by iteratively estimating the identifiable
parameters of reaction rates and latent cell-specific variables, i.e. transcriptional state and cell-internal latent time.
More precisely, we explicitly model four transcriptional states to account for all possible configurations
of gene activity: two dynamic transient states (induction and repression) and two steady states (active and inactive)
potentially reached after each dynamic transition. For each observation per state a model-optimal latent time is computed
to obtain a mapping onto the learned curve of unspliced/spliced dynamics. The cell-to-curve mapping then yields
likelihoods of respective state membership, and individual reaction rate parameters via maximizing the overall likelihood.

This yields more consistent velocity estimates and better identification of transcriptional states.
The model further enables to systemically identifies dynamics-driving genes in a likelihood-based way,
thereby finding the key drivers that govern cell fate transitions. Moreover, the dynamical model infers a universal
cell-internal latent time shared across genes that enables relating genes and identifying regimes of transcriptional changes (e.g. branching points).

For best results we obviously recommend using the superior dynamical model.
If runtime matters, we recommend using the stochastic model, which is designed to approximate the dynamical model.
The stochastic model takes less than a minute on 30k cells.
The dynamical yet can take up to one hour, however, enhancing efficiency is in progress.

Beyond RNA velocity
~~~~~~~~~~~~~~~~~~~
There are multiple extensions that can be easily explored, including terminal states (root and end points),
pseudotemporal ordering based on velocities, infer directionality in abstracted graph and many more. A more detailed description will follow very soon.
