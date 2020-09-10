About scVelo
------------

Measuring gene activity in individual cells requires destroying these cells to read out their content, making it
challenging to study dynamic processes and to learn about cellular decision making. The introduction of RNA velocity by
`La Manno et al. (Nature, 2018) <https://doi.org/10.1038/s41586-018-0414-6>`_ has
enabled the recovery of directed dynamic information by leveraging the fact that newly
transcribed, unspliced pre-mRNAs and mature, spliced mRNAs can be distinguished in common single-cell RNA-seq protocols,
the former detectable by the presence of introns.
This concept of measuring not only gene activity, but also their changes in individual cells (RNA velocity),
has opened up new ways of studying cellular differentiation. The originally proposed framework obtains velocities as the deviation of the observed ratio of spliced and unspliced
mRNA from an inferred steady state. Errors in velocity estimates arise if the central assumptions of a common splicing
rate and the observation of the full splicing dynamics with steady-state mRNA levels are violated.

With scVelo, developed by `Bergen et al. (Nature Biotechnology, 2020) <https://doi.org/10.1038/s41587-020-0591-3>`_,
these restrictions are addressed by solving the full transcriptional dynamics of splicing kinetics using
a likelihood-based dynamical model. This generalizes RNA velocity to a wide variety of systems comprising transient
cell states, which are common in development and in response to perturbations.
Further, scVelo infers gene-specific rates of transcription, splicing and degradation, and recovers the latent time of the underlying
cellular processes. This latent time represents the cell’s internal clock and approximates the real time experienced by
cells as they differentiate, based only on its transcriptional dynamics.
Moreover, scVelo identifies regimes of regulatory changes such as stages of cell fate commitment and, therein,
systematically detects putative driver genes.


RNA velocity models
~~~~~~~~~~~~~~~~~~~
With RNA velocity, inference of directional trajectories is explored by connecting measurements to the underlying mRNA splicing kinetics:
Transcriptional induction for a particular gene results in an increase of (newly transcribed) precursor unspliced mRNAs
while, conversely, repression or absence of transcription results in a decrease of unspliced mRNAs.
Hence, by distinguishing unspliced from spliced mRNA, the change of mRNA abundance (RNA velocity) can be approximated.
The combination of velocities across mRNAs can then be used to estimate the future state of an individual cell.

RNA velocity estimation can currently be tackled with three existing approaches:

- steady-state / deterministic model (using steady-state residuals)
- stochastic model (using second-order moments),
- dynamical model (using a likelihood-based framework).

The **steady-state / deterministic model**, as being used in velocyto, estimates velocities as follows: Under the assumption
that transcriptional phases (induction and repression) last sufficiently long to reach a steady-state equilibrium
(active and inactive), velocities are quantified as the deviation of the observed ratio from its steady-state ratio.
The equilibrium mRNA levels are approximated with a linear regression on the presumed steady states in the lower and upper quantiles.
This simplification makes two fundamental assumptions: a common splicing rate across genes and steady-state mRNA levels to be
reflected in the data. It can lead to errors in velocity estimates and cellular states as the assumptions are often
violated, in particular when a population comprises multiple heterogeneous subpopulation dynamics.

The **stochastic model** aims to better capture the steady states. By treating transcription, splicing and degradation
as probabilistic events, the resulting Markov process is approximated by moment equations.
By including second-order moments, it exploits not only the balance of unspliced to spliced
mRNA levels but also their covariation. It has been demonstrated on the endocrine pancreas that
stochasticity adds valuable information, overall yielding higher consistency than the deterministic
model, while remaining as efficient in computation time.

The **dynamical model** (most powerful while computationally most expensive) solves the full dynamics of splicing kinetics
for each gene. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations,
as does not rely on the restrictions of a common splicing rate or steady states to be sampled.

The splicing dynamics

.. math::
   \begin{align}
   \frac{du(t)}{dt}=&~ \alpha_k(t) - \beta u(t),\\
   \frac{ds(t)}{dt}=&~ \beta u(t) - \gamma s(t),
   \end{align}

is solved in a likelihood-based expectation-maximization framework, by iteratively estimating the
parameters of reaction rates and latent cell-specific variables, i.e. transcriptional state *k* and cell-internal latent time *t*.

It thereby aims to learn the unspliced/spliced phase trajectory.
Four transcriptional states are modeled to account for all possible configurations of gene activity:
two dynamic transient states (induction and repression) and two steady states (active and inactive)
potentially reached after each dynamic transition.

In the expectation step, for a given model estimate of the unspliced/spliced phase trajectory,
a latent time is assigned to an observed mRNA value by minimizing its distance to the phase trajectory.
The transcriptional states are then assigned by associating a likelihood to the respective segments on the phase trajectory
(induction, repression, active and inactive steady states).
In the maximization step, the overall likelihood is then optimized by updating the parameters of reaction rates.

The model yields more consistent velocity estimates and better identification of transcriptional states.
It further enables the systematic identification of dynamics-driving genes in a likelihood-based way,
thereby finding the key drivers that govern cell fate transitions. Moreover, the dynamical model infers a universal
cell-internal latent time shared across genes that enables relating genes and identifying regimes of transcriptional changes.

For best results and above-described additional insights, we recommend using the dynamical model.
If runtime is of importance, the stochastic model is advised to be used as it very efficiently approximates the dynamical model,
taking few minutes on 30k cells. The dynamical yet can take up to one hour, however, enhancing efficiency is in progress.

See `Bergen et al. (2020) <https://doi.org/10.1038/s41587-020-0591-3>`_ for a detailed exposition of the methods.
