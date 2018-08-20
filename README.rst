scvelo – single cell dynamics with stochastic RNA velocity
==========================================================

**scvelo** is a scalable toolkit for estimating and analyzing stochastic RNA velocities.

RNA velocity is the time derivative of mRNA abundance obtained by distinguishing unspliced (precursor) from spliced
(mature) mRNA, and serves as a predictive indicator for the future state of an individual cell. The main principles
of RNA velocity estimation have been presented in velocyto_ (La Manno et al., 2017) and are based on a deterministic
model of transcriptional dynamics. scvelo uses a stochastic formulation and incorporates intrinsic expression variability.

It is compatible with scanpy_ (Wolf et al., 2017). Making use of sparse implementation, multiprocessing,
iterative neighbor search and closed-form solutions, it is efficient in terms of memory and runtime (< 1 min for 20k cells).

Read the documentation_.

.. _velocyto: http://velocyto.org/
.. _scanpy: https://github.com/theislab/scanpy
.. _documentation: https://scvelo.readthedocs.io