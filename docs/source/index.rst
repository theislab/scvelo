scvelo – stochastic single cell RNA velocity
==========================================================

**scvelo** is a scalable toolkit for estimating and analyzing stochastic RNA velocities in single cells.

RNA velocity is the time derivative of mRNA abundance obtained by distinguishing unspliced (precursor) from spliced
(mature) mRNA, and serves as a predictive indicator for the future state of an individual cell. The main principles
of RNA velocity estimation have been presented in
velocyto_ (`La Manno et al., 2018 <https://doi.org/10.1038/s41586-018-0414-6>`_)
and are based on a deterministic model of transcriptional dynamics. scvelo uses a stochastic formulation and
incorporates intrinsic expression variability.

It is compatible with scanpy_ (`Wolf et al., 2018 <https://doi.org/10.1186/s13059-017-1382-0>`_).

Read the documentation_.

.. _velocyto: http://velocyto.org/
.. _scanpy: https://github.com/theislab/scanpy
.. _documentation: https://scvelo.readthedocs.io

Report issues and see the code on `GitHub <https://github.com/VolkerBergen/scRNAvelo>`__.

.. include:: release_notes.rst

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   usage_principles
   references