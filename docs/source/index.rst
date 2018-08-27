scvelo – stochastic single cell RNA velocity
==========================================================

|PyPI| |Docs| |travis|

**scvelo** is a scalable toolkit for estimating and analyzing stochastic RNA velocities in single cells.

RNA velocity is the time derivative of mRNA abundance obtained by distinguishing unspliced (precursor) from spliced
(mature) mRNA, and serves as a predictive indicator for the future state of an individual cell. The main principles
of RNA velocity estimation have been presented in
velocyto_ (`La Manno et al., 2018 <https://doi.org/10.1038/s41586-018-0414-6>`_)
and are based on a deterministic model of transcriptional dynamics. scvelo uses a stochastic formulation and
incorporates intrinsic expression variability.

It is compatible with scanpy_ (`Wolf et al., 2018 <https://doi.org/10.1186/s13059-017-1382-0>`_). Making use of sparse
implementation, iterative neighbors search and other techniques, it is remarkably efficient in terms of memory and
runtime without loss in accuracy (<1GB and <1min for 30,000 cells on a MacBook Pro 2017 with 2.3 GHz i5).

.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
    :target: https://pypi.org/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |travis| image:: https://travis-ci.org/theislab/scvelo.svg?branch=master
   :target: https://travis-ci.org/theislab/scvelo

.. _velocyto: http://velocyto.org/
.. _scanpy: https://github.com/theislab/scanpy

Report issues and see the code on `GitHub <https://github.com/theislab/scvelo>`__.

.. include:: release_notes.rst

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   usage_principles
   api/index
   references