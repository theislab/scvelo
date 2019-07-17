|PyPI| |Docs| |travis|

scVelo - stochastic single cell RNA velocity
==========================================================

.. image:: https://drive.google.com/uc?export=view&id=1rcgHou-YFTJCKDR-Vd37zQ_AvLiaHLut
   :width: 600px
   :align: center

**scVelo** is a scalable toolkit for estimating and analyzing RNA velocities in single cells.

RNA velocity, the time derivative of mRNA abundance, enables you to infer directionality in your data by superimposing
splicing information. The main principles have been presented in
`La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`_,
and are based on a deterministic steady-state model of transcriptional dynamics.
scVelo uses a stochastic dynamical model that captures the complete molecular kinetics, thereby recovering the systematic
relationship of unspliced to spliced mRNA. That enables you to estimate velocities even in non-stationary populations.

It is compatible with scanpy_ (`Wolf et al., 2018 <https://doi.org/10.1186/s13059-017-1382-0>`_).
Making use of sparse implementation, iterative neighbors search and other techniques, it is remarkably efficient in
terms of memory and runtime without loss in accuracy and runs easily on your local machine (30k cells in a few minutes).

.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
    :target: https://pypi.org/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |travis| image:: https://travis-ci.org/theislab/scvelo.svg?branch=master
   :target: https://travis-ci.org/theislab/scvelo

.. _velocyto: http://velocyto.org/
.. _scanpy: https://github.com/theislab/scanpy
.. _notebooks: https://nbviewer.jupyter.org/github/theislab/scvelo_notebooks/tree/master/

Report issues and see the code on `GitHub <https://github.com/theislab/scvelo>`__.


.. toctree::
   :maxdepth: 1
   :hidden:

   getting_started
   Example <https://scvelo-notebooks.readthedocs.io/en/latest/DentateGyrus.html>
   api
   release_notes
   references
