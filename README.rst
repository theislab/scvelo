|PyPI| |Docs| |travis|

scVelo - single-cell RNA velocity generalized to transient cell states
======================================================================

.. image:: https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png
   :width: 90px
   :align: left

**scVelo** is a scalable toolkit for estimating and analyzing RNA velocities in single cells using dynamical modeling.

The methods used herein are based on our preprint `Bergen et al. (2019) <https://doi.org/10.1101/820936>`_.

RNA velocity, the time derivative of mRNA abundance, enables you to infer directionality in your data by superimposing
splicing information. The main principles have been presented in
`La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`_,
and are based on a deterministic steady-state model of transcriptional dynamics.
scVelo provides two extensions as described in `Bergen et al. (2019) <https://doi.org/10.1101/820936>`_:
A stochastic model that incorporates second-order moments, and a dynamical model that captures the full splicing
kinetics. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations.

It is compatible with scanpy_ (`Wolf et al., 2018 <https://doi.org/10.1186/s13059-017-1382-0>`_).
Making use of sparse implementation, iterative neighbors search and other techniques, it is remarkably efficient in
terms of memory and runtime without loss in accuracy and runs easily on your local machine (30k cells in a few minutes).

Install from PyPI::

    pip install -U scvelo

See the the documentation at `<https://scvelo.org>`_ for details, which includes:

- **How scVelo estimates RNA velocity:** `<https://scvelo.org/about.html>`_
- **How to install and getting started:** `<https://scvelo.org/getting_started.html>`_
- **About the API:** `<https://scvelo.org/api.html>`_
- **Tutorials:** `<http://tutorials.scvelo.org/Pancreas.html>`_

Your feedback, in particular any issue you stumble upon, is highly appreciated and addressed to `feedback@scvelo.org <mailto:feedback@scvelo.org>`_.


.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
    :target: https://pypi.org/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |travis| image:: https://travis-ci.org/theislab/scvelo.svg?branch=master
   :target: https://travis-ci.org/theislab/scvelo

.. _scanpy: https://github.com/theislab/scanpy
.. _documentation: https://scvelo.org
