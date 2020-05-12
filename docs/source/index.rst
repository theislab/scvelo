|PyPI| |PyPIDownloads| |Docs|

scVelo - RNA velocity generalized through dynamical modeling
============================================================

.. image:: https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png
   :width: 300px
   :align: left

**scVelo** is a scalable toolkit for RNA velocity analysis in single cells. |br|
The methods are based on our preprint
`Bergen et al. (2019) <https://doi.org/10.1101/820936>`_.

RNA velocity enables the recovery of directed dynamic information by leveraging splicing information.
scVelo generalizes the concept of RNA velocity (`La Manno et al., 2018 <https://doi.org/10.1038/s41586-018-0414-6>`_)
by relaxing previously made assumptions with a stochastic and a dynamical model that solves the full
transcriptional dynamics. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations.

scVelo is compatible with scanpy_ and hosts efficient implementations of all RNA velocity models.

scVelo's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^
- estimate RNA velocity to study cellular dynamics.
- identify putative driver genes and regimes of regulatory changes.
- infer a latent time to reconstruct the temporal sequence of transcriptomic events.
- estimate reaction rates of transcription, splicing and degradation.
- use statistical tests, e.g., to detect different kinetics regimes.

Reference
^^^^^^^^^
Bergen et al. (2019), *Generalizing RNA velocity to transient cell states through dynamical modeling*,
`biorxiv <https://doi.org/10.1101/820936>`_.

Support
^^^^^^^
Feel free to submit an `issue <https://github.com/theislab/scvelo/issues/new/choose>`_
or send us an `email <mailto:mail@scvelo.org>`_.
Your help to improve scVelo is highly appreciated.


.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:

   about
   installation
   api
   release_notes
   references

.. toctree::
   :caption: Tutorials
   :maxdepth: 1
   :hidden:

   Getting Started <https://scvelo.readthedocs.io/getting_started.html>
   RNA Velocity Basics <https://scvelo-notebooks.readthedocs.io/VelocityBasics.html>
   Dynamical Modeling <https://scvelo-notebooks.readthedocs.io/DynamicalModeling.html>
   Differential Kinetics <https://scvelo-notebooks.readthedocs.io/DifferentialKinetics.html>

.. toctree::
   :caption: Example Datasets
   :maxdepth: 1
   :hidden:

   Endocrine Pancreas <https://scvelo-notebooks.readthedocs.io/Pancreas.html>
   Dentate Gyrus <https://scvelo-notebooks.readthedocs.io/DentateGyrus.html>

.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
   :target: https://pypi.org/project/scvelo

.. |PyPIDownloads| image:: https://pepy.tech/badge/scvelo
   :target: https://pepy.tech/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |travis| image:: https://travis-ci.org/theislab/scvelo.svg?branch=master
   :target: https://travis-ci.org/theislab/scvelo

.. _scanpy: https://scanpy.readthedocs.io

.. |br| raw:: html

  <br/>