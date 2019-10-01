|PyPI| |Docs| |travis|

scVelo - RNA velocity using dynamical modeling
==============================================

.. image:: https://user-images.githubusercontent.com/31883718/65906280-8cfc0d00-e3c2-11e9-94ee-bb74d3da15e2.png
   :width: 90px
   :align: left

**scVelo** is a scalable toolkit for estimating and analyzing RNA velocities in single cells using dynamical modeling.

RNA velocity, the time derivative of mRNA abundance, enables you to infer directionality in your data by superimposing
splicing information. The main principles have been presented in
`La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`_,
and are based on a deterministic steady-state model of transcriptional dynamics.
scVelo provides two extensions: A stochastic model that incorporates second-order moments,
and a dynamical model that captures the full splicing kinetics. It thereby adapts RNA velocity to widely varying
specifications such as non-stationary populations.

It is compatible with scanpy_ (`Wolf et al., 2018 <https://doi.org/10.1186/s13059-017-1382-0>`_).
Making use of sparse implementation, iterative neighbors search and other techniques, it is remarkably efficient in
terms of memory and runtime without loss in accuracy and runs easily on your local machine (30k cells in a few minutes).

Install from PyPI::

    pip install -U scvelo

See the documentation_ for details, which includes:

- `How scVelo estimates RNA velocity <https://scvelo.readthedocs.io/en/latest/about.html>`_
- `How to install and getting started <https://scvelo.readthedocs.io/en/latest/getting_started.html>`_
- `About the API <https://scvelo.readthedocs.io/en/latest/api.html>`_

Your feedback, in particular any issue you stumble upon, is highly appreciated and addressed to `feedback@scvelo.de <mailto:feedback@scvelo.de>`_.


.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
    :target: https://pypi.org/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |travis| image:: https://travis-ci.org/theislab/scvelo.svg?branch=master
   :target: https://travis-ci.org/theislab/scvelo

.. _scanpy: https://github.com/theislab/scanpy
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _documentation: https://scvelo.readthedocs.io
.. _`velocyto pipeline`: http://velocyto.org/velocyto.py/tutorial/cli.html
.. _`kallisto pipeline`: https://pachterlab.github.io/kallisto/about