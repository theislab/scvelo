|PyPI| |PyPIDownloads| |CI|

scVelo - RNA velocity generalized through dynamical modeling
============================================================

.. raw:: html

    <a href="https://scvelo.org">
    <img src="https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png" width="400px" align="left">
    </a>

**scVelo** is a scalable toolkit for RNA velocity analysis in single cells. |br|
The methods are based on
`Bergen et al. (Nature Biotech, 2020) <https://doi.org/10.1038/s41587-020-0591-3>`_.

RNA velocity enables the recovery of directed dynamic information by leveraging splicing information.
scVelo generalizes the concept of RNA velocity (`La Manno et al., Nature, 2018 <https://doi.org/10.1038/s41586-018-0414-6>`_)
by relaxing previously made assumptions with a stochastic and a dynamical model that solves the full
transcriptional dynamics. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations.
|br|

scVelo is compatible with scanpy_ and hosts efficient implementations of all RNA velocity models.

See `<https://scvelo.org>`_ for documentation and tutorials.

scVelo's key applications
-------------------------
- estimate RNA velocity to study cellular dynamics.
- identify putative driver genes and regimes of regulatory changes.
- infer a latent time to reconstruct the temporal sequence of transcriptomic events.
- estimate reaction rates of transcription, splicing and degradation.
- use statistical tests, e.g., to detect different kinetics regimes.

Upcoming talks
--------------
- `Sep 15: Temporal Single-Cell Analysis (SCOG) <https://twitter.com/fabian_theis/status/1305621028056465412>`_
- `Nov 12: Single Cell Biology (SCB) <https://coursesandconferences.wellcomegenomecampus.org/our-events/single-cell-biology-2020/>`_

Reference
---------
Bergen *et al.* (2020), Generalizing RNA velocity to transient cell states through dynamical modeling,
`Nature Biotechnology <https://doi.org/10.1038/s41587-020-0591-3>`_.
|dim|


Support
-------
Feel free to submit an `issue <https://github.com/theislab/scvelo/issues/new/choose>`_
or send us an `email <mailto:mail@scvelo.org>`_.
Your help to improve scVelo is highly appreciated.
As GitHub requests currently exceeds our limited time capacities, please excuse potential delays.
For urgent matters, you may schedule a zoom meeting `here <https://calendly.com/scvelo>`_.


.. |PyPI| image:: https://img.shields.io/pypi/v/scvelo.svg
   :target: https://pypi.org/project/scvelo

.. |PyPIDownloads| image:: https://pepy.tech/badge/scvelo
   :target: https://pepy.tech/project/scvelo

.. |Docs| image:: https://readthedocs.org/projects/scvelo/badge/?version=latest
   :target: https://scvelo.readthedocs.io

.. |CI| image:: https://img.shields.io/github/workflow/status/theislab/scvelo/CI/master
   :target: https://github.com/theislab/scvelo/actions?query=workflow%3ACI

.. _scanpy: https://scanpy.readthedocs.io

.. |br| raw:: html

  <br/>

.. |dim| raw:: html

   <span class="__dimensions_badge_embed__" data-doi="10.1101/820936" data-style="small_rectangle"></span>
   <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
