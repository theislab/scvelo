[![PyPi][badge-pypi]][link-pypi]
[![PyPIDownloads][badge-pypidownloads]][link-pypidownloads]
[![CI][badge-ci]][link-ci]

[badge-pypi]: https://img.shields.io/pypi/v/scvelo.svg
[link-pypi]: https://pypi.org/project/scvelo
[badge-pypidownloads]: https://pepy.tech/badge/scvelo
[link-pypidownloads]: https://pepy.tech/project/scvelo
[badge-ci]: https://img.shields.io/github/actions/workflow/status/theislab/scvelo/ci.yml?branch=main
[link-ci]: https://github.com/theislab/scvelo/actions/workflows/ci.yml

# scVelo - RNA velocity generalized through dynamical modeling

<img src="https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png" width="400px" align="left">

**scVelo** is a scalable toolkit for RNA velocity analysis in single cells; RNA velocity
enables the recovery of directed dynamic information by leveraging splicing kinetics
<sup>[1](https://doi.org/10.1038/s41586-018-0414-6)</sup>. scVelo collects different
methods for inferring RNA velocity using an expectation-maximization framework
<sup>[2](https://doi.org/10.1038/s41587-020-0591-3)</sup>, deep generative modeling
<sup>[3](https://doi.org/10.1038/s41592-023-01994-w)</sup>,
or metabolically labeled transcripts<sup>[4](https://doi.org/10.1101/2023.07.19.549685)</sup>.

## scVelo's key applications

-   estimate RNA velocity to study cellular dynamics.
-   identify putative driver genes and regimes of regulatory changes.
-   infer a latent time to reconstruct the temporal sequence of transcriptomic events.
-   estimate reaction rates of transcription, splicing and degradation.
-   use statistical tests, e.g., to detect different kinetics regimes.

## Citing scVelo

If you include or rely on scVelo when publishing research, please adhere to the
following citation guide:

### EM and steady-state model

If you use the _EM_ (_dynamical_) or _steady-state model_, cite

```bibtex
@article{Bergen2020,
  title = {Generalizing RNA velocity to transient cell states through dynamical modeling},
  volume = {38},
  ISSN = {1546-1696},
  url = {http://dx.doi.org/10.1038/s41587-020-0591-3},
  DOI = {10.1038/s41587-020-0591-3},
  number = {12},
  journal = {Nature Biotechnology},
  publisher = {Springer Science and Business Media LLC},
  author = {Bergen, Volker and Lange, Marius and Peidli, Stefan and Wolf, F. Alexander and Theis, Fabian J.},
  year = {2020},
  month = aug,
  pages = {1408–1414}
}
```

### veloVI

If you use _veloVI_ (_VI model_), cite

```bibtex
@article{Gayoso2023,
  title = {Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells},
  ISSN = {1548-7105},
  url = {http://dx.doi.org/10.1038/s41592-023-01994-w},
  DOI = {10.1038/s41592-023-01994-w},
  journal = {Nature Methods},
  publisher = {Springer Science and Business Media LLC},
  author = {Gayoso, Adam and Weiler, Philipp and Lotfollahi, Mohammad and Klein, Dominik and Hong, Justin and Streets, Aaron and Theis, Fabian J. and Yosef, Nir},
  year = {2023},
  month = sep
}
```

### RNA velocity inference through metabolic labeling information

If you use the implemented method for estimating RNA velocity from metabolic labeling
information, cite

```bibtex
@article{Weiler2023,
  title = {Unified fate mapping in multiview single-cell data},
  url = {http://dx.doi.org/10.1101/2023.07.19.549685},
  DOI = {10.1101/2023.07.19.549685},
  publisher = {Cold Spring Harbor Laboratory},
  author = {Weiler, Philipp and Lange, Marius and Klein, Michal and Pe’er, Dana and Theis, Fabian J.},
  year = {2023},
  month = jul
}
```

## Support

Found a bug or would like to see a feature implemented? Feel free to submit an
[issue](https://github.com/theislab/scvelo/issues/new/choose).
Have a question or would like to start a new discussion? Head over to
[GitHub discussions](https://github.com/theislab/scvelo/discussions).
Your help to improve scVelo is highly appreciated.
For further information visit [scvelo.org](https://scvelo.org).
