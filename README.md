[![PyPi][badge-pypi]][link-pypi]
[![PyPIDownloads][badge-pypidownloads]][link-pypidownloads]
[![CI][badge-ci]][link-ci]

[badge-pypi]: https://img.shields.io/pypi/v/scvelo.svg
[link-pypi]: https://pypi.org/project/scvelo
[badge-pypidownloads]: https://pepy.tech/badge/scvelo
[link-pypidownloads]: https://pepy.tech/project/scvelo
[badge-ci]: https://img.shields.io/github/actions/workflow/status/theislab/scvelo/ci.yml?branch=master
[link-ci]: https://github.com/theislab/scvelo/actions/workflows/ci.yml

# scVelo - RNA velocity generalized through dynamical modeling

<img src="https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png" width="400px" align="left">

**scVelo** is a scalable toolkit for RNA velocity analysis in single cells, based on
[Bergen et al. (Nature Biotech, 2020](https://doi.org/10.1038/s41587-020-0591-3).

RNA velocity enables the recovery of directed dynamic information by leveraging splicing kinetics.
scVelo generalizes the concept of RNA velocity
[La Manno et al., Nature, 2018](https://doi.org/10.1038/s41586-018-0414-6>)
by relaxing previously made assumptions with a stochastic and a dynamical model that solves the full
transcriptional dynamics. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations.

scVelo is compatible with scanpy\_ and hosts efficient implementations of all RNA velocity models.

## scVelo's key applications

-   estimate RNA velocity to study cellular dynamics.
-   identify putative driver genes and regimes of regulatory changes.
-   infer a latent time to reconstruct the temporal sequence of transcriptomic events.
-   estimate reaction rates of transcription, splicing and degradation.
-   use statistical tests, e.g., to detect different kinetics regimes.

scVelo has, for instance, recently been used to study immune response in COVID-19
patients and dynamic processes in human lung regeneration. Find out more in this list of
[application examples](https://scholar.google.com/scholar?cites=18195185735875895912).

## Latest news

-   Aug/2021: [Perspectives paper out in MSB](https://doi.org/10.15252/msb.202110282)
-   Feb/2021: scVelo goes multi-core
-   Dec/2020: Cover of [Nature Biotechnology](https://www.nature.com/nbt/volumes/38)
-   Nov/2020: Talk at [Single Cell Biology](https://coursesandconferences.wellcomegenomecampus.org/our-events/single-cell-biology-2020/)
-   Oct/2020: [Helmholtz Best Paper Award](https://twitter.com/ICBmunich/status/1318611467722199041)
-   Oct/2020: Map cell fates with [CellRank](https://cellrank.org)
-   Sep/2020: Talk at [Single Cell Omics](https://twitter.com/fabian_theis/status/1305621028056465412)
-   Aug/2020: [scVelo out in Nature Biotech](https://www.helmholtz-muenchen.de/en/aktuelles/latest-news/press-information-news/article/48658/index.html)

## References

[La Manno _et al._ (2018), RNA velocity of single cells, Nature](https://doi.org/10.1038/s41586-018-0414-6).

[Bergen _et al._ (2020), Generalizing RNA velocity to transient cell states through dynamical modeling, Nature Biotech](https://doi.org/10.1038/s41587-020-0591-3).

[Bergen _et al._ (2021), RNA velocity - current challenges and future perspectives, Molecular Systems Biology](https://doi.org/10.15252/msb.202110282).

## Support

Found a bug or would like to see a feature implemented? Feel free to submit an
[issue](https://github.com/theislab/scvelo/issues/new/choose).
Have a question or would like to start a new discussion? Head over to
[GitHub discussions](https://github.com/theislab/scvelo/discussions).
In either case, you can also always send us an [email](mailto:mail@scvelo.org).
Your help to improve scVelo is highly appreciated.
For further information visit [scvelo.org](https://scvelo.org).
