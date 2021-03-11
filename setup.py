from pathlib import Path

from setuptools import find_packages, setup

setup(
    name="scvelo",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    python_requires=">=3.6",
    install_requires=[
        library.strip()
        for library in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    extras_require=dict(
        louvain=["python-igraph", "louvain"],
        hnswlib=["pybind11", "hnswlib"],
        dev=[
            library.strip()
            for library in Path("requirements-dev.txt").read_text("utf-8").splitlines()
            if not library.startswith("-")
        ],
        docs=[r for r in Path("docs/requirements.txt").read_text("utf-8").splitlines()],
    ),
    packages=find_packages(),
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description="RNA velocity generalized through dynamical modeling",
    long_description=Path("pypi.rst").read_text("utf-8"),
    license="BSD",
    url="https://github.com/theislab/scvelo",
    download_url="https://github.com/theislab/scvelo",
    keywords=[
        "RNA",
        "velocity",
        "single cell",
        "transcriptomics",
        "stochastic",
        "dynamical",
    ],
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
