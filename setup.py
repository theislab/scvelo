from pathlib import Path

from setuptools import find_packages, setup


def read_requirements(req_path):
    """Read abstract requirements."""
    requirements = Path(req_path).read_text("utf-8").splitlines()
    return [r.strip() for r in requirements if not r.startswith("-")]


setup(
    name="scvelo",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    python_requires=">=3.6",
    install_requires=read_requirements("requirements.txt"),
    extras_require=dict(
        louvain=["python-igraph", "louvain"],
        hnswlib=["pybind11", "hnswlib"],
        dev=read_requirements("requirements-dev.txt"),
        docs=read_requirements("docs/requirements.txt"),
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
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
