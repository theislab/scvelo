from setuptools import setup, find_packages
import numpy as np
from version import __version__

setup(
    name="scvelo",
    version=__version__,
    install_requires=['scanpy>=1.3.3', 'anndata>=0.6.18', 'loompy>=2.0.12', 'numpy>=1.14', 'pandas>=0.23.0',
                      'scipy>=1.0', 'scikit-learn>=0.19.1', 'matplotlib>=2.2', 'setuptools'],
    packages=find_packages(),
    include_dirs=[np.get_include()],
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='Stochastic RNA velocity for inferring single cell dynamics',
    license='BSD',
    url="https://github.com/theislab/scvelo",
    download_url=f"https://github.com/theislab/scvelo",
    keywords=["RNAseq", "singlecell", "stochastic", "velocity", "transcriptomics"]
    )
