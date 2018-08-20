from setuptools import setup, find_packages
import numpy as np
import versioneer

setup(
    name="scvelo",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=['anndata',
                      'scanpy',
                      'pandas',
                      'numpy',
                      'scipy',
                      'numba',
                      'matplotlib',
                      'scikit-learn',
                      'h5py',
                      'loompy'],
    packages=find_packages(),
    include_dirs=[np.get_include()],
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='Stochastic RNA velocity for inferring single cell dynamics',
    license='BSD',
    url="https://github.com/VolkerBergen/scvelo",
    download_url=f"https://github.com/VolkerBergen/scvelo",
    keywords=["RNAseq", "singlecell", "stochastic", "velocity", "transcriptomics"],
    )
