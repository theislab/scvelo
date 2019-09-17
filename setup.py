from setuptools import setup, find_packages
from pathlib import Path

import sys
if sys.version_info < (3, 6):
    sys.exit('scvelo requires Python >= 3.6')

setup(
    name="scvelo",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    python_requires='>=3.6',
    install_requires=['scanpy>=1.4.0', 'anndata>=0.6.18', 'numpy>=1.17.0', 'pandas>=0.23.0', 'matplotlib>=2.2',
                      'scipy>=1.0', 'scikit-learn>=0.19.1, != 0.21.0, != 0.21.1',  # exclude buggy versions
                      'joblib', 'umap-learn>=0.3.0', 'loompy>=2.0.12'],
    packages=find_packages(),
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='Stochastic RNA velocity for inferring single cell dynamics',
    long_description=Path('README.rst').read_text('utf-8'),
    license='BSD',
    url="https://github.com/theislab/scvelo",
    download_url=f"https://github.com/theislab/scvelo",
    keywords=["RNAseq", "singlecell", "stochastic", "velocity", "transcriptomics"],
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization']
    )
