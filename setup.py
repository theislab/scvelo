from setuptools import setup, find_packages
from pathlib import Path

setup(
    name="scvelo",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    python_requires='>=3.6',
    install_requires=['scanpy>=1.4',
                      'anndata>=0.7',
                      'numpy>=1.17',
                      'scipy>=1.3',
                      'pandas>=0.23',
                      'scikit-learn>=0.21.2',
                      'matplotlib>=3.1.2',
                      'umap-learn>=0.3.10',
                      'loompy>=2.0.12'
                      ],
    packages=find_packages(),
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='RNA velocity generalized through dynamical modeling',
    long_description=Path('README.rst').read_text('utf-8'),
    license='BSD',
    url="https://github.com/theislab/scvelo",
    download_url="https://github.com/theislab/scvelo",
    keywords=["RNAseq", "single cell", "stochastic", "dynamical", "velocity", "transcriptomics"],
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
