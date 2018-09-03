from setuptools import setup, find_packages
from pathlib import Path
import numpy as np
import versioneer

req_path = Path('requirements.txt')
requires = [r for r in req_path.read_text().strip().split('\n')]

with open('README.rst') as readme_f:
    readme = readme_f.read()

setup(
    name="scvelo",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=requires,
    python_requires='>=3.6',
    packages=find_packages(),
    include_dirs=[np.get_include()],
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='stochastic single cell RNA velocity ',
    long_description=readme,
    license='BSD',
    url="https://github.com/theislab/scvelo",
    download_url=f"https://github.com/theislab/scvelo",
    keywords=["RNAseq", "singlecell", "stochastic", "velocity", "transcriptomics"]
    )
