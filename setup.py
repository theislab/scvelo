from setuptools import setup, find_packages
from pathlib import Path
import numpy as np
import versioneer

req_path = Path('requirements.txt')
if not req_path.is_file():
    req_path = Path('scvelo.egg-info') / req_path
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

setup(
    name="scvelo",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=requires,
    packages=find_packages(),
    include_dirs=[np.get_include()],
    author="Volker Bergen",
    author_email="volker.bergen@helmholtz-muenchen.de",
    description='Stochastic RNA velocity for inferring single cell dynamics',
    license='BSD',
    url="https://github.com/VolkerBergen/scvelo",
    download_url=f"https://github.com/VolkerBergen/scvelo",
    keywords=["RNAseq", "singlecell", "stochastic", "velocity", "transcriptomics"]
    )
