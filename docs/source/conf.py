# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import logging
import sys
from datetime import datetime
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import scvelo

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

logger = logging.getLogger(__name__)
# -- Project information -----------------------------------------------------

project = "scVelo"
author = "Volker Bergen, Philipp Weiler"
version = scvelo.__version__.replace(".dirty", "")
copyright = f"{datetime.now():%Y}, Theislab"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinxcontrib.bibtex",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.githubpages",
    "edit_on_github",
    "sphinx_autodoc_typehints",
    "nbsphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/latest/", None),
    "cellrank": ("https://cellrank.readthedocs.io/en/latest/", None),
}
master_doc = "index"
pygments_style = "sphinx"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

nitpicky = True

# bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

# -- Options for HTML output ----------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {"navigation_depth": 1, "titles_only": True}
github_repo = "scvelo"
github_nb_repo = "scvelo_notebooks"
html_static_path = ["_static"]

# -- Basic notebooks and those stored under /vignettes and /perspectives --

notebooks_url = "https://github.com/theislab/scvelo_notebooks/raw/master/"
notebooks = []
notebook = [
    "VelocityBasics.ipynb",
    "DynamicalModeling.ipynb",
    "DifferentialKinetics.ipynb",
]
notebooks.extend(notebook)

notebook = [
    "Pancreas.ipynb",
    "DentateGyrus.ipynb",
    "NatureBiotechCover.ipynb",
    "Fig1_concept.ipynb",
    "Fig2_dentategyrus.ipynb",
    "Fig3_pancreas.ipynb",
    "FigS9_runtime.ipynb",
    "FigSuppl.ipynb",
]
notebooks.extend([f"vignettes/{nb}" for nb in notebook])

notebook = ["Perspectives.ipynb", "Perspectives_parameters.ipynb"]
notebooks.extend([f"perspectives/{nb}" for nb in notebook])

# -- Retrieve all notebooks --

for nb in notebooks:
    url = notebooks_url + nb
    try:
        urlretrieve(url, nb)
    except URLError as e:
        logger.error(f"Unable to retrieve notebook: `{url}`. Reason: `{e}`")
