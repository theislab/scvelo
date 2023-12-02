from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

from scanpy import read

from scvelo.core import cleanup
from scvelo.read_load import load

url_datadir = "https://github.com/theislab/scvelo_notebooks/raw/master/"


def bonemarrow(
    file_path: Optional[
        Union[str, Path]
    ] = "data/BoneMarrow/human_cd34_bone_marrow.h5ad"
):
    """Human bone marrow.

    Data from `Setty et al. (2019) <https://doi.org/10.1038/s41587-019-0068-4>`__.

    The bone marrow is the primary site of new blood cell production or haematopoiesis.
    It is composed of hematopoietic cells, marrow adipose tissue, and supportive stromal
    cells.

    This dataset served to detect important landmarks of hematopoietic differentiation,
    to identify key transcription factors that drive lineage fate choice and to closely
    track when cells lose plasticity.

    .. image:: https://user-images.githubusercontent.com/31883718/118402252-68bd7480-b669-11eb-9ef3-5f992b74a2d3.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "https://ndownloader.figshare.com/files/27686835"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def dentategyrus(file_path: Optional[Union[str, Path]] = None, adjusted=True):
    """Dentate Gyrus neurogenesis.

    Data from `Hochgerner et al. (2018) <https://doi.org/10.1038/s41593-017-0056-2>`__.

    Dentate gyrus (DG) is part of the hippocampus involved in learning, episodic memory
    formation and spatial coding. The experiment from the developing DG comprises two
    time points (P12 and P35) measured using droplet-based scRNA-seq
    (10x Genomics Chromium).

    The dominating structure is the granule cell lineage, in which neuroblasts develop
    into granule cells. Simultaneously, the remaining population forms distinct cell
    types that are fully differentiated (e.g. Cajal-Retzius cells) or cell types that
    form a sub-lineage (e.g. GABA cells).

    .. image:: https://user-images.githubusercontent.com/31883718/79433223-255b8700-7fcd-11ea-8ecf-3dc9eb1a6159.png
       :width: 600px

    Arguments:
    ---------
    file_path
        Path where to save dataset and read it from.

    Returns
    -------
    Returns `adata` object
    """
    if file_path is None and adjusted:
        file_path = "data/DentateGyrus/10X43_1.h5ad"
    elif file_path is None:
        file_path = "data/DentateGyrus/10X43_1.loom"

    if adjusted:
        url = f"{url_datadir}data/DentateGyrus/10X43_1.h5ad"
        adata = read(file_path, backup_url=url, sparse=True, cache=True)
    else:
        url = "http://pklab.med.harvard.edu/velocyto/DG1/10X43_1.loom"
        adata = read(file_path, backup_url=url, cleanup=True, sparse=True, cache=True)
        cleanup(adata, clean="all", keep={"spliced", "unspliced", "ambiguous"})

        url_louvain = f"{url_datadir}data/DentateGyrus/DG_clusters.npy"
        url_umap = f"{url_datadir}data/DentateGyrus/DG_umap.npy"

        adata.obs["clusters"] = load(
            "./data/DentateGyrus/DG_clusters.npy", url_louvain, allow_pickle=True
        )
        adata.obsm["X_umap"] = load(
            "./data/DentateGyrus/DG_umap.npy", url_umap, allow_pickle=True
        )

        adata.obs["clusters"] = pd.Categorical(adata.obs["clusters"])

    return adata


def dentategyrus_lamanno(
    file_path: Optional[Union[str, Path]] = "data/DentateGyrus/DentateGyrus.loom"
):
    """Dentate Gyrus neurogenesis.

    From `La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`__.

    The experiment from the developing mouse hippocampus comprises two time points
    (P0 and P5) and reveals the complex manifold with multiple branching lineages
    towards astrocytes, oligodendrocyte precursors (OPCs), granule neurons and pyramidal
    neurons.

    .. image:: https://user-images.githubusercontent.com/31883718/118401264-49bce380-b665-11eb-8678-e7570ede13d6.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "http://pklab.med.harvard.edu/velocyto/DentateGyrus/DentateGyrus.loom"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    adata.obsm["X_tsne"] = np.column_stack([adata.obs["TSNE1"], adata.obs["TSNE2"]])
    adata.obs["clusters"] = adata.obs["ClusterName"]
    cleanup(adata, clean="obs", keep=["Age", "clusters"])

    adata.uns["clusters_colors"] = {
        "RadialGlia": [0.95, 0.6, 0.1],
        "RadialGlia2": [0.85, 0.3, 0.1],
        "ImmAstro": [0.8, 0.02, 0.1],
        "GlialProg": [0.81, 0.43, 0.72352941],
        "OPC": [0.61, 0.13, 0.72352941],
        "nIPC": [0.9, 0.8, 0.3],
        "Nbl1": [0.7, 0.82, 0.6],
        "Nbl2": [0.448, 0.85490196, 0.95098039],
        "ImmGranule1": [0.35, 0.4, 0.82],
        "ImmGranule2": [0.23, 0.3, 0.7],
        "Granule": [0.05, 0.11, 0.51],
        "CA": [0.2, 0.53, 0.71],
        "CA1-Sub": [0.1, 0.45, 0.3],
        "CA2-3-4": [0.3, 0.35, 0.5],
    }
    return adata


def forebrain(file_path: Union[str, Path] = "data/ForebrainGlut/hgForebrainGlut.loom"):
    """Developing human forebrain.

    From `La Manno et al. (2018) <https://doi.org/10.1038/s41586-018-0414-6>`__.

    Forebrain tissue of a human week 10 embryo, focusing on glutamatergic neuronal
    lineage, obtained from elective routine abortions (10 weeks post-conception).

    Arguments:
    ---------
    file_path
        Path where to save dataset and read it from.

    Returns
    -------
    Returns `adata` object
    """
    url = "http://pklab.med.harvard.edu/velocyto/hgForebrainGlut/hgForebrainGlut.loom"
    adata = read(file_path, backup_url=url, cleanup=True, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def gastrulation(
    file_path: Optional[Union[str, Path]] = "data/Gastrulation/gastrulation.h5ad"
):
    """Mouse gastrulation.

    Data from `Pijuan-Sala et al. (2019) <https://doi.org/10.1038/s41586-019-0933-9>`__.

    Gastrulation represents a key developmental event during which embryonic pluripotent
    cells diversify into lineage-specific precursors that will generate the adult
    organism.

    This data contains the erythrocyte lineage from Pijuan-Sala et al. (2019).
    The experiment reveals the molecular map of mouse gastrulation and early
    organogenesis. It comprises transcriptional profiles of 116,312 single cells from
    mouse embryos collected at nine sequential time points ranging from 6.5 to 8.5 days
    post-fertilization. It served to explore the complex events involved in the
    convergence of visceral and primitive streak-derived endoderm.

    .. image:: https://user-images.githubusercontent.com/31883718/130636066-3bae153e-1626-4d11-8f38-6efab5b81c1c.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "https://ndownloader.figshare.com/files/28095525"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def gastrulation_e75(
    file_path: Optional[Union[str, Path]] = "data/Gastrulation/gastrulation_e75.h5ad"
):
    """Mouse gastrulation subset to E7.5.

    Data from `Pijuan-Sala et al. (2019) <https://doi.org/10.1038/s41586-019-0933-9>`__.

    Gastrulation represents a key developmental event during which embryonic pluripotent
    cells diversify into lineage-specific precursors that will generate the adult
    organism.

    .. image:: https://user-images.githubusercontent.com/31883718/130636292-7f2a599b-ded4-4616-99d7-604d2f324531.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "https://ndownloader.figshare.com/files/30439878"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def gastrulation_erythroid(
    file_path: Optional[Union[str, Path]] = "data/Gastrulation/erythroid_lineage.h5ad"
):
    """Mouse gastrulation subset to erythroid lineage.

    Data from `Pijuan-Sala et al. (2019) <https://doi.org/10.1038/s41586-019-0933-9>`__.

    Gastrulation represents a key developmental event during which embryonic pluripotent
    cells diversify into lineage-specific precursors that will generate the adult
    organism.

    .. image:: https://user-images.githubusercontent.com/31883718/118402002-40814600-b668-11eb-8bfc-dbece2b2b34e.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "https://ndownloader.figshare.com/files/27686871"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def pancreas(file_path: Union[str, Path] = "data/Pancreas/endocrinogenesis_day15.h5ad"):
    """Pancreatic endocrinogenesis.

    Data from `Bastidas-Ponce et al. (2019) <https://doi.org/10.1242/dev.173849>`__.

    Pancreatic epithelial and Ngn3-Venus fusion (NVF) cells during secondary transition
    with transcriptome profiles sampled from embryonic day 15.5.

    Endocrine cells are derived from endocrine progenitors located in the pancreatic
    epithelium. Endocrine commitment terminates in four major fates: glucagon- producing
    α-cells, insulin-producing β-cells, somatostatin-producing δ-cells and
    ghrelin-producing ε-cells.

    .. image:: https://user-images.githubusercontent.com/31883718/67709134-a0989480-f9bd-11e9-8ae6-f6391f5d95a0.png
       :width: 600px

    Arguments:
    ---------
    file_path
        Path where to save dataset and read it from.

    Returns
    -------
    Returns `adata` object
    """
    url = f"{url_datadir}data/Pancreas/endocrinogenesis_day15.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


def pbmc68k(file_path: Optional[Union[str, Path]] = "data/PBMC/pbmc68k.h5ad"):
    """Peripheral blood mononuclear cells.

    Data from `Zheng et al. (2017) <https://doi.org/10.1038/ncomms14049>`__.

    This experiment contains 68k peripheral blood mononuclear cells (PBMC) measured
    using 10X.

    PBMCs are a diverse mixture of highly specialized immune cells.
    They originate from hematopoietic stem cells (HSCs) that reside in the bone marrow
    and give rise to all blood cells of the immune system (hematopoiesis).
    HSCs give rise to myeloid (monocytes, macrophages, granulocytes, megakaryocytes,
    dendritic cells, erythrocytes) and lymphoid (T cells, B cells, NK cells) lineages.

    .. image:: https://user-images.githubusercontent.com/31883718/118402351-e1243580-b669-11eb-8256-4a49c299da3d.png
       :width: 600px

    Returns
    -------
    Returns `adata` object
    """
    url = "https://ndownloader.figshare.com/files/27686886"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata


# TODO: Remove function and add subsetting functionality for each dataset
def toy_data(
    file_path: Union[str, Path] = "data/DentateGyrus/10X43_1.h5ad", n_obs=None
):
    """Randomly sampled from the Dentate Gyrus dataset.

    Arguments:
    ---------
    file_path
        Path where to save dataset and read it from.
    n_obs: `int` (default: `None`)
        Size of the sampled dataset

    Returns
    -------
    Returns `adata` object
    """
    adata_dg = dentategyrus(file_path=file_path)

    if n_obs is not None:
        indices = np.random.choice(adata_dg.n_obs, n_obs)
        adata = adata_dg[indices]
    else:
        adata = adata_dg
    adata.obs_names_make_unique()
    return adata.copy()
