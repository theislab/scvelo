from typing import Optional

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from anndata import AnnData

from scvelo import logging as logg


def kinetic_clusters(
    adata: AnnData, n_clusters: int = 4, copy: bool = False
) -> Optional[AnnData]:
    """
    Cluster genes based on their kinetic parameters.

    Identifies functional gene groups (e.g. 'Fast Response', 'Transient',
    'Accumulating') based on the differential equations learned by recover_dynamics.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n_clusters
        Number of kinetic regimes to find.
    copy
        Return a copy instead of writing to ``adata``.

    Returns
    -------
    Returns or updates ``adata`` with the column ``kinetic_cluster`` in ``adata.var``.
    """
    logg.info("clustering genes by kinetics", r=True)

    # 1. Validation
    if "fit_alpha" not in adata.var.keys():
        raise ValueError(
            "Kinetic parameters not found. Run scv.tl.recover_dynamics(adata) first."
        )

    # 2. Extract Data (Filter valid fits)
    valid_mask = (adata.var["fit_r2"] > 0) & (adata.var["fit_alpha"].notnull())

    if valid_mask.sum() < n_clusters:
        raise ValueError(
            f"Not enough fitted genes ({valid_mask.sum()}) to perform clustering."
        )

    features = adata.var.loc[valid_mask, ["fit_alpha", "fit_beta", "fit_gamma"]]

    # 3. Log Transform & Scale
    X = np.log1p(features.values)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # 4. Clustering (KMeans)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)

    # 5. Save Results
    col_name = "kinetic_cluster"
    adata.var[col_name] = "nan"
    adata.var.loc[valid_mask, col_name] = clusters.astype(str)

    logg.info("    finished", time=True, end=" ")
    logg.hint(f"added \n    {col_name!r}, clusters of kinetic parameters (adata.var)")

    return adata if copy else None


def score_kinetic_clusters(adata: AnnData, copy: bool = False) -> Optional[AnnData]:
    """
    Score cells based on kinetic cluster activity.

    Calculates a module score for each kinetic cluster across all cells.
    Allows visualization of where specific kinetic regimes are active.

    Parameters
    ----------
    adata
        Annotated data matrix.
    copy
        Return a copy instead of writing to ``adata``.

    Returns
    -------
    Returns or updates ``adata`` with columns ``Kinetic_Cluster_0``, etc. in ``adata.obs``.
    """
    import scanpy as sc

    if "kinetic_cluster" not in adata.var.keys():
        raise ValueError("Please run scv.tl.kinetic_clusters(adata) first.")

    clusters = adata.var["kinetic_cluster"].unique()
    # Filter out 'nan' string
    clusters = [c for c in clusters if str(c) != "nan"]

    logg.info("scoring kinetic clusters", r=True)

    for c in clusters:
        gene_list = adata.var[adata.var["kinetic_cluster"] == c].index.tolist()
        if len(gene_list) > 0:
            score_name = f"Kinetic_Cluster_{c}"
            sc.tl.score_genes(adata, gene_list=gene_list, score_name=score_name)

    logg.info("    finished", time=True, end=" ")
    logg.hint("added scores for each kinetic cluster to adata.obs")

    return adata if copy else None
