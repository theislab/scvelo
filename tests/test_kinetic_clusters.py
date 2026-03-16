import pytest

import scanpy as sc

import scvelo as scv


def test_kinetic_clusters():
    """Test the full pipeline: Preprocessing -> Dynamics -> Clustering -> Plotting."""
    adata = scv.datasets.pancreas()
    sc.pp.filter_genes(adata, min_counts=20)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # small number of genes to keep the test fast
    sc.pp.highly_variable_genes(adata, n_top_genes=500)
    scv.pp.moments(adata, n_pcs=10, n_neighbors=10)

    # Calculate UMAP (Required for 'embedding' mode plotting)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    scv.tl.recover_dynamics(adata, max_iter=5, n_jobs=1)

    scv.tl.kinetic_clusters(adata, n_clusters=3)

    # Verify: Did it create the column?
    assert "kinetic_cluster" in adata.var.keys()
    # Verify: Are there actual clusters found (not just 'nan')?
    assert adata.var["kinetic_cluster"].value_counts().sum() > 0

    try:
        scv.pl.kinetic_clusters(adata, mode="pca", show=False)
    except Exception as e:
        pytest.fail(f"Plotting mode='pca' crashed: {e}")

    # This implicitly tests score_kinetic_clusters() as well
    try:
        scv.pl.kinetic_clusters(adata, mode="embedding", show=False)
    except Exception as e:
        pytest.fail(f"Plotting mode='embedding' crashed: {e}")


def test_kinetic_clusters_errors():
    """Test that the tool correctly raises errors when prerequisites are missing."""
    adata = scv.datasets.pancreas()
    sc.pp.filter_genes(adata, min_counts=20)

    # CASE 1: Run clustering without recover_dynamics (Should fail)
    # The tool needs 'fit_alpha' in adata.var
    with pytest.raises(ValueError, match="Kinetic parameters not found"):
        scv.tl.kinetic_clusters(adata)

    # CASE 2: Run plotting without running clustering first (Should fail)
    # manually add 'fit_alpha' to bypass the first check,
    # so we can test the specific check for 'kinetic_cluster' column.
    adata.var["fit_alpha"] = 1.0

    with pytest.raises(ValueError, match="Please run scv.tl.kinetic_clusters"):
        scv.pl.kinetic_clusters(adata, show=False)
