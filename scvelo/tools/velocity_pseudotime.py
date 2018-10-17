import numpy as np


def principal_curve(data, basis='pca', n_comps=4, clusters_list=None, copy=False):
    """Computes the principal curve
    
    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'pca'`)
        Basis to use for computing the principal curve.
    n_comps: `int` (default: 4)
        Number of pricipal components to be used.
    copy: `bool`, (default: `False`)
        Return a copy instead of writing to adata.
    
    Returns
    -------
    Returns or updates `adata` with the attributes
    principal_curve: `.uns`
        dictionary containing `projections`, `ixsort` and `arclength`
    """
    adata = data.copy() if copy else data
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr

    if clusters_list is not None:
        cell_subset = np.array([label in clusters_list for label in adata.obs['clusters']])
        X_emb = adata[cell_subset].obsm['X_' + basis][:, :n_comps]
    else:
        cell_subset = None
        X_emb = adata.obsm['X_' + basis][:, :n_comps]

    n_obs, n_dim = X_emb.shape

    # convert array to R matrix
    xvec = robjects.FloatVector(X_emb.T.reshape((X_emb.size)))
    X_R = robjects.r.matrix(xvec, nrow=n_obs, ncol=n_dim)

    fit = importr("princurve").principal_curve(X_R)

    adata.uns['principal_curve'] = dict()
    adata.uns['principal_curve']['ixsort'] = ixsort = np.array(fit[1])-1
    adata.uns['principal_curve']['projections'] = np.array(fit[0])[ixsort]
    adata.uns['principal_curve']['arclength'] = np.array(fit[2])
    adata.uns['principal_curve']['cell_subset'] = cell_subset

    return adata if copy else None


# def pseudotime(adata, iroot=229, clusters_list=None, copy=False):
#     # iroot = adata.obsm['X_umap'][adata.obs['clusters']=='neuroblast 1'][:,1].argmax()
#     root_cell = adata.obs_names[iroot]
#     if clusters_list is not None:  # ['granule mature', 'neuroblast 1', 'neuroblast 2', 'granule immature']
#         cell_subset = np.array([label in clusters_list for label in adata.obs['clusters']])
#         adata_subset = adata.copy()
#         adata_subset._inplace_subset_obs(cell_subset)
#         adata_subset.uns['iroot'] = np.where(adata_subset.obs_names == root_cell)[0][0]
#     dpt(adata_subset, n_branchings=0)
#
#     adata.obs['dpt_pseudotime'] = np.zeros(adata.n_obs)
#     adata.obs['dpt_pseudotime'][cell_subset] = adata_subset.obs['dpt_pseudotime']
#     adata.obs['dpt_pseudotime'][~cell_subset] = -.5
#
#     return adata if copy else None
