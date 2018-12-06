from ..preprocessing import filter_and_normalize, moments
from . import velocity, velocity_graph, velocity_embedding


def run_all(data, basis=None, mode=None, min_counts=30, min_counts_u=20, n_top_genes=3000, n_pcs=30, n_neighbors=30,
            copy=False):
    from time import time
    start = time()

    adata = data.copy() if copy else data
    filter_and_normalize(adata, min_counts=min_counts, min_counts_u=min_counts_u, n_top_genes=n_top_genes)
    print("Number of genes to be used:", adata.X.shape[1])
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    velocity(adata, mode=mode)
    print("Number of genes to be used:", adata.var.velocity_genes.sum())

    velocity_graph(adata)
    if basis is not None: velocity_embedding(adata, basis=basis)

    print('Total time (seconds): ' + str(round(time() - start, 2)))

    return adata if copy else None


def convert_to_adata(vlm, basis=None):
    from collections import OrderedDict
    from .. import AnnData
    X = vlm.S_norm.T if hasattr(vlm, 'S_norm') else vlm.S_sz.T if hasattr(vlm, 'S_sz') else vlm.S.T  # .sparse().T.tocsr()

    layers = OrderedDict()
    layers['spliced'] = vlm.S_sz.T if hasattr(vlm, 'S_sz') else vlm.S.T
    layers['unspliced'] = vlm.U_sz.T if hasattr(vlm, 'U_sz') else vlm.U.T
    if hasattr(vlm, 'A') and (vlm.A.T.shape == layers['spliced'].shape): layers['ambiguous'] = vlm.A.T

    if hasattr(vlm, 'velocity'): layers['velocity'] = vlm.velocity.T
    if hasattr(vlm, 'Sx_sz'): layers['Ms'] = vlm.Sx_sz.T
    if hasattr(vlm, 'Ux_sz'): layers['Mu'] = vlm.Ux_sz.T

    obs = dict(vlm.ca)
    if 'CellID' in obs.keys(): obs['obs_names'] = obs.pop('CellID')

    var = dict(vlm.ra)
    if 'Gene' in var.keys(): var['var_names'] = var.pop('Gene')
    if hasattr(vlm, 'q'): var['velocity_offset'] = vlm.q
    if hasattr(vlm, 'gammas'): var['velocity_gamma'] = vlm.gammas
    if hasattr(vlm, 'R2'): var['velocity_r2'] = vlm.R2

    adata = AnnData(X, obs=obs, var=var, layers=layers)

    if hasattr(vlm, 'corrcoef'): adata.uns['velocity_graph'] = vlm.corrcoef

    if basis is not None:
        if hasattr(vlm, 'ts'): adata.obsm['X_' + basis] = vlm.ts
        if hasattr(vlm, 'delta_embedding'): adata.obsm['velocity_' + basis] = vlm.delta_embedding

    return adata


def convert_to_loom(adata, basis=None):
    from scipy.sparse import issparse
    import numpy as np
    import velocyto

    class VelocytoLoom(velocyto.VelocytoLoom):
        def __init__(self, adata, basis=None):
            self.S = adata.layers['spliced'].T
            self.U = adata.layers['unspliced'].T
            self.S = self.S.A if issparse(self.S) else self.S
            self.U = self.U.A if issparse(self.U) else self.U

            if 'ambiguous' in adata.layers.keys():
                self.A = adata.layers['ambiguous'].T
                if issparse(self.A): self.A = self.A.A

            self.initial_cell_size = self.S.sum(0)
            self.initial_Ucell_size = self.U.sum(0)

            self.ca = dict()
            self.ra = dict()

            self.ca['CellID'] = np.array(adata.obs_names)
            self.ra['Gene'] = np.array(adata.var_names)

            for key in adata.obs.keys():
                self.ca[key] = np.array(adata.obs[key].values)
            for key in adata.var.keys():
                self.ra[key] = np.array(adata.var[key].values)
            if isinstance(basis, str) and 'X_' + basis in adata.obsm.keys():
                if 'X_' + basis in adata.obsm.keys(): self.ts = adata.obsm['X_' + basis]
                if 'velocity_' + basis in adata.obsm.keys(): self.delta_embedding = adata.obsm['velocity_' + basis]

            if 'clusters' in self.ca:
                self.set_clusters(self.ca['clusters'])
            elif 'louvain' in self.ca:
                self.set_clusters(self.ca['louvain'])

        def filter_and_normalize(self, min_counts=3, min_counts_u=3, min_cells=0, min_cells_u=0,
                                 n_top_genes=None, max_expr_avg=20):
            # counterpart to scv.pp.filter_and_normalize()
            self.score_detection_levels(min_expr_counts=min_counts, min_expr_counts_U=0,
                                        min_cells_express=min_cells, min_cells_express_U=0)
            self.filter_genes(by_detection_levels=True)

            if n_top_genes is not None and n_top_genes < self.S.shape[0]:
                self.score_cv_vs_mean(n_top_genes, max_expr_avg=max_expr_avg)
                self.filter_genes(by_cv_vs_mean=True)

            self.score_detection_levels(min_expr_counts=0, min_cells_express=0,
                                        min_expr_counts_U=min_counts_u, min_cells_express_U=min_cells_u)
            self.filter_genes(by_detection_levels=True)

            self._normalize_S(relative_size=self.initial_cell_size, target_size=np.median(self.initial_cell_size))
            self._normalize_U(relative_size=self.initial_Ucell_size, target_size=np.median(self.initial_Ucell_size))

            #self.normalize_by_total()
            print("Number of genes to be used:", self.S.shape[0])

        def impute(self, n_pcs=30, n_neighbors=30, balanced=False, renormalize=False):
            # counterpart to scv.pp.moments()
            self.perform_PCA(n_components=n_pcs)
            k = n_neighbors
            self.knn_imputation(n_pca_dims=n_pcs, k=k, balanced=balanced, b_sight=k * 8, b_maxl=k * 4)
            if renormalize: self.normalize_median()

        def velocity_estimation(self, limit_gamma=False, fit_offset=False):
            self.fit_gammas(limit_gamma=limit_gamma, fit_offset=fit_offset)
            self.filter_genes_good_fit()

            self.predict_U()
            self.calculate_velocity()
            self.calculate_shift()
            self.extrapolate_cell_at_t()
            print("Number of genes to be used:", self.S.shape[0])

        def velocity_graph(self, n_neighbors=100, transform='linear', sampled_fraction=.5, expression_scaling=False,
                           sigma_corr=.05, calculate_randomized=False):
            if not hasattr(self, 'ts'):
                raise ValueError('Compute embedding first.')
            else:
                # counterpart to scv.tl.velocity_graph()
                self.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform=transform, n_neighbors=n_neighbors,
                                              knn_random=True, sampled_fraction=sampled_fraction, calculate_randomized=calculate_randomized)

                # counterpart to scv.tl.velocity_embedding()
                self.calculate_embedding_shift(sigma_corr=sigma_corr, expression_scaling=expression_scaling)

        def velocity_embedding(self, smooth=0.5, steps=(50, 50), n_neighbors=100):
            self.calculate_grid_arrows(smooth=smooth, steps=steps, n_neighbors=n_neighbors)

        def run_all(self, min_counts=3, min_counts_u=3, n_pcs=30, n_neighbors=30, n_neighbors_graph=100,
                    n_top_genes=None, fit_offset=False, limit_gamma=False, transform='linear', expression_scaling=False,):
            from time import time
            start = time()

            self.filter_and_normalize(min_counts=min_counts, min_counts_u=min_counts_u, n_top_genes=n_top_genes)
            print('Preprocessing: ' + str(round(time() - start, 2)))
            timestamp = time()

            self.impute(n_pcs=n_pcs, n_neighbors=n_neighbors)
            print('Imputation: ' + str(round(time() - timestamp, 2)))
            timestamp = time()

            self.velocity_estimation(limit_gamma=limit_gamma, fit_offset=fit_offset)
            print('Velocity Estimation: ' + str(round(time() - timestamp, 2)))
            timestamp = time()

            self.velocity_graph(n_neighbors=n_neighbors_graph, transform=transform, expression_scaling=expression_scaling)
            print('Velocity Graph: ' + str(round(time() - timestamp, 2)))

            print('Total: ' + str(round(time() - start, 2)))

    return VelocytoLoom(adata, basis)
