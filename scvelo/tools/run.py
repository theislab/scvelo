from ..preprocessing import filter_and_normalize, moments
from . import velocity, velocity_graph, velocity_embedding


def run_all(data, basis=None, mode='deterministic', min_counts=3, n_pcs=30, n_neighbors=30, copy=False):
    from time import time
    start = time()

    adata = data.copy() if copy else data
    filter_and_normalize(adata, min_counts=min_counts)
    print("Number of genes to be used:", adata.X.shape[1])
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    velocity(adata, mode=mode)
    print("Number of genes to be used:", adata.var.velocity_genes.sum())

    velocity_graph(adata)
    if basis is not None: velocity_embedding(adata, basis=basis)

    print('Total time (seconds): ' + str(round(time() - start, 2)))

    return adata if copy else None


def convert_to_loom(adata, basis=None):
    from scipy.sparse import issparse
    import numpy as np
    import velocyto

    class VelocytoLoom(velocyto.VelocytoLoom):
        def __init__(self, adata, basis=None):
            self.S = adata.layers['spliced'].T
            self.U = adata.layers['unspliced'].T
            self.A = adata.layers['ambiguous'].T

            if issparse(self.S):
                self.S = self.S.A
                self.U = self.U.A
                self.A = self.A.A

            self.initial_cell_size = self.S.sum(0)
            self.initial_Ucell_size = self.U.sum(0)

            self.ca = dict()
            self.ra = dict()
            for key in adata.obs.keys():
                self.ca[key] = np.array(adata.obs[key].values)
            for key in adata.var.keys():
                self.ra[key] = np.array(adata.var[key].values)
            if basis is not None:
                self.ts = adata.obsm['X_' + basis]

            if 'clusters' in self.ca:
                self.set_clusters(self.ca['clusters'])
            elif 'louvain' in self.ca:
                self.set_clusters(self.ca['louvain'])

        def filter_and_normalize(self, min_counts=3, min_counts_u=3, n_top_genes=None):
            # counterpart to scv.pp.filter_and_normalize()
            self.score_detection_levels(min_expr_counts=min_counts, min_expr_counts_U=min_counts_u,
                                        min_cells_express=0, min_cells_express_U=0)
            self.filter_genes(by_detection_levels=True)

            self.normalize('both', size=True, log=False)
            if n_top_genes is not None:
                self.score_cv_vs_mean(n_top_genes)
                self.filter_genes(by_cv_vs_mean=True)

            self.normalize_by_total()
            print("Number of genes to be used:", self.S.shape[0])

        def impute(self, n_pcs=30, n_neighbors=30):
            # counterpart to scv.pp.moments()
            self.perform_PCA(n_components=n_pcs)
            self.knn_imputation(k=n_neighbors, n_pca_dims=n_pcs)
            self.normalize_median()

        def velocity(self, fit_offset=False):
            self.fit_gammas(fit_offset=fit_offset)
            self.filter_genes_good_fit()

            self.predict_U()
            self.calculate_velocity()
            self.calculate_shift()
            self.extrapolate_cell_at_t()
            print("Number of genes to be used:", self.S.shape[0])

        def velocity_graph(self, n_neighbors_graph=100):
            if not hasattr(self, 'ts'):
                raise ValueError('Compute embedding first.')
            else:
                # counterpart to scv.tl.velocity_graph()
                self.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="linear",
                                              n_neighbors=n_neighbors_graph, knn_random=True, sampled_fraction=1)

                # counterpart to scv.tl.velocity_embedding()
                self.calculate_embedding_shift(sigma_corr=0.1, expression_scaling=True)

        def run_all(self, min_counts=3, min_counts_u=3, n_pcs=30, n_neighbors=30, n_neighbors_graph=100,
                    n_top_genes=None, fit_offset=False):
            from time import time
            start = time()

            self.filter_and_normalize(min_counts=min_counts, min_counts_u=min_counts_u, n_top_genes=n_top_genes)
            print('Preprocessing: ' + str(round(time() - start, 2)))
            timestamp = time()

            self.impute(n_pcs=n_pcs, n_neighbors=n_neighbors)
            print('Imputation: ' + str(round(time() - timestamp, 2)))
            timestamp = time()

            self.velocity(fit_offset=fit_offset)
            print('Velocity Estimation: ' + str(round(time() - timestamp, 2)))
            timestamp = time()

            self.velocity_graph(n_neighbors_graph=n_neighbors_graph)
            print('Velocity Graph: ' + str(round(time() - timestamp, 2)))

            print('Total: ' + str(round(time() - start, 2)))

    return VelocytoLoom(adata, basis)
