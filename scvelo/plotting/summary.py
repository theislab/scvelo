import numpy as np
import pandas as pd

from scvelo.tools.dynamical_model import latent_time
from scvelo.tools.rank_velocity_genes import rank_velocity_genes
from scvelo.tools.score_genes_cell_cycle import score_genes_cell_cycle
from scvelo.tools.velocity_pseudotime import velocity_pseudotime
from .gridspec import GridSpec
from .utils import make_unique_list


def summary(adata, basis="umap", color="clusters", n_top_genes=12, var_names=None):
    if "fit_alpha" in adata.var.keys():
        tkey = "latent_time"
        if tkey not in adata.obs.keys():
            latent_time(adata)
        top_genes = adata.var.sort_values("fit_likelihood", ascending=False)
        top_genes = top_genes.index[:n_top_genes]
    else:
        tkey = "velocity_pseudotime"
        velocity_pseudotime(adata)
        if "rank_velocity_genes" not in adata.uns.keys():
            rank_velocity_genes(adata, groupby=color)
        top_genes = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
        top_genes = top_genes.iloc[0][:n_top_genes]
    nrows = 1 + int(np.ceil(len(top_genes) / 4))

    if "S_score" not in adata.obs.keys() or "G2M_score" not in adata.obs.keys():
        score_genes_cell_cycle(adata)
    adata.obs["S"] = np.clip(adata.obs["S_score"], 0.1, None)
    adata.obs["G2M"] = np.clip(adata.obs["G2M_score"], 0.1, None)

    if var_names is not None:
        var_names = make_unique_list(var_names)
        var_names = [name for name in var_names if name in adata.var_names]
        nrows += int(np.ceil(len(var_names) / 4))

    kwargs = dict(
        c=color, legend_loc_lines="none", add_outline="fit_diff_kinetics", frameon=False
    )
    with GridSpec(ncols=4, nrows=nrows) as pl:
        pl.scatter(adata, basis=basis, c=color)
        pl.velocity_embedding_stream(adata, c=tkey)
        pl.scatter(adata, color_gradients=["root_cells", "end_points"], palette="jet")
        pl.scatter(adata, color_gradients=["S", "G2M"])

        if top_genes is not None:
            for gene in top_genes:
                text = r"$L=$ " + f"{np.round(adata[:, gene].var.fit_likelihood[0], 2)}"
                pl.scatter(adata, basis=gene, add_text=text, **kwargs)

        if var_names is not None:
            for gene in var_names:
                text = r"$L=$ " + f"{np.round(adata[:, gene].var.fit_likelihood[0], 2)}"
                pl.scatter(adata, basis=gene, add_text=text, **kwargs)
