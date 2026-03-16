from typing import Optional, Union

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from anndata import AnnData

from scvelo.plotting.scatter import scatter
from scvelo.plotting.utils import savefig_or_show


def kinetic_clusters(
    adata: AnnData,
    mode: str = "pca",
    basis: str = "umap",
    show: bool = True,
    save: Optional[str] = None,
    **kwargs,
) -> Union[Axes, None]:
    """
    Visualize kinetic clusters.

    Modes:
    - 'pca': Scatter plot of Log-Transformed kinetic parameters.
    - 'embedding': Project kinetic scores onto the cell embedding (e.g. UMAP).

    Parameters
    ----------
    adata
        Annotated data matrix.
    mode
        'pca' for parameter space, 'embedding' for cell space.
    basis
        Basis for embedding (e.g. 'umap', 'tsne'). Used if mode='embedding'.
    show
        Whether to display the plot.
    save
        Path to save the plot, or None.
    **kwargs
        Arguments passed to matplotlib or scvelo.pl.scatter.
    """
    if "kinetic_cluster" not in adata.var.keys():
        raise ValueError("Please run scv.tl.kinetic_clusters(adata) first.")

    # MODE 1: THE MATH PLOT (Log-Log Space)
    if mode == "pca":
        mask = adata.var["kinetic_cluster"] != "nan"
        df = adata.var.loc[mask].copy()

        clusters = df["kinetic_cluster"].unique()
        try:
            clusters = sorted(clusters, key=int)
        except ValueError:
            clusters = sorted(clusters)

        fig, ax = plt.subplots(figsize=(6, 5))

        for cluster_id in clusters:
            subset = df[df["kinetic_cluster"] == cluster_id]

            x_vals = np.log1p(subset["fit_gamma"])
            y_vals = np.log1p(subset["fit_alpha"])

            ax.scatter(
                x_vals,
                y_vals,
                label=f"Cluster {cluster_id}",
                s=kwargs.get("s", 20),
                alpha=kwargs.get("alpha", 0.7),
                edgecolors="none",
            )

        min_val = min(np.log1p(df["fit_gamma"].min()), np.log1p(df["fit_alpha"].min()))
        max_val = max(np.log1p(df["fit_gamma"].max()), np.log1p(df["fit_alpha"].max()))
        ax.plot(
            [min_val, max_val],
            [min_val, max_val],
            "--",
            color="grey",
            alpha=0.3,
            label="Steady State",
        )

        ax.set_xlabel("Log Degradation Rate (Gamma)")
        ax.set_ylabel("Log Transcription Rate (Alpha)")
        ax.set_title("Gene Kinetic Regimes (Log Scale)")
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.3)

        savefig_or_show(save=save, show=show)
        return ax if not show else None

    # MODE 2: THE BIOLOGICAL PLOT (UMAP)
    elif mode == "embedding":
        if not any("Kinetic_Cluster_" in col for col in adata.obs.columns):
            from scvelo.tools.kinetic_clusters import score_kinetic_clusters

            score_kinetic_clusters(adata)

        clusters = adata.var["kinetic_cluster"].unique()
        clusters = sorted([c for c in clusters if str(c) != "nan"])
        keys = [f"Kinetic_Cluster_{c}" for c in clusters]

        ax = scatter(
            adata,
            basis=basis,
            color=keys,
            cmap="viridis",
            title=[f"Regime {c}" for c in clusters],
            frameon=False,
            show=show,
            save=save,
            **kwargs,
        )
        return ax

    else:
        raise ValueError(f"Unknown mode '{mode}'. Use 'pca' or 'embedding'.")
