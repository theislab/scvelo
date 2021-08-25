import numpy as np
import pandas as pd

from scvelo import logging as logg

# fmt: off
s_genes_list = \
    ['Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
     'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Mlf1ip', 'Hells', 'Rfc2',
     'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
     'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
     'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8']

g2m_genes_list = \
    ['Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
     'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
     'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
     'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
     'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
     'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
     'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa']
# fmt: on


def get_phase_marker_genes(adata):
    """\
    Return S and G2M phase marker genes.

    Parameters
    ----------
    adata
        The annotated data matrix.
    phase

    Returns
    -------
    List of S and/or G2M phase marker genes.
    """

    name, gene_names = adata.var_names[0], adata.var_names
    up, low = name.isupper(), name.islower()
    s_genes_list_ = [
        x.upper() if up else x.lower() if low else x.capitalize() for x in s_genes_list
    ]
    g2m_genes_list_ = [
        x.upper() if up else x.lower() if low else x.capitalize()
        for x in g2m_genes_list
    ]
    s_genes = np.array([x for x in s_genes_list_ if x in gene_names])
    g2m_genes = np.array([x for x in g2m_genes_list_ if x in gene_names])
    return s_genes, g2m_genes


def score_genes_cell_cycle(adata, s_genes=None, g2m_genes=None, copy=False, **kwargs):
    """\
    Score cell cycle genes.

    Calculates scores and assigns a cell cycle phase (G1, S, G2M) using the list of cell
    cycle genes defined in Tirosh et al, 2015 (https://doi.org/10.1126/science.aad0501).

    Parameters
    ----------
    adata
        The annotated data matrix.
    s_genes
        List of genes associated with S phase.
    g2m_genes
        List of genes associated with G2M phase.
    copy
        Copy `adata` or modify it inplace.
    **kwargs
        Are passed to :func:`~scanpy.tl.score_genes`. `ctrl_size` is not
        possible, as it's set as `min(len(s_genes), len(g2m_genes))`.

    Returns
    -------
    S_score: `adata.obs`, dtype `object`
        The score for S phase for each cell.
    G2M_score: `adata.obs`, dtype `object`
        The score for G2M phase for each cell.
    phase: `adata.obs`, dtype `object`
        The cell cycle phase (`S`, `G2M` or `G1`) for each cell.
    """

    logg.info("calculating cell cycle phase")
    from scanpy.tools._score_genes import score_genes

    adata = adata.copy() if copy else adata

    s_genes_, g2m_genes_ = get_phase_marker_genes(adata)
    if s_genes is None:
        s_genes = s_genes_
    if g2m_genes is None:
        g2m_genes = g2m_genes_

    ctrl_size = min(len(s_genes), len(g2m_genes))

    kwargs.update({"ctrl_size": ctrl_size})
    score_genes(adata, gene_list=s_genes, score_name="S_score", **kwargs)
    score_genes(adata, gene_list=g2m_genes, score_name="G2M_score", **kwargs)
    scores = adata.obs[["S_score", "G2M_score"]]

    phase = pd.Series("S", index=scores.index)  # default phase is S
    phase[scores.G2M_score > scores.S_score] = "G2M"  # G2M, if G2M is higher than S
    phase[np.all(scores < 0, axis=1)] = "G1"  # G1, if all scores are negative

    adata.obs["phase"] = phase
    logg.hint("    'S_score' and 'G2M_score', scores of cell cycle phases (adata.obs)")
    return adata if copy else None
