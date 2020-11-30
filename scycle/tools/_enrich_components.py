import numpy as np
from ..data import (
    g2m_markers,
    g1s_markers,
    g2m_inhibitory_markers,
    histone_markers,
)


def enrich_components(adata, verbose=True):
    """Pathway enrichment for components from the dimensionality reduction

    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        tl.dimensionality_reduction.
    verbose: bool
        If True, messages about function progress will be printed.

    Returns
    ------------------
    `adata` will be updated with the results of the component enrichment and
    automatic selection attempt using the function parameters.
    """

    # -- Get components
    assert "dimRed" in adata.uns, "Dimensionality reduction must be performed first."

    g1s_scores = _compute_scores(adata, g1s_markers)
    g2m_scores = _compute_scores(adata, g2m_markers)
    g2m_inhibitory_scores = _compute_scores(adata, g2m_inhibitory_markers)
    histone_scores = _compute_scores(adata, histone_markers)

    # -- Return
    adata.uns["scycle"]["enrich_components"] = {
        "G1/S": np.argmax(g1s_scores),
        "G2/M+": np.argmax(g2m_scores),
        "G2/M-": np.argmax(g2m_inhibitory_scores),
        "Histone": np.argmax(histone_scores),
    }


def _compute_scores(adata, marker_genes):
    marker_genes = [g for g in marker_genes if g in adata.var_names]
    sel = [(g in marker_genes) for g in adata.var_names]
    return np.mean(adata.uns["dimRed"].S_[:, sel], axis=0)
