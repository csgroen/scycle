import numpy as np
from sklearn.decomposition import PCA
from ..data import (
    g2m_markers,
    g1s_markers,
    g2m_inhibitory_markers,
    histone_markers,
)


def enrich_components(adata, thr=3, verbose=True):
    """Pathway enrichment for components from the dimensionality reduction

    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        tl.dimensionality_reduction.
    thr: float
      Score threshold for IC recognition
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
    g2mi_scores = _compute_scores(adata, g2m_inhibitory_markers)
    histone_scores = _compute_scores(adata, histone_markers)

    g2m_idx = np.argmax(g2m_scores)
    g1s_scores[g2m_idx] = 0
    g2mi_scores[g2m_idx] = 0
    histone_scores[g2m_idx] = 0

    g1s_idx = np.argmax(g1s_scores)
    g2mi_scores[g1s_idx] = 0
    histone_scores[g1s_idx] = 0

    g2mi_idx = np.argmax(g2mi_scores)
    histone_scores[g2mi_idx] = 0

    histone_idx = np.argmax(histone_scores)

    g1s_score = g1s_scores[g1s_idx]
    g2m_score = g2m_scores[g2m_idx]
    g2mi_score = g2mi_scores[g2mi_idx]
    histone_score = histone_scores[histone_idx]

    if verbose:
        print("G1/S  score: %f" % g1s_score)
        print("G2/M+ score: %f" % g2m_score)
        print("G2/M- score: %f" % g2mi_score)
        print("HIST  score: %f" % histone_score)

    if g1s_score < thr or g2m_score < thr or histone_score < thr:
        print("ERROR: low cell-cycle related signatures. Try decreasing the threshold.")
        return

    g2mi_found = g2mi_score >= thr
    if not g2mi_found:
        print(
            "Warning: only three cell-cycle signatures found (out of 4), G2/M inhibition missing (score=%f < %f)."
            % (g2mi_score, thr)
        )

    # -- Print if verbose
    if verbose:
        if g2mi_found:
            print(
                "--- Selected components:",
                "G1/S: %i (score=%f)" % (g1s_idx, g1s_score),
                "G2/M: %i (score=%f)" % (g2m_idx, g2m_score),
                "G2/M-: %i (score=%f)" % (g2mi_idx, g2mi_score),
                "Histones: %i (score=%f)" % (histone_idx, histone_score),
                sep="\n",
            )
        else:
            print(
                "--- Selected components:",
                "G1/S: %i (score=%f)" % (g1s_idx, g1s_score),
                "G2/M: %i (score=%f)" % (g2m_idx, g2m_score),
                "Histones: %i (score=%f)" % (histone_idx, histone_score),
                sep="\n",
            )

    # -- Update dimRed, pc3
    if g2mi_found:
        indices = [g1s_idx, g2m_idx, g2mi_idx, histone_idx]
    else:
        indices = [g1s_idx, g2m_idx, histone_idx]

    xdr = adata.obsm["X_dimRed"][:, indices]
    xdr3d = PCA(n_components=3).fit_transform(xdr)

    adata.obsm["X_cc"] = xdr
    adata.obsm["X_pca_scycle"] = xdr3d

    adata.uns["scycle"]["enrich_components"] = {}

    adata.uns["scycle"]["enrich_components"]["G1/S"] = g1s_idx
    adata.uns["scycle"]["enrich_components"]["G2/M+"] = g2m_idx
    if g2mi_found:
        adata.uns["scycle"]["enrich_components"]["G2/M-"] = g2mi_idx
    adata.uns["scycle"]["enrich_components"]["Histone"] = histone_idx


def _compute_scores(adata, marker_genes):
    marker_genes = adata.var_names.intersection(marker_genes)
    positions = [adata.var_names.get_loc(g) for g in marker_genes]
    return adata.uns["dimRed"].S_[:, positions].mean(axis=1)
