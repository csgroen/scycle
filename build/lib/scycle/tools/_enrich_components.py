import numpy as np
from sklearn.decomposition import PCA
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
    
    #-- Print if verbose
    if verbose:
        print('--- Selected components:',
              'G1/S: ' + str(g1s_idx),
              'G2/M: ' + str(g2m_idx),
              'G2/M-: ' + str(g2mi_idx),
              'Histones: ' + str(histone_idx), sep = '\n')
    
    #-- Update dimRed, dimRed2d
    xdr = adata.obsm['X_dimRed']
    xdr = xdr[:,[g1s_idx, g2m_idx, g2mi_idx, histone_idx]]
    xdr2d = PCA(n_components = 3).fit_transform(xdr)
    
    adata.obsm['X_dimRed'] = xdr
    adata.obsm['X_dimRed3d'] = xdr2d
    
    # -- Return
    adata.uns["scycle"]["enrich_components"] = {
        "G1/S": g1s_idx,
        "G2/M+": g2m_idx,
        "G2/M-": g2mi_idx,
        "Histone": histone_idx,
    }


def _compute_scores(adata, marker_genes):
    marker_genes = [g for g in marker_genes if g in adata.var_names]
    sel = [(g in marker_genes) for g in adata.var_names]
    return np.mean(adata.uns["dimRed"].S_[:, sel], axis=1)
