import numpy as np
from sklearn.decomposition import PCA
from ..data import (
    g2m_markers,
    g1s_markers,
    g2m_inhibitory_markers,
    histone_markers,
)
import warnings


def find_cc_components(adata, thr=3, verbose=True):
    """Find cell cycle related components using cell cycle signature annotations

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
    
    #-- Get scores
    g1s_scores = _compute_scores(adata, g1s_markers)
    g2m_scores = _compute_scores(adata, g2m_markers)
    g2mi_scores = _compute_scores(adata, g2m_inhibitory_markers)
    histone_scores = _compute_scores(adata, histone_markers)
    
    #-- Save scores before finding optimal indices
    cc_scores = {
        'G1-S': g1s_scores,
        'G2-M': g2m_scores,
        'G2-M-': g2mi_scores,
        'Histone': histone_scores
        }
    
    #-- Get idx of max scores
    
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
        
    #-- Check threshold
    g1s_fail = g1s_score < thr
    g2m_fail = g2m_score < thr
    g2mi_fail = g2mi_score < thr
    hist_fail = histone_score < thr

    if g1s_fail or g2m_fail:
        raise Exception("We couldn't find cell-cycle related components. Try decreasing the threshold.")
        
    indices = [g1s_idx, g2m_idx]
    idx_dict = {'G1-S': g1s_idx, 'G2-M': g2m_idx}
    if g2mi_fail:
        print(
            "Warning: G2-M inhibition component missing (score=%f < %f)."
            % (g2mi_score, thr)
        )
    else: 
        indices.append(g2mi_idx)
        idx_dict['G2-M-'] = g2mi_idx
    if hist_fail:
        print(
            "Warning: Histone component missing (score=%f < %f)."
            % (histone_score, thr)
        )
    else: 
        indices.append(histone_idx)
        idx_dict['Histone'] = histone_idx

    # -- Print if verbose
    if verbose:
        print("--- Selected components:",
                "G1-S: %i (score=%f)" % (g1s_idx, g1s_score),
                "G2-M: %i (score=%f)" % (g2m_idx, g2m_score), 
                sep ='\n')
        if not g2mi_fail:
            print("G2-M-: %i (score=%f)" % (g2mi_idx, g2mi_score))
        if not hist_fail:
            print("Histones: %i (score=%f)" % (histone_idx, histone_score))
            
    #-- Get top components from scores
    xdr = adata.obsm["X_dimRed"][:, indices]
    xdr3d = PCA(n_components=3).fit_transform(xdr)
    
    #-- Overwrite
    adata.obsm["X_cc"] = xdr
    adata.obsm["X_pca_scycle"] = xdr3d
    adata.uns["scycle"]["find_cc_components"] = {
        'scores': cc_scores,
        'indices': idx_dict
        }

def _compute_scores(adata, marker_genes):
    marker_genes = adata.var_names.intersection(marker_genes)
    positions = [adata.var_names.get_loc(g) for g in marker_genes]
    return adata.uns["dimRed"].S_[:, positions].mean(axis=1)

def enrich_components(adata, thr=3, verbose=True):
    """DEPRECATED: use find_cc_components
    """
    warnings.warn("enrich_components is deprecated; use find_cc_components", DeprecationWarning)
    return (find_cc_components(adata, thr, verbose))
    
