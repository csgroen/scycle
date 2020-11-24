#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from anndata import AnnData
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

def select_components (adata: AnnData, comps = None):
    """ Choose components for 2D projection
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object to be selected. Must be previously evaluated by
        tl.dimensionality_reduction.
    comps: list
        A list of integers representing components to be projected by PCA into
        2d for finding the cell cycle trajectory. If None, will use the suggestions
        from the automatic selection made by tl.enrich_components.
        
    Returns
    -------------
    `adata` will be modified to recompute to do a PCA for projection using only
    the chosen components.
    """
    if comps is None:
        comps = adata.uns['scycle']['enrich_components']['suggested_comps']
    
    dr = adata.obsm['X_dimRed'][:,comps]
    pca = PCA(n_components = 2).fit(dr)

    adata.obsm['X_dimRed2d'] = pca.transform(dr)
    
    if adata.uns['scycle']['dimRed']['method'] == 'ica':
        ica = adata.uns['dimRed']
        adata.uns['dimRed'] = _update_ica(ica, comps)
        adata.obsm['X_dimRed'] = dr
        adata.uns['scycle']['dimRed']['select_comps'] = comps
    
def _update_ica(ica, comps):
    ica.components_ = ica.components_[:,comps]    
    ica.whitening_ = ica.whitening_[comps,:]
    ica.n_components = len(comps)
    ica.mixing_ = ica.mixing_[:,comps]
    return ica