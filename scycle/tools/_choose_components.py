#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from anndata import AnnData
def choose_components (adata: AnnData, comps: list=[0,1]):
    """ Choose components for 2D projection
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object to be selected.
    comps: list
        A list of two integers representing the two components to use by default
        for the projection plots and to perform the principal circle computation.
        
    Returns
    -------------
    `adata` will be modified to use the chosen components in the default 2D
    projections, and for principal circle calculation.
    """
    if (len(comps) != 2):
        raise Exception('`comps` must be a list with 2 integers.')
    adata.obsm['X_dimRed2d'] = adata.obsm['X_dimRed'][:,comps]