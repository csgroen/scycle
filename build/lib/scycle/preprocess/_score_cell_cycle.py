#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional

import scanpy as sc
import numpy as np
from anndata import AnnData
from ..annot import cellcycle_signatures


def score_cell_cycle(adata: AnnData, sigs: Optional[list]=None, 
                     scale: bool=True, verbose: bool=True):
    if sigs != None:
        if type(sigs) != dict:
            raise TypeError('sigs must be a dictionary')
        if (('S-phase' not in sigs) | ('G2-M' not in sigs)):
            raise Exception('sigs must have at least two signatures: one named "S-phase" and the other named "G2-M)"')
    
    _score_cell_cycle(adata, sigs = sigs, scale = scale, verbose = verbose)
    
def _score_cell_cycle(adata, sigs = None, scale = True, verbose = True):
    #-- Get signatures
    if sigs == None:
        cc_sigs = cellcycle_signatures()
    else:
        cc_sigs = sigs
        
    if verbose: print('-- Scoring G1 genes...')
    sc.tl.score_genes(adata, cc_sigs['G1'], score_name = 'G1')
    if verbose: print("-- Scoring S-phase...")
    sc.tl.score_genes(adata,  cc_sigs['S-phase'], score_name = 'S-phase')

    if verbose: print("-- Scoring G2-M...")
    sc.tl.score_genes(adata, cc_sigs['G2-M'], score_name = 'G2-M')
    
    if verbose: print("-- Scoring Histones...")
    adata.obs['Histones'] = _calc_histone_score(adata, verbose = verbose)
    
    if scale:
        if verbose: print("-- Scalling signatures...")
        sigs = ['G2-M', 'S-phase', 'Histones']
        for sig in sigs: 
            adata.obs[sig] = adata.obs[sig]/np.max(adata.obs[sig])
            
def _calc_histone_score(adata2k, verbose = True):
    #-- Find histone names
    hist_inds = adata2k.var_names.str.match('^H[1-4]')
    histone_names = adata2k.var_names[hist_inds]
    if verbose: print('Found histone genes:',*histone_names)
    #-- Calculate scores
    matrix = adata2k.to_df().to_numpy()
    matrix_sel = matrix[:,hist_inds]
    scores = np.mean(matrix_sel,axis=1)
    return scores