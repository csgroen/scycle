#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional

import scanpy as sc
import numpy as np
from anndata import AnnData
from ..data import (
    g2m_markers,
    g1s_markers,
    histone_markers,
    g2m_markers_short,
    g1s_markers_short,
    G2M_genes_Freeman,
    G1S_genes_Freeman
)

def score_cell_cycle(adata: AnnData, verbose: bool=True):
    if verbose:
        print('Scoring cell cycle...')
    _score_cell_cycle(adata, g1s_markers, "G1S_Tirosh")
    _score_cell_cycle(adata, g2m_markers, "G2M_Tirosh")
    _score_cell_cycle(adata, G1S_genes_Freeman, "G1S_Freeman")
    _score_cell_cycle(adata, G2M_genes_Freeman, "G2M_Freeman")
    _score_cell_cycle(adata, g1s_markers_short, "G1S_short")
    _score_cell_cycle(adata, g2m_markers_short, "G2M_short")
    _score_cell_cycle(adata, histone_markers, "Histones")

    adata.obs['G1-S'] = adata.obs['G1S_Tirosh']
    adata.obs['G2-M'] = adata.obs['G2M_Tirosh']

def _score_cell_cycle(adata, markers, score_name):
    sc.tl.score_genes(adata, markers, score_name=score_name)
    adata.obs[score_name] = zscore(adata.obs[score_name])

# def _score_cell_cycle(adata, sigs = None, scale = True, verbose = True):
#     #-- Get signatures
#     if sigs == None:
#         cc_sigs = cellcycle_signatures()
#     else:
#         cc_sigs = sigs
#
#     if verbose: print('-- Scoring G1 genes...')
#     sc.tl.score_genes(adata, cc_sigs['G1'], score_name = 'G1')
#     if verbose: print("-- Scoring S-phase...")
#     sc.tl.score_genes(adata,  cc_sigs['S-phase'], score_name = 'S-phase')
#
#     if verbose: print("-- Scoring G2-M...")
#     sc.tl.score_genes(adata, cc_sigs['G2-M'], score_name = 'G2-M')
#
#     if verbose: print("-- Scoring Histones...")
#     adata.obs['Histones'] = _calc_histone_score(adata, verbose = verbose)
#
#     if scale:
#         if verbose: print("-- Scalling signatures...")
#         sigs = ['G2-M', 'S-phase', 'Histones']
#         for sig in sigs:
#             adata.obs[sig] = adata.obs[sig]/np.max(adata.obs[sig])
#
# def _calc_histone_score(adata2k, verbose = True):
#     #-- Find histone names
#     hist_inds = adata2k.var_names.str.match('^H[1-4]')
#     histone_names = adata2k.var_names[hist_inds]
#     if verbose: print('Found histone genes:',*histone_names)
#     #-- Calculate scores
#     matrix = adata2k.to_df().to_numpy()
#     matrix_sel = matrix[:,hist_inds]
#     scores = np.mean(matrix_sel,axis=1)
#     return scores
