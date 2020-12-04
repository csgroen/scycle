#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:50:39 2020

@author: clarice
"""
import numpy as np
from ._prep_simple import quality_control
from ._prep_pooling import prep_pooling
from ..tools import  dimensionality_reduction, enrich_components, principal_circle

def normalize_by_partition(adata, 
                           rerun_pc = True, 
                           n_ref_parts = 10, 
                           verbose = True):
    
    """ Normalize samples by median partition library size
     
     This procedure improves the discovery of the cell division moment.
    
    Parameters
    ------------
    adata_src: AnnData
        The analysis object to be evaluated. Should be previously evaluated by
    rerun_pc: bool
        If True, `tl.principal circle` is re-run with n_ref_parts as n_nodes.
        Using fewer partitions for re-normalization than the default number used
    n_ref_parts: int
        How many partitions to be used for library size normalization.
    verbose: bool
        If True, the function will print messages.

    """ 
    adata_src = adata.copy()

    if 'principal_circle' not in adata.uns['scycle'].keys() | rerun_pc:
        if verbose: print("-- Running `tl.principal_circle` with n_ref_parts...")
        principal_circle(adata, n_nodes = n_ref_parts, verbose = False)
               

    #-- Run QC
    params = adata.uns['scycle']
    pp_params = params['preprocess']
    dr_params = params['dimRed']

    quality_control(adata_src,
                    min_counts = pp_params['min_counts'], 
                    max_counts = pp_params['max_counts'], 
                    max_mt_ratio = pp_params['max_mt_ratio'],
                    verbose = False)
    old_totals = adata_src.obs['total_counts']
    
    #--- Apply filter
    if verbose: print('Normalizing by partition...')

    #---- Get partitions and re-noralize
    prt = adata.obs['partition']
    gexp = adata_src.X
    
    npart = np.max(prt)+1
    new_gexp = np.empty(gexp.shape)
    for p in range(npart):
        sidx = prt == p # sample index
        totals = np.sum(gexp[sidx,:], axis = 1) # total counts per sample in group
        median = np.median(totals) # median counts for samples in group
        new_gexp[sidx,:] = gexp[sidx,:] / totals[:,None] * median
        
    #---- Re-run procedure
    adata_src.X = new_gexp
    
    if verbose: print('Re-running pooling...')
    
    prep_pooling(adata_src, 
                 filter_cells = False,
                 embed_n_comps=pp_params['embed_n_comps'],
                 min_counts = pp_params['min_counts'], 
                 max_counts = pp_params['max_counts'], 
                 max_mt_ratio = pp_params['max_mt_ratio'],
                 n_neighbors = pp_params['n_neighbors'],
                 normalize_counts = False,
                 filter_var_genes = pp_params['filter_var_genes'],
                 log_transform= pp_params['log_transform'],
                 n_top_genes = pp_params['n_top_genes'],
                 verbose = False)
    
    if verbose: print('Re-running dimensionality reduction..')
    dimensionality_reduction(adata_src,
                             method = dr_params['method'],
                             n_comps = dr_params['n_comps'],
                             seed = dr_params['seed'],
                             verbose = False)
    if 'enrich_components' in params.keys():    
        if verbose: print('Re-running component enrichment...')
        enrich_components(adata_src, verbose = False)
    
    if verbose: print('Finding the principal circle...')
    principal_circle(adata_src, verbose = False)
    
    adata_src.obs['total_counts'] = np.sum(adata_src.X, axis = 1)
    adata_src.obs['total_counts_raw'] = old_totals