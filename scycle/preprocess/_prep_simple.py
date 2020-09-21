#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
from anndata import AnnData

def prep_simple(adata: AnnData, min_counts: int= 1, target_sum: int=10000, 
                filter_var_genes: bool=True, n_top_genes: int=2000, 
                n_bins: int=20, for_pooling: bool=False, verbose: bool=True):
    """Pre-processes AnnData without pooling
    
    Parameters
    ----------
    adata: AnnData
        The AnnData object to be pre-processed. This should already have been
        processed to remove "bad cells" (high mitochondrial percentage,
        aberrant total counts).
    min_counts: int
        Minimum number of counts required for a gene to pass filtering. Passed
        to sc.pp.filter_genes
    target_sum: int
        Target sum of counts for library normalization. Passed to
        sc.pp.normalize_total.
    filter_var_genes: bool
        If True, only `n_top_genes` highly variable genes are kept.
    n_top_genes: int
        Number of genes to keep after highly variable filter. Used if 
        `filter_var_genes` is True. Passed to sc.pp.highly_variable_genes.
    n_bins: int
        Number of bins for binning the mean gene expression. Normalization is 
        done with respect to each bin. If just a single gene falls into a bin, 
        the normalized dispersion is artificially set to 1. 
        Passed to sc.pp.highly_variable_genes.
    for_pooling: bool
        Set to True if the function is called by the `prep_pooling` function.
        Changes the return object parameters.
    verbose: bool
        If True, messages about function progress will be printed.
    
    Returns
    ----------
    None
        
    """
    sc.pp.filter_genes(adata, min_counts=min_counts)
    sc.pp.normalize_total(adata, target_sum = target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = n_top_genes, n_bins = n_bins)
    ind_genes = np.where(adata.var['highly_variable'])[0]
    adata = adata[:,ind_genes] if filter_var_genes else adata
    if not for_pooling: 
        adata.uns['scycle'] = {'preprocess': {'method': 'simple', 'min_counts': min_counts, 
                                              'target_sum': target_sum, 'filter_var_genes': filter_var_genes, 
                                              'n_top_genes': n_top_genes, 'n_bins': n_bins}}