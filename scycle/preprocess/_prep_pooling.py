#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
from sklearn.neighbors import NearestNeighbors
from anndata import AnnData

def prep_pooling(adata: AnnData, dim_red_method_pooling: str='pca', 
                 n_neighbors: int=10, min_counts: int=1, target_sum: int=10000, 
                  filter_var_genes: bool=True, n_top_genes: int=2000,
                 n_bins: int=20, embed_n_comps: int=30,  verbose: bool=True):
    """Pre-processes AnnData without pooling
    
    Parameters
    ----------
    adata: AnnData
        The AnnData object to be pre-processed. This should already have been
        processed to remove "bad cells" (high mitochondrial percentage,
        aberrant total counts).
    dim_red_method_pooling: str
        Method to use for dimensionality reduction to do the pooling procedure.
        Default: 'pca'. TO-DO: support 'ica' and other?
    n_neighbors: int
        Number of nearest neighbors to use for pooling.
    min_counts: int
        Minimum number of counts required for a gene to pass filtering. Passed
        to `pp.prep_simple`.
    target_sum: int
        Target sum of counts for library normalization. Passed to
        `pp.prep_simple`.
    filter_var_genes: bool
        If True, only `n_top_genes` highly variable genes are kept. Passed to
        `pp.prep_simple`.
    n_top_genes: int
        Number of genes to keep after highly variable filter. Used if 
        `filter_var_genes` is True. Passed to sc.pp.highly_variable_genes.
    n_bins: int
        Number of bins for binning the mean gene expression. Normalization is 
        done with respect to each bin. If just a single gene falls into a bin, 
        the normalized dispersion is artificially set to 1. 
        Passed to`pp.prep_simple`.
    embed_n_comps: int
        Number of components to use for the embedding to do the pooling.
    verbose: bool
        If True, messages about function progress will be printed.
    
    Returns
    ----------
    None
    """
    
    if 'scycle' in adata.uns:
        raise Exception('Data has already been pre-processed')
    
    if verbose: print("Preparing embedding...")
    adata_simple = adata.copy()
    prep_simple(adata_simple, verbose = False)

    if verbose: print("Embedding for pooling...")
    X_embed = _embed_for_pooling(adata_simple, dim_red_method_pooling, n_comps = embed_n_comps)
        
    if verbose: print("Pooling", str(X_embed.shape[0]), "samples...")
    _smooth_adata_by_pooling(adata, X_embed, n_neighbours = n_neighbors)
    prep_simple(adata, min_counts = min_counts, target_sum = target_sum, 
                n_top_genes = n_top_genes, n_bins = n_bins, 
                filter_var_genes = filter_var_genes, verbose = verbose)
    adata.uns['scycle'] = {'preprocess': {'method': 'pooling', 'n_neighbors': n_neighbors,
                                          'min_counts': min_counts, 'target_sum': target_sum,
                                          'n_top_genes': n_top_genes, 'n_bins': n_bins,
                                          'embed_n_comps': embed_n_comps}}
    
def _embed_for_pooling(obj, dim_red, n_comps):
    if(dim_red == 'pca'):
        sc.tl.pca(obj, n_comps = n_comps)
        X_embed = obj.obsm['X_pca']
        return(X_embed)
    
def _smooth_adata_by_pooling(adata, X_embed, n_neighbours=10, copy = False):
    # adata_pooled = adata.copy if copy else adata
    nbrs = NearestNeighbors(n_neighbors=n_neighbours).fit(X_embed)
    distances, indices = nbrs.kneighbors(X_embed)    
    adata.X = _smooth_matrix_by_pooling(_get_nd_array(adata.X),indices)
    if 'matrix' in adata.layers:
        adata.layers['matrix'] = _smooth_matrix_by_pooling(_get_nd_array(adata.layers['matrix']),indices)
    if 'spliced' in adata.layers:
        adata.layers['spliced'] = _smooth_matrix_by_pooling(_get_nd_array(adata.layers['spliced']),indices)
    if 'unspliced' in adata.layers:
        adata.layers['unspliced'] = _smooth_matrix_by_pooling(_get_nd_array(adata.layers['unspliced']),indices)
    
    adata.uns['scycle'] = True
    

def _smooth_matrix_by_pooling(matrix,indices):
    matrix_pooled = matrix.copy()
    for i in range(len(indices)):
        matrix_pooled[i,:] = np.mean(matrix[indices[i],:],axis=0)
    return matrix_pooled

def _get_nd_array(arr):
    x = None
    if str(type(arr)):
        x = arr
    else:
        x = arr.toarray()
    return x
