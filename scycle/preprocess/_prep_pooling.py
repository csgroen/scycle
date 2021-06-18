#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
from sklearn.neighbors import NearestNeighbors
from anndata import AnnData
from ._prep_simple import prep_simple
import gc


def prep_pooling(
    adata: AnnData,
    dim_red_method_pooling: str = "pca",
    n_neighbors: int = 5,
    embed_n_comps: int = 20,
    normalize_counts: bool = True,
    filter_var_genes: bool = True,
    n_top_genes: int = 10000,
    log_transform: bool = True,
    division_factor: float = 1,
    score_cell_cycle: bool = True,
    verbose: bool = True,
):
    """Pre-processes AnnData with cell pooling

    Parameters
    ----------
    adata: AnnData
        The AnnData object to be pre-processed. This should already have been
        processed to remove "bad cells" (high mitochondrial percentage,
        aberrant total counts). See pp.filter_cells
    dim_red_method_pooling: str
        Method to use for dimensionality reduction to do the pooling procedure.
        Default: 'pca'. TO-DO: support 'ica' and other?
    n_neighbors: int
        Number of nearest neighbors to use for pooling.
    embed_n_comps: int
        Number of components to use for the embedding to do the pooling.
    normalize_counts: bool
        Set it to False if library does not need normalization
    filter_var_genes: bool
        If True, only `n_top_genes` highly variable genes are kept.
    n_top_genes: int
        Number of genes to keep after highly variable filter. Used if
        `filter_var_genes` is True. Passed to sc.pp.highly_variable_genes.
    log_transform: bool
        Set it to false if you do not want values to be log-transformed.
    division_factor: int
        Scaling factor, divides the counts matrix by this value.
    score_cell_cycle: bool
        Should cell cycle scores be added?
    verbose: bool
        If True, messages about function progress will be printed.

    Returns
    ----------
    None
    """
    if "scycle" in adata.uns:
        if "preprocess" in adata.uns["scycle"].keys():
            raise Exception("Data has already been pre-processed")

    if verbose:
        print("Preparing embedding...")

    assert division_factor != 0, "Null division factor. Terminating..."
    assert "float" in str(adata.X.dtype), "adata.X dtype must be float."

    # -- sparse -> array
    if 'ndarray' not in str(type(adata.X)):
        adata.X = adata.X.toarray()
    # np.divide(adata.X, division_factor, out=adata.X)

    if len(adata.layers.keys()) == 0:  # keep "raw" data
        adata.layers["matrix"] = adata.X

    adata_simple = adata.copy()
    prep_simple(
        adata_simple,
        normalize_counts,
        filter_var_genes,
        n_top_genes,
        True,
        log_transform,
        1,
        score_cell_cycle,
        False,
    )

    if verbose:
        print("Embedding for pooling...")
    X_embed = _embed_for_pooling(
        adata_simple, dim_red_method_pooling, n_comps=embed_n_comps
    )
    del adata_simple

    if verbose:
        print("Pooling", str(X_embed.shape[0]), "cells...")
    _smooth_adata_by_pooling(adata, X_embed, n_neighbours=n_neighbors)

    prep_simple(
        adata,
        normalize_counts,
        filter_var_genes,
        n_top_genes,
        False,
        log_transform,
        1,
        score_cell_cycle,
        verbose,
    )

    if "scycle" not in adata.uns.keys():
        adata.uns["scycle"] = {}

    adata.uns["scycle"]["preprocess"] = {
        "method": "pooling",
        "n_neighbors": n_neighbors,
        "normalize_counts": normalize_counts,
        "filter_var_genes": filter_var_genes,
        "division_factor": division_factor,
        "log_transform": log_transform,
        "n_top_genes": n_top_genes,
        "embed_n_comps": embed_n_comps,
    }
    gc.collect()


def _embed_for_pooling(adata, dim_red, n_comps):
    if dim_red == "pca":
        sc.tl.pca(adata, n_comps=n_comps)
        X_embed = adata.obsm["X_pca"]
        return X_embed.copy()


def _smooth_adata_by_pooling(adata, X_embed, n_neighbours=10, copy=False):
    # adata_pooled = adata.copy if copy else adata
    nbrs = NearestNeighbors(n_neighbors=n_neighbours).fit(X_embed)
    distances, indices = nbrs.kneighbors(X_embed)
    adata.X = _smooth_matrix_by_pooling(_get_nd_array(adata.X), indices)
    if "matrix" in adata.layers:
        adata.layers["matrix"] = _smooth_matrix_by_pooling(
            _get_nd_array(adata.layers["matrix"]), indices
        )
    if "spliced" in adata.layers:
        adata.layers["spliced"] = _smooth_matrix_by_pooling(
            _get_nd_array(adata.layers["spliced"]), indices
        )
    if "unspliced" in adata.layers:
        adata.layers["unspliced"] = _smooth_matrix_by_pooling(
            _get_nd_array(adata.layers["unspliced"]), indices
        )

    adata.uns["scycle"] = True


def _smooth_matrix_by_pooling(matrix, indices):
    matrix_pooled = matrix.copy()
    for i in range(len(indices)):
        matrix_pooled[i, :] = np.mean(matrix[indices[i], :], axis=0)
    return matrix_pooled


def _get_nd_array(arr):
    if 'ndarray' not in str(type(arr)):
        x = arr.toarray()
    else:
        x = arr
    return x
