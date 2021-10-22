#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
from anndata import AnnData
import numpy as np
from ..data import (
    g2m_markers,
    g1s_markers,
    histone_markers,
)
from scipy.stats import zscore
import gc

def prep_simple(
    adata: AnnData,
    normalize_counts: bool = True,
    filter_var_genes: bool = True,
    n_top_genes: int = 10000,
    for_pooling: bool = True,
    log_transform: bool = True,
    division_factor: float = 1,
    verbose: bool = True,
):
    """Pre-processes AnnData without pooling. Should be done only once.

    Parameters
    ----------
    adata: AnnData
        The raw AnnData object to be pre-processed
    normalize_counts: bool
        Set it to False if library does not need normalization
    filter_var_genes: bool
        If True, only `n_top_genes` highly variable genes are kept.
    n_top_genes: int
        Number of genes to keep after highly variable filter. Used if
        `filter_var_genes` is True. Passed to sc.pp.highly_variable_genes.
    for_pooling: bool
        Set to True if the function is called by the `prep_pooling` function.
        Changes the return object parameters.
    log_transform: bool
        Set it to false if you do not want values to be log-transformed.
    division_factor: int
        Scaling factor, divides the counts matrix by this value.
    verbose: bool
        If True, messages about function progress will be printed.

    Returns
    ----------
    - `None`
    """

    assert division_factor != 0, "Null division factor. Terminating..."
    if 'ndarray' not in str(type(adata.X)):
        adata.X = adata.X.toarray()
    # np.divide(adata.X, division_factor, out=adata.X)

    if "total_counts" not in adata.obs.keys():
        adata.obs["total_counts"] = adata.X.sum(1)

    # Normalization step
    if normalize_counts:
        sc.pp.normalize_total(adata, target_sum=np.median(adata.obs["total_counts"]))

    if verbose:
        print("Scoring cell cycle...")
    _score_cell_cycle(adata, g1s_markers, "G1-S")
    _score_cell_cycle(adata, g2m_markers, "G2-M")
    _score_cell_cycle(adata, histone_markers, "Histones")

    # Highly variable genes filtering
    if filter_var_genes:
        variances = np.var(adata.X, axis=0)
        inds = np.flip(np.argsort(variances))
        ind_genes = inds[0:n_top_genes]
        if 0 in variances[ind_genes]:
            ind_first_zero = np.argwhere(variances[ind_genes] == 0)[0][0]
            ind_genes = ind_genes[0:ind_first_zero]
        adata._inplace_subset_var(ind_genes)

    # Logarithmization
    if log_transform:
        sc.pp.log1p(adata, base=10)

    if not for_pooling:
        adata.uns["scycle"] = {
            "preprocess": {
                "method": "simple",
                "n_top_genes": n_top_genes,
                "normalize_counts": normalize_counts,
                "filter_var_genes": filter_var_genes,
                "division_factor": division_factor,
                "log_transform": log_transform,
                "n_top_genes": n_top_genes,
            }
        }
    gc.collect()


def _score_cell_cycle(adata, markers, score_name):
    sc.tl.score_genes(adata, markers, score_name=score_name)
    adata.obs[score_name] = zscore(adata.obs[score_name])
