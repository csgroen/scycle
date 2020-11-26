#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
from anndata import AnnData
import numpy as np


def prep_simple(
    adata: AnnData,
    filter_cells: bool = True,
    min_counts: int = 10000,
    max_counts: int = 40000,
    max_mt_ratio: int = 20,
    normalize_counts: bool = True,
    target_sum: int = 10000,
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
        The raw AnnData object to be pre-processed.
    filter_cells: bool
        Set it to False if bad quality cells were already filtered
    min_counts: int
        Minimum number of counts required for a cell to pass filtering.
    max_counts: int
        Maximum number of counts required for a cell to pass filtering.
    max_mt_ratio: int
        Maximum proportion of mitochondrial genes in a cell to pass
        filtering.
    normalize_counts: bool
        Set it to False if library does not need normalization
    target_sum: int
        Target sum of counts for library normalization. Passed to
        sc.pp.normalize_total.
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
    None
    """

    assert division_factor != 0, "Null division factor. Terminating..."
    adata.X = adata.X / division_factor

    # Quality control
    if filter_cells:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        inds1 = np.where(
            (adata.obs["total_counts"] > min_counts)
            & (adata.obs["total_counts"] < max_counts)
        )
        inds2 = np.where(adata.obs["pct_counts_mt"] < max_mt_ratio)
        if verbose:
            print(len(inds1[0]), "samples pass the count filter")
            print(len(inds2[0]), " samples pass the mt filter")
        ind_samples = np.intersect1d(inds1[0], inds2[0])
        if verbose:
            print("Samples selected", len(ind_samples))
        adata._inplace_subset_obs(ind_samples)

    # Normalization step
    if normalize_counts:
        sc.pp.normalize_total(adata, target_sum=target_sum)

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
        sc.pp.log1p(adata)

    if not for_pooling:
        adata.uns["scycle"] = {
            "preprocess": {
                "method": "simple",
                "min_counts": min_counts,
                "target_sum": target_sum,
                "filter_var_genes": filter_var_genes,
                "n_top_genes": n_top_genes,
            }
        }
