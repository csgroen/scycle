#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scanpy as sc
import scrublet as scr
import numpy as np
from anndata import AnnData
import gc

def filter_cells(
    adata: AnnData,
    min_counts: int = -1,
    max_counts: int = -1,
    max_mt_ratio: int = 20,
    # doublet_detection: bool = False,
    # scrublet_kwargs: dict = {
    #     "total_counts": None,
    #     "sim_doublet_ratio": 2.0,
    #     "n_neighbors": None,
    #     "expected_doublet_rate": 0.1,
    #     "stdev_doublet_rate": 0.02,
    #     "random_state": 0,
    # },
    verbose=True,
):
    """Filter problematic cells in an AnnData

    Args:
      adata(AnnData): The AnnData object to be pre-processed.
      min_counts(int): Minimum number of counts required for a cell to pass filtering.
      `-1` -> median(counts) - std(counts)
      max_counts(int): Maximum number of counts required for a cell to pass filtering.
      `-1` -> median(counts) + std(counts)
      max_mt_ratio(int): Maximum proportion of mitochondrial genes in a cell to pass
      filtering.
      verbose: (Default value = True)

    Returns:
    * Sets
    """
    # doublet_detection(bool): Uses doublet detection instead of max counts to remove doublets
    # scrublet_kwargs(dict): Arguments passed to Scrublet for doublet detection

    # -- sparse -> array
    if 'ndarray' not in str(type(adata.X)):
        adata.X = adata.X.toarray()

    # -- Mitochondrial content
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # -- min/max suggestion
    counts = adata.X.sum(axis=1)
    md = np.median(counts)
    sd = np.std(counts)
    if min_counts == -1:
        min_counts = max(0, md - sd)
    if max_counts == -1:
        max_counts = md + sd
    # # -- Doublet detection
    # if doublet_detection:
    #     scrub = scr.Scrublet(
    #         adata.X,
    #         total_counts=scrublet_kwargs["total_counts"],
    #         sim_doublet_ratio=scrublet_kwargs["sim_doublet_ratio"],
    #         n_neighbors=scrublet_kwargs["n_neighbors"],
    #         expected_doublet_rate=scrublet_kwargs["expected_doublet_rate"],
    #         stdev_doublet_rate=scrublet_kwargs["stdev_doublet_rate"],
    #         random_state=scrublet_kwargs["random_state"],
    #     )
    #     (
    #         adata.obs["doublet_scores"],
    #         adata.obs["predicted_doublets"],
    #     ) = scrub.scrub_doublets()
    #     inds1 = np.where(
    #         (~adata.obs["predicted_doublets"].values)
    #         & (adata.obs["total_counts"] < max_counts)
    #         & (adata.obs["total_counts"] > min_counts)
    #     )
    #     del scrub
    # else:
    inds1 = np.where(
        (adata.obs["total_counts"] > min_counts) &
        (adata.obs["total_counts"] < max_counts))
    inds2 = np.where(adata.obs["pct_counts_mt"] < max_mt_ratio)
    if verbose:
        # if doublet_detection:
        #     print(np.sum(adata.obs["predicted_doublets"]), "doublets encountered")
        #     print(len(inds1[0]), "cells pass the doublet and counts filters.")
        # else:
        print(len(inds1[0]), "cells pass the count filter")
        print(len(inds2[0]), " cells pass the mt filter")
    ind_cells = np.intersect1d(inds1[0], inds2[0])
    if verbose:
        print("Cells selected", len(ind_cells))
    adata._inplace_subset_obs(ind_cells)
    gc.collect()
