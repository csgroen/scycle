#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:50:39 2020

@author: clarice
"""
import numpy as np
from ._prep_simple import quality_control
from ._prep_pooling import prep_pooling
from ..tools import (
    dimensionality_reduction,
    enrich_components,
    principal_circle,
    pseudotime,
)


def normalize_by_partition(
    adata_src, adata_ref=None, rerun_pc=True, n_ref_parts=10, verbose=True
):

    """Normalize samples by median partition library size

     This procedure improves the discovery of the cell division moment.

    Parameters
    ------------
    adata_src: AnnData
        Unprocessed AnnData
    adata_ref: AnnData
        Reference AnnData object, that has been processed. If None, src will be
        automatically copied and preprocessed.
    rerun_pc: bool
        If True, `tl.principal circle` is re-run with n_ref_parts as n_nodes.
        Using fewer partitions for re-normalization than the default number used
    n_ref_parts: int
        How many partitions to be used for library size normalization.
    verbose: bool
        If True, the function will print messages.

    """
    if adata_ref is None:
        adata_ref = adata_src.copy()

    if "scycle" not in adata_ref.uns:
        print(
            "Reference adata should be preprocessed... Doing it with default values..."
        )
        prep_pooling(adata_ref, verbose=verbose)
        dimensionality_reduction(adata_ref, method="ica", verbose=verbose)

    pc_ran = "principal_circle" not in adata_ref.uns["scycle"].keys()

    if pc_ran | rerun_pc:
        if verbose:
            print("-- Running `tl.principal_circle` with n_ref_parts...")
        principal_circle(adata_ref, n_nodes=n_ref_parts, verbose=False)
        pseudotime(adata_ref, scale=False)

    # -- Run QC
    params = adata_ref.uns["scycle"]
    pp_params = params["preprocess"]
    dr_params = params["dimRed"]

    quality_control(
        adata_src,
        min_counts=pp_params["min_counts"],
        max_counts=pp_params["max_counts"],
        max_mt_ratio=pp_params["max_mt_ratio"],
        verbose=False,
    )
    old_totals = adata_src.obs["total_counts"]

    # --- Apply filter
    if verbose:
        print("Normalizing by partition...")

    # ---- Get partitions and re-noralize

    prt = adata_ref.obs["partition"]
    gexp = adata_src.X

    npart = np.max(prt) + 1
    new_gexp = np.empty(gexp.shape)
    # computing median per partition
    median_per_partition = np.array([0] * npart)
    for p in range(int(npart)):
        sidx = prt == p  # sample index
        totals = np.sum(gexp[sidx, :], axis=1)  # total counts per sample in group
        median_per_partition[p] = np.median(
            totals
        )  # median counts for samples in group

    totals = np.sum(gexp, axis=1).reshape(len(gexp), 1)
    next_prt = (prt + 1) % npart
    offsets = adata_ref.obs["pseudotime"] - adata_ref.obs["pseudotime"].astype(int)
    offsets = np.clip(offsets, 0, 1)
    medians = (1 - offsets) * median_per_partition[
        prt
    ] + offsets * median_per_partition[next_prt]

    adata_src.X = gexp / totals * np.array(medians).reshape(len(gexp), 1)

    adata_src.obs["total_counts"] = np.sum(adata_src.X, axis=1)
    adata_src.obs["total_counts_raw"] = old_totals
