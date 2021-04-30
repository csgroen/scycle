#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:50:39 2020

@author: clarice
"""
from anndata import AnnData
import numpy as np
from ._prep_pooling import prep_pooling
from ..tools import (
    dimensionality_reduction,
    find_cc_components,
    self_consistent_trajectory,
    pseudotime
)

def normalize_by_partition(
    adata_src: AnnData, 
    adata_ref: AnnData = None, 
    rerun_pc: bool = False, 
    n_ref_parts: int = 20,
    verbose: bool = True
):

    """Normalize samples by median partition library size

     This procedure improves the discovery of the cell division node.

    Parameters
    ------------
    adata_src: AnnData
        AnnData object that has been through quality control and pre-processing
    adata_ref: AnnData
        Reference AnnData object. If not provided, trajectory and partitions are
        calculated with default values from adata_src.
    rerun_pc: bool
        If True, `tl.principal circle` is re-run with n_ref_parts as n_nodes
    n_ref_parts: int
        How many partitions to be used for library size normalization.
    verbose: bool
        If True, the function will print messages.

    """
    if adata_ref is None:
        adata_ref = adata_src.copy()
    else:
        adata_ref = adata_ref.copy()

    if "scycle" not in adata_ref.uns:
        print("Preprocessing reference with default values...")
        prep_pooling(adata_ref, verbose=False)
        dimensionality_reduction(adata_ref, method="ica", verbose=False)
        find_cc_components(adata_ref, verbose=False)
    elif 'dimRed' not in adata_ref.uns['scycle'].keys():
        print("Preprocessing reference with default values...")
        dimensionality_reduction(adata_ref, method="ica", verbose=False)
        find_cc_components(adata_ref, verbose=False)
    elif 'find_cc_components' not in adata_ref.uns['scycle'].keys():
        if adata_ref.uns['scycle']['dimRed']['method'] == 'ica':
            print("Preprocessing reference with default values...")
            find_cc_components(adata_ref, verbose=False)

    sct_didntrun = 'self-consistent_trajectory' not in adata_ref.uns["scycle"].keys()
    if sct_didntrun | rerun_pc:
        if verbose:
            print("-- Running `tl.self_consistent_trajectory` with n_ref_parts...")
        self_consistent_trajectory(adata_ref, n_nodes=n_ref_parts, verbose=False)
        pseudotime(adata_ref, scale = False, verbose = False)

    old_totals = adata_ref.obs["total_counts"]

    # ---- Get partitions and re-normalize
    prt = adata_ref.obs["partition"]
    
    gexp_raw = adata_src.layers['matrix'] # "raw" data
    gexp = adata_src.X

    npart = np.max(prt) + 1
    # computing median per partition
    median_per_partition = np.zeros(npart)
    rmedian_per_partition = np.zeros(npart)
    for p in range(int(npart)):
        sidx = prt == p  # sample index
        totals = np.sum(gexp[sidx, :], axis=1)  # total counts per sample in group
        median_per_partition[p] = np.median(
            totals
        )  # median counts for samples in group
        
        rtotals = gexp_raw[sidx,:].sum(1)
        rmedian_per_partition[p] = np.median(rtotals)

    totals = np.sum(gexp, axis=1).reshape(len(gexp), 1)
    rtotals = gexp_raw.sum(1).reshape(len(gexp),1)
    next_prt = (prt + 1) % npart
    offsets = adata_ref.obs["pseudotime"] - adata_ref.obs["pseudotime"].astype(int)
    offsets = np.clip(offsets, 0, 1)
    medians = (1 - offsets) * median_per_partition[
        prt
    ] + offsets * median_per_partition[next_prt]
    rmedians = (1 - offsets) * rmedian_per_partition[
        prt
    ] + offsets * rmedian_per_partition[next_prt]
    
    #-- Normalize counts
    if verbose: print("Re-normalizing counts by partition...")
    part_meds = np.array(medians).reshape(len(gexp), 1)
    part_rmeds = np.array(rmedians).reshape(len(gexp), 1)
    adata_src.layers['matrix'] = gexp_raw / rtotals * part_rmeds
    adata_src.layers['unnorm'] = gexp
    adata_src.X = gexp / totals * part_meds
    
    new_totals = adata_src.layers['matrix'].sum(1)
    
    adata_src.obs["total_counts"] = new_totals
    adata_src.obs["total_counts_raw"] = old_totals
    