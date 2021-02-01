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
    enrich_components,
    principal_circle,
    pseudotime,
)

def normalize_by_partition(
    adata_src: AnnData, 
    adata_ref: AnnData = None, 
    rerun_pc: bool = False, 
    n_ref_parts: int = 10,
    verbose: bool = True
):

    """Normalize samples by median partition library size

     This procedure improves the discovery of the cell division moment.

    Parameters
    ------------
    adata_src: AnnData
        Unprocessed AnnData
    adata_ref: AnnData
        Reference AnnData object, that has been processed.
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
        enrich_components(adata_ref, verbose=False)
    elif 'dimRed' not in adata_ref.uns['scycle'].keys():
        print("Preprocessing reference with default values...")
        dimensionality_reduction(adata_ref, method="ica", verbose=False)
        enrich_components(adata_ref, verbose=False)
    elif 'enrich_components' not in adata_ref.uns['scycle'].keys():
        if adata_ref.uns['scycle']['dimRed']['method'] == 'ica':
            print("Preprocessing reference with default values...")
            enrich_components(adata_ref, verbose=False)

    pc_didntrun = "principal_circle" not in adata_ref.uns["scycle"].keys()
    if pc_didntrun | rerun_pc:
        if verbose:
            print("-- Running `tl.principal_circle` with n_ref_parts...")
        principal_circle(adata_ref, n_nodes=n_ref_parts, verbose=False)
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
    adata_src.X = gexp / totals * part_meds
    
    new_totals = adata_src.layers['matrix'].sum(1)
    
    adata_src.obs["total_counts"] = new_totals
    adata_src.obs["total_counts_raw"] = old_totals
    