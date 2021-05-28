#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from ._subtract_cc import subtract_cc

#-- test
import anndata

adata = anndata.read_h5ad('/home/clarice/Desktop/test.h5ad')
#---

def cell_cycle_genes(adata, r2_thrshold = 0.6):
    """Uses the cell cycle trajectory to shortlist cell cycle genes

    Parameters
    ------------
    adata: AnnData
        The analysis object to be evaluated. Must first be evaluated by
        `tl.pseudotime`.
    r2_threshold: float
        R^2 threshold for cell-cycle related genes when regressed against the
        trajectory.

    Returns
    ------------
    `adata` will be updated with the partitions calculated from the self-consistent
    trajectory
    """
    assert 'egr' in adata.uns, "Pseudotime must be computed first."

    Xtemp = adata.X
    cc.tl.subtract_cc(adata)
    adata.X = Xtemp
    r2scores = adata.var['r2_scores']

    r2_filtered = r2scores[r2scores > 0.5].sort_values(ascending = False)
    r2_df = pd.DataFrame({'genes': r2_filt.index.tolist(), 'r2_cc_trajectory': r2_filt})

    return(r2_df)

def cc_comp_genes(adata, nsd = 3):
    """Returns annotation for the genes used to compute the cell cycle space

    Parameters
    ------------
    adata: AnnData
        The analysis object to be evaluated. Must first have been evaluated by
        `tl.find_cc_components`
    nsd: int
        Minimum number of standard deviations of the gene weight to be included
        in the annotation. The higher this number, the more stringent the annotation.

    Returns
    ------------
    A pandas.DataFrame containing the information on the genes used to compute the
    cell cycle components
    """

    #-- Get
    assert (
    "find_cc_components" in adata.uns['scycle']
    ), "Cell cycle components missing.\
    Dimensionality reduction and finding cell cycle components must be computed first."

    cc_comps = adata.uns['scycle']['find_cc_components']['indices']
    idx = [int(i) for i in cc_comps.values()]
    cc_Pmat = adata.varm['P_dimRed'][:,idx]

    genes = np.array(adata.var.index)

    cc_gene_list = []
    cc_weight_list = []
    cc_gene_dir = []
    cc_comp_list = []
    sds = np.std(cc_Pmat, axis = 0)
    thr = sds * nsd


    #-- Ugly loop
    for i in range(len(cc_comps.keys())):
        ic_name = list(cc_comps.keys())[i]
        c_weights = cc_Pmat[:,i]
        idx_pos = c_weights > thr[i]
        idx_neg = c_weights < -thr[i]

        #-- List extends
        plen = np.sum(idx_pos); nlen = np.sum(idx_neg)
        cc_gene_list.extend(genes[idx_pos])
        cc_gene_list.extend(genes[idx_neg])

        cc_weight_list.extend(c_weights[idx_pos])
        cc_weight_list.extend(c_weights[idx_neg])

        cc_gene_dir.extend([1] * plen)
        cc_gene_dir.extend([-1] * nlen)

        cc_comp_list.extend([ic_name]*(plen+nlen))

        res_tb = pd.DataFrame({
        'component': cc_comp_list,
        'gene': cc_gene_list,
        'direction': cc_gene_dir,
        'weight': cc_weight_list
        })

    res_tb = res_tb.sort_values('component').sort_values('direction').sort_values('weight', ascending=False).reset_index(drop = True)
    res_tb

    return res_tb
