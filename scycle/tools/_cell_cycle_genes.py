#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from ._subtract_cc import subtract_cc

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
    Xtemp = adata.X
    cc.tl.subtract_cc(adata)
    adata.X = Xtemp
    r2scores = adata.var['r2_scores']

    r2_filtered = r2scores[r2scores > 0.5].sort_values(ascending = False)
    r2_df = pd.DataFrame({'genes': r2_filt.index.tolist(), 'r2_cc_trajectory': r2_filt})

    return(r2_df)
