#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from plotnine import aes, geom_smooth, geom_vline, annotate, labs
from ._scatter_pseudotime import scatter_pseudotime

def scatter_cell_cycle (adata, size = 1.5, alpha = 1, y_annot = 0):
    """Plots cell cycle signatures vs pseudotime
    
    Parameters
    ----------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.cell_cycle_phase`.
    size: float
        Controls the point size of the plot.
    alpha: float
        A value between 0 and 1. Controls point transparency.
    y_annot: float
        Controls the position of the cell cycle phase annotation.
    
    Returns
    --------------
    A plotnine scatter plot of pseudotime vs 3 cell cycle signatures.
    
    """
    if not 'cell_cycle_division' in adata.uns['scycle']:
        raise Exception('You need to run `tl.cell_cycle_phase` to compute the cell cycle phase divisions.')
    y = ['S-phase', 'G2-M', 'Histones']
    cc_divs = adata.uns['scycle']['cell_cycle_division']
    
    time_scatter = scatter_pseudotime(adata, y = y, size = size, alpha = alpha)
    
    #-- Add cell cycle annotations
    g1_labpos = np.mean((0, cc_divs['sphase_start']))
    s_labpos = np.mean((cc_divs['sphase_start'], cc_divs['g2_start']))
    g2_labpos = np.mean((cc_divs['g2_start'], cc_divs['m_start']))
    m_labpos = np.mean((cc_divs['m_start'], 1))
    labpos_df = pd.DataFrame({'label': ['G1', 'S', 'G2', 'M'],
                              'labpos': [g1_labpos, s_labpos, g2_labpos, m_labpos]})
        
    cell_cycle_plt = (time_scatter
    + geom_smooth(aes(group = 'signature'), method = 'loess', color = 'black', linetype = 'dashed')
    + geom_vline(xintercept = cc_divs['sphase_start'], linetype = 'dotted', size = 1)
    + geom_vline(xintercept = cc_divs['g2_start'], linetype = 'dotted', size = 1)
    + geom_vline(xintercept = cc_divs['m_start'], linetype = 'dotted', size = 1)
    + annotate('text', x = labpos_df['labpos'], y = y_annot, label = labpos_df['label'])
    + labs(x = 'Pseudotime', y = 'Signature scores', color = 'Signature')
    )

    return cell_cycle_plt