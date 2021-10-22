#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional
import pandas as pd
from anndata import AnnData

def annotate_cell_cycle (adata: AnnData, pr_start: Optional[float]=None,
                         rep_start: Optional[float]=None):
                         # m_start: Optional[float]=None):
    """Annotates each cell with its cell cycle phase

    Parameters
    --------------
    adata: AnnData
         The analysis object to be evaluated. Must first be evaluated by
        `tl.cell_cycle_division`.
    pr_start: Optional float
        The pseudotime of the start of the S-phase. If None, the S-phase start
        calculated by `tl.cell_cycle_division` is applied.
    rep_start: Optional float
        The pseudotime of the start of the G2 phase. If None, the G2 start
        calculated by `tl.cell_cycle_division` is applied.

    Returns
    -------------
    The `adata` object will be updated with assignment of cell cycle phases
    to each cell.
    """
    if ((pr_start == None) | (rep_start == None)) & ('cell_cycle_division' not in adata.uns['scycle']):
        raise Exception('This object must be evaluated by `tl.cell_cycle_division` before `tl.annotate_cell_cycle`')

    if pr_start == None:
        pr_start = adata.uns['scycle']['cell_cycle_division']['pr_start']
    else:
        adata.uns['scycle']['cell_cycle_division']['pr_start'] = pr_start
    if rep_start == None:
        rep_start = adata.uns['scycle']['cell_cycle_division']['rep_start']
    else:
        adata.uns['scycle']['cell_cycle_division']['rep_start'] =  rep_start
    # if m_start == None:
    #     m_start = adata.uns['scycle']['cell_cycle_division']['m_start']
    # else:
    #     adata.uns['scycle']['cell_cycle_division']['m_start'] = m_start

    cell_times = adata.obs['pseudotime'].values
    cycle_phase = [_get_cycle_phase(cell_time, pr_start, rep_start) for cell_time in cell_times]
    # cycle_phase = [_get_cycle_phase(cell_time, pr_start, rep_start, m_start) for cell_time in cell_times]

    adata.obs['cell_cycle_phase'] = cycle_phase
    adata.obs['cell_cycle_phase'] = pd.Categorical(adata.obs['cell_cycle_phase'],
                                                   categories = ['G1 post-mitotic', 'G1 pre-replication', 'S/G2/M'])

def _get_cycle_phase(cell_time, pr_start, rep_start):
    if cell_time < pr_start:
        return 'G1 post-mitotic'
    elif (cell_time >= pr_start) & (cell_time < rep_start):
        return 'G1 pre-replication'
    # elif (cell_time >= rep_start) & (cell_time < m_start):
    #     return 'G2'
    else:
        return 'S/G2/M'
        # return 'M'
