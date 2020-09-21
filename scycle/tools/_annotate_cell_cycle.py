#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional
import pandas as pd
from anndata import AnnData

def annotate_cell_cycle (adata: AnnData, s_start: Optional[float]=None, 
                         g2_start: Optional[float]=None, 
                         m_start: Optional[float]=None):
    """Annotates each cell with its cell cycle phase
    
    Parameters
    --------------
    adata: AnnData
         The analysis object to be evaluated. Must first be evaluated by
        `tl.cell_cycle_division`.
    s_start: Optional float
        The pseudotime of the start of the S-phase. If None, the S-phase start
        calculated by `tl.cell_cycle_division` is applied.
    g2_start: Optional float
        The pseudotime of the start of the G2 phase. If None, the G2 start
        calculated by `tl.cell_cycle_division` is applied.
    m_start: Optional float
        The pseudotime of the start of the M-phase. If None, the M-phase start
        calculated by `tl.cell_cycle_division` is applied.
    
    Returns
    -------------
    The `adata` object will be updated with assignment of cell cycle phases
    to each cell.
    """
    if ((s_start == None) | (g2_start == None) | (m_start == None)) & ('cell_cycle_division' not in adata.uns['scycle']):
        raise Exception('This object must be evaluated by `tl.cell_cycle_division` before `tl.annotate_cell_cycle`')
    
    if s_start == None:
        s_start = adata.uns['scycle']['cell_cycle_division']['sphase_start']
    if g2_start == None:
        g2_start = adata.uns['scycle']['cell_cycle_division']['g2_start']
    if m_start == None:
        m_start = adata.uns['scycle']['cell_cycle_division']['m_start']
        
    cell_times = adata.obs['pseudotime'].values
    cycle_phase = [_get_cycle_phase(cell_time, s_start, g2_start, m_start) for cell_time in cell_times]
    
    adata.obs['cell_cycle_phase'] = cycle_phase
    adata.obs['cell_cycle_phase'] = pd.Categorical(adata.obs['cell_cycle_phase'], 
                                                   categories = ['G1', 'S', 'G2', 'M'])

def _get_cycle_phase(cell_time, s_start, g2_start, m_start):
    if cell_time < s_start:
        return 'G1'
    elif (cell_time >= s_start) & (cell_time < g2_start):
        return 'S'
    elif (cell_time >= g2_start) & (cell_time < m_start):
        return 'G2'
    else:
        return 'M'  