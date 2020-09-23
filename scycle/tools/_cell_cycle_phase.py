#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from skmisc.loess import loess

def cell_cycle_phase(adata, n_points = 100, 
                     ref_gradient = 0.01, 
                     plateau_score_minquant = 0.65, 
                     verbose = True):
    """Estimates the phase of cell cycle for each cell based on pseudotime
    
    Parameters
    -----------
    adata: AnnData
        The analysis object to be evaluated. Must first be evaluated by
        `tl.pseudotime`.
    n_points: int
        Number of points to use for the local regressions used to estimate
        the signatures.
    ref_gradient: float
        ref_gradient controls the mimimal gradient for the slope of the S-phase
        signature for start of S-phase and the gradient below which the histones
        plateau starts for passing from S-phase into G2. The higher the number
        of nodes, the smaller this number should be.
    plateau_score_minquant: float
        Minimal quantile of histones score for passage of S-phase to G2.
    verbose: bool
        If True, the function will print messages.
        
    Returns
    ------------
    The `adata` object will be updated with the points of suggested start of 
    cell cycle phases. To assign cell cycle phase to each cell, you need to
    run `tl.annotate_cell_cycle`.
    """
        
    #-- Divide cell cycle
    pseudotime = adata.obs['pseudotime']
    sphase = adata.obs['S-phase']
    hists = adata.obs['Histones']     
    g2m = adata.obs['G2-M']
    
    #-- Get y's (time points)
    time_points = np.array(range(0,n_points))/n_points
    
    if not adata.uns['scycle']['pseudotime']['scale']:
        n_nodes = adata.uns['scycle']['principal_circle']['n_nodes']
        time_points = time_points * n_nodes
    
    #-- Get loess models    
    sphase_pred = _signature_loess(pseudotime, sphase, time_points)
    hist_pred = _signature_loess(pseudotime, hists, time_points)
    g2m_pred = _signature_loess(pseudotime, g2m, time_points)
    
    #-- Get gradient changes
    sphase_grad = np.gradient(sphase_pred)
    hist_grad = np.gradient(hist_pred)
    g2m_grad = np.gradient(g2m_pred)
    
    #-- Find limits of phases
    sphase_start = np.min(np.where(sphase_grad > ref_gradient))
    
    hist_plateau = (hist_grad < ref_gradient) & [i > sphase_start for i in range(n_points)] & (hist_pred > np.quantile(hist_pred, plateau_score_minquant))
    g2_start = np.min(np.where(hist_plateau))
    
    m_increase = [i > g2_start for i in range(n_points)] & (g2m_grad < ref_gradient)
    m_start = np.min(np.where(m_increase)) 
    
    #-- Inform
    if verbose:
        print('-- Suggested cell cycle division:')
        print('G1:', ' 0  ', '-', time_points[sphase_start])
        print('S: ', time_points[sphase_start], '-', time_points[g2_start])
        print('G2:', time_points[g2_start], '-', time_points[m_start])
        print('M: ', time_points[m_start], '-   1')
    
    #-- Return
    adata.uns['scycle']['cell_cycle_division'] = {'sphase_start': time_points[sphase_start], 
                                                  'g2_start': time_points[g2_start], 
                                                  'm_start': time_points[m_start]}
    
def _signature_loess(x, y, xpreds):
    model = loess(x,y)
    preds = model.predict(xpreds).values
    return preds
