#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks


def cell_cycle_phase(adata,
                     smoothing_factor = 10,
                     verbose = True):
    """Estimates the phase of cell cycle for each cell based on the curvature
    of the trajectory embedded in the cell cycle space
    
    Parameters
    -----------
    adata: AnnData
        The analysis object to be evaluated. Must first be evaluated by
        `tl.pseudotime`.
    smoothing_factor: int
        A smoothing factor used for calculating the Univariate Spline used to
        estimate the curvature of the trajectory
    verbose: bool
        If True, the function will print messages.
        
    Returns
    ------------
    The `adata` object will be updated with the points of suggested start of 
    cell cycle phases. To assign cell cycle phase to each cell, you need to
    run `tl.annotate_cell_cycle`.
    """
    
    #-- Get node positions
    node_coords = adata.uns['princirc_gr']['node_coords']
    idx = np.array(['dim' in cname for cname in node_coords.columns])
    
    nodep = np.array(node_coords.iloc[:,idx])
    nnodes = nodep.shape[0]
    
    #-- Calculate curvature and find peaks
    x, curv = _calc_curvature(nodep, smoothing_factor)
    peaks = find_peaks(curv)[0]
    peak_times = peaks/nnodes
    
    #-- Use peaks to find phase divisions
    if len(peaks) == 3:
        s_start = peak_times[0]
        g2_start = peak_times[1]
        m_start = peak_times[2]
    if len(peaks) == 4:
        s_start = peak_times[1]
        g2_start = peak_times[2]
        m_start = peak_times[3]
    elif len(peaks) == 5:
        s_start = peak_times[1]
        g2_start = peak_times[3]
        m_start = peak_times[4]
    else:
        pkorder = np.argsort(curv[peaks])
        max_pktimes = [peak_times[pkorder[-1]], peak_times[pkorder[-2]]]
        s_start = np.min(max_pktimes) # one of the maximum peaks, earlier
        g2_start = np.max(max_pktimes) # one of the maximum peaks, later
        m_start = np.max(peak_times) # last peak

    #-- Save curvature info
    curv_data = pd.DataFrame(dict(x = x, pseudotime = x/nnodes, curvature = curv)).merge(
        pd.DataFrame(dict(x = peaks, ispeak = 'peak')),
        on = 'x', how = 'left')
    
    #-- Inform
    if verbose:
        print('-- Suggested cell cycle division:')
        print('G1:', ' 0  ', '-', s_start)
        print('S: ', s_start, '-', g2_start)
        print('G2:', g2_start, '-', m_start)
        print('M: ', m_start, '-   1')
    
    #-- Return
    adata.uns['scycle']['cell_cycle_division'] = {'s_start': s_start, 
                                                  'g2_start': g2_start, 
                                                  'm_start': m_start,
                                                  'curvature': curv_data}
    
    
def _calc_curvature(nodep,smoothing_factor=10):
    splines = []
    derivs = []
    for i in range(nodep.shape[1]):
        spline = UnivariateSpline(np.linspace(0,nodep.shape[0]-1,nodep.shape[0]), nodep[:,i],s=0,k=3)
        splines.append(spline)
        derivs.append(spline.derivative(n=2))
    n_points = nodep.shape[0]
    curv = np.zeros(nodep.shape[0])
    x = np.linspace(0,nodep.shape[0]-1,n_points)
    temp = np.zeros((n_points,nodep.shape[1]))
    print(temp.shape)
    for i in range(len(derivs)):
        temp[:,i] = derivs[i](x)
    curv = np.sqrt(np.sum(temp**2,axis=1))
    curv_spl = UnivariateSpline(x,curv,s=np.var(curv)*smoothing_factor,k=3)
    return x, curv_spl(x)

# def _calc_curvature_cyclic(nodep,smoothing_factor=5):
#     nodep3 = np.concatenate((nodep,nodep,nodep),axis=0)
#     x,curv = _calc_curvature(nodep3,smoothing_factor=smoothing_factor)
#     return x[nodep.shape[0]:2*nodep.shape[0]]-nodep.shape[0], curv[nodep.shape[0]:2*nodep.shape[0]]