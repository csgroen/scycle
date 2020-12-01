#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.linear_model import LinearRegression
from anndata import AnnData

def celldiv_moment (adata: AnnData, var = 'total_counts', verbose: bool= True):
    """Find moment of cell division
    
    Parameters
    ---------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.principal.circle`
    var: string
        Name of a variable in adata.obs used to compute the moment of cell division.
    verbose:
        If True, the function will print messages.
    
    Returns
    --------------
    `adata` will be modified with an estimation of the edge of the principal
    circle where cell division occurs and the direction of the cell cycle in
    the circle nodes (1 = clockwise, -1 = counter-clockwise). By default, moment
    of cell division is estimated as the edge with the maximum difference of 
    total counts between nodes; and the direction of cell cycle is estimated from
    the direction of the slope of total counts by nid.
    """
    pc_gr = adata.uns['princirc_gr']
    edges = pc_gr['edges'] 
    #-- Get var per node / diff_var 
    
    edges['mean_var'] = adata.obs.groupby('partition').mean()[var]
    diff_var = np.zeros(len(edges))
    for i in range(len(edges)-1):
        diff_var[i] = np.abs(edges['mean_var'][i+1] - edges['mean_var'][i])
        
    edges['diff_var'] = diff_var
        
    #-- Find moment of division (max abs diff)
    edge_to_max = np.argmax(edges['diff_var'])
    sugg_edge = edges.iloc[edge_to_max][['e1', 'e2']].values.astype(int)    
        
    if sugg_edge[0] != 0:
        coef1 = _edgeRegCoefficient(edges, end_edge = sugg_edge[0])
        coef2 = _edgeRegCoefficient(edges, start_edge = sugg_edge[0])
        
        direction = -1 if (coef1 < 0)[0] | (coef2 < 0)[0] else 1
    else:
        coef = _edgeRegCoefficient(edges)
        direction = 1 if (coef > 0)[0] else -1
    
    if verbose: 
        print('Suggested moment of cell division:', sugg_edge)
        print('Direction of cell cycle:', str(direction))
        
    #-- Add data to object
    adata.uns['scycle']['cell_div_moment'] = {'ref_var': var,
                                              'cell_div_edge': sugg_edge, 
                                              'cell_cycle_direction': direction}
    
def _edgeRegCoefficient(edges, start_edge = 0, end_edge = None):
    if end_edge is None:
        end_edge = len(edges)
    na_idx = np.invert(edges['mean_var'].apply(np.isnan))
    idx = (edges['e1'] < end_edge) & (edges['e1'] >= start_edge)
    idx = idx[na_idx]
    edges2 = edges[na_idx]
    
    x = np.array(edges2['e1'][idx]).reshape((-1,1))
    y = np.array(edges2['mean_var'][idx])
    model = LinearRegression().fit(x, y)
    
    return model.coef_
    