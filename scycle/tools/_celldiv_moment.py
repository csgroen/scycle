#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.linear_model import LinearRegression
from anndata import AnnData

def celldiv_moment (adata: AnnData, verbose: bool= True):
    """Find moment of cell division
    
    Parameters
    ---------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    verbose:
        If True, the function will print messages.
    
    Returns
    --------------
    `adata` will be modified with an estimation of the edge of the principal
    circle where cell division occurs and the direction of the cell cycle in
    the circle nodes (1 = clockwise, -1 = counter-clockwise)
    """
    pc_gr = adata.uns['princirc_gr']
    edges = pc_gr['edges'] 
    #-- Find moment of division (max abs diff)
    edge_to_max = np.argmax(edges['dif_counts'])
    sugg_edge = edges.iloc[edge_to_max][['e1', 'e2']].values.astype(int)
    
    if sugg_edge[0] != 0:
        coef1 = _edgeRegCoefficient(edges, end_edge = sugg_edge[0])
        coef2 = _edgeRegCoefficient(edges, start_edge = sugg_edge[0])
        
        direction = -1 if (coef1 < 0)[0] & (coef2 < 0)[0] else 1
    else:
        coef = _edgeRegCoefficient(edges)
        direction = 1 if (coef > 0)[0] else -1
    
    if verbose: 
        print('Suggested moment of cell division:', sugg_edge)
        print('Direction of cell cycle:', str(direction))
        
    #-- Add data to object
    adata.uns['scycle']['cell_div_moment'] = {'cell_div_edge': sugg_edge, 
                                              'cell_cycle_direction': direction}
    
def _edgeRegCoefficient(edges, start_edge = 0, end_edge = None):
    if end_edge is None:
        end_edge = len(edges)
    idx = (edges['e1'] < end_edge) & (edges['e1'] >= start_edge)
    
    x = np.array(edges['e1'][idx]).reshape((-1,1))
    y = np.array(edges['mean_counts'][idx])
    model = LinearRegression().fit(x, y)
    
    return model.coef_
    
    