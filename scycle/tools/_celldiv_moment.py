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
      
    #-- Find direction of cell division
    x = np.array(edges['e1']).reshape((-1,1))
    y = np.array(edges['mean_counts'])
    model = LinearRegression().fit(x, y)
    direction = 1 if (model.coef_ > 0)[0] else -1
    
    if verbose: 
        print('Suggested moment of cell division:', sugg_edge)
        print('Direction of cell cycle:', str(direction))
        
    #-- Add data to object
    adata.uns['scycle']['cell_div_moment'] = {'cell_div_edge': sugg_edge, 
                                              'cell_cycle_direction': direction}