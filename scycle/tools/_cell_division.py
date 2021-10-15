#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.linear_model import LinearRegression
from anndata import AnnData
from ._remap_nodes import remap_nodes
import warnings

def cell_division (adata: AnnData,
                   var = 'total_counts',
                   dir_var = 'total_counts',
                   remap = True,
                   verbose: bool= True):
    """Find moment of cell division

    Parameters
    ---------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.principal.circle`
    var: string
        Name of a variable in adata.obs used to compute the moment of cell division.
        Generally, either 'total_counts' or 'G2-M'
    dir_var: string
        Name of a variable used to determine the direction of cell cycle.
    remap: bool
        If True, nodes will be remaped using the suggested edge
    verbose: bool
        If True, the function will print messages.

    Returns
    --------------
    `adata` will be modified with an estimation of the edge of the principal
    circle where cell division occurs and the direction of the cell cycle in
    the circle nodes (1 = clockwise, -1 = counter-clockwise).

    If var='total_counts' or another variable that countains transcript counts,
    moment of cell division is estimated as the edge with
    the maximum difference of total counts between nodes; and the direction of
    cell cycle is estimated from the direction of the slope of total counts by nid.

    If var='G2-M', moment of cell division is determined as the node where
    the decrease of the G2-M signatures starts.
    """
    pc_gr = adata.uns['princirc_gr']
    edges = pc_gr['edges']
    #-- Get var per node / diff_var
    edges['mean_var'] = adata.obs.groupby('partition').mean()[var]
    edges['mean_dirvar'] = adata.obs.groupby('partition').mean()[dir_var]

    #-- Get direction
    diff_dirvar = np.zeros(len(edges))
    for i in range(len(edges)-1):
        diff_dirvar[i] = np.abs(edges['mean_dirvar'][i+1] - edges['mean_dirvar'][i])
    edges['diff_dirvar'] = diff_dirvar
    edge_to_max = np.argmax(edges['diff_dirvar'])
    dirvar_edge = edges.iloc[edge_to_max][['e1', 'e2']].values.astype(int)

    coef = _edgeRegCoefficient(edges, var = 'mean_dirvar',
    start_edge = dirvar_edge[1], end_edge = dirvar_edge[0])
    direction = 1 if (coef > 0) else -1

    if var == 'G2-M':
        diff_var = np.zeros(len(edges))
        for i in range(len(edges)-1):
            diff_var[i] = edges['mean_var'][i+1] - edges['mean_var'][i]
        edges['diff_var'] = diff_var if direction == 1 else -diff_var

        div_edge = edges['diff_var'].idxmin()
        sugg_edge = edges.iloc[div_edge][['e1', 'e2']].values.astype(int)

    else:
        edges['mean_var'] = adata.obs.groupby('partition').mean()[var]
        diff_var = np.zeros(len(edges))
        for i in range(len(edges)-1):
            diff_var[i] = np.abs(edges['mean_var'][i+1] - edges['mean_var'][i])

        edges['diff_var'] = diff_var

        #-- Find moment of division (max abs diff)
        edge_to_max = np.argmax(edges['diff_var'])
        sugg_edge = edges.iloc[edge_to_max][['e1', 'e2']].values.astype(int)

    if verbose:
        print('Suggested moment of cell division:', sugg_edge)
        print('Direction of cell cycle:', str(direction))

    #-- Add data to object
    adata.uns['scycle']['cell_div_moment'] = {'ref_var': var,
                                              'cell_div_edge': sugg_edge,
                                              'cell_cycle_direction': direction}
    if remap:
        remap_nodes(adata, verbose = verbose)

def _edgeRegCoefficient(edges, var = 'mean_counts', start_edge = 0, end_edge = None):
    nedge = edges.shape[0]
    x = np.array([float(i) for i in range(nedge)])
    y = []

    for i in range(nedge):
        edge1 = start_edge + i if start_edge + i < nedge else i - (nedge-start_edge)
        idx = edges['e1'] == edge1
        y.append(edges[idx][var].values[0])

    na_idx = ~np.isnan(y)
    x = x[na_idx]
    y = np.array(y)[na_idx]
    model = LinearRegression().fit(x.reshape(nedge,-1),y)

    return model.coef_[0]

def celldiv_moment(adata: AnnData, method = 'counts', var = 'total_counts', remap = True, verbose: bool= True):
    """DEPRECATED: use cell_division
    """
    warnings.warn("celldiv_moment is deprecated; use cell_division", DeprecationWarning)
    return (cell_division(adata, method, var, remap, verbose))
