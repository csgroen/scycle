#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import elpigraph
from anndata import AnnData


def pseudotime(adata: AnnData, scale: bool=True, verbose: bool=True):
    """ Calculate pseudotime from the principal circle
    
    Parameters
    ------------
    adata: AnnData
        The analysis object to be evaluated. Must be previously evaluated by
        `tl.principal_circle` and `tl.remap_nodes`.
    scale: bool
        If True, pseudotime values are given in a scale of 0 to 1.
    verbose: bool
        If True, the function will print messages.

    """                
    #-- Get input data
    pr_gr = adata.uns['princirc_gr']
    X_emb = adata.obsm['X_dimRed2d']
    n_dims = X_emb.shape[1]
    node_p = np.array(pr_gr['node_coords'].iloc[:,0:n_dims])
    n_nodes = len(node_p)
    edges = np.array(pr_gr['edge_tree']); edges = np.vstack((edges, [0, n_nodes-1]))
    partition = adata.obs['partition']
    
    if verbose: print('Calculating pseudotimes for each cell...')
    
    #-- Project onto graph
    ProjStruct = elpigraph.src.reporting.project_point_onto_graph(X_emb,
                                     NodePositions = node_p,
                                     Edges = edges,
                                     Partition = partition)
    edgeid = ProjStruct['EdgeID'].astype(int)
    prj_vals = ProjStruct['ProjectionValues']        
    
    #-- Get time points
    node_times = np.array(range(n_nodes))
    
    time_point = np.array([_point_time(X_emb[i,:], prj_vals[i], edgeid[i], node_p, node_times, n_nodes) 
                           for i in range(len(edgeid))])
    # time_point = np.array([_point_time(X_emb[i,:], prj_vals[i], partition[i], node_p, node_times, n_nodes) for i in range(len(partition))])

    time_point[time_point < 0] = n_nodes - time_point[time_point < 0]
    time_point[time_point > n_nodes] = time_point[time_point > n_nodes] - n_nodes
    
    #-- Add to adata
    pseudotimes = np.array(time_point)/(n_nodes-1) if scale else np.array(time_point)
    adata.obs['pseudotime'] = pseudotimes
    adata.uns['scycle']['pseudotime'] = {'scale': scale}
    
def _point_time(Xi, prj_val, part, node_p, node_times, n_nodes):
    if prj_val < 0:
        prj_val_x = 0
    elif prj_val > 1:
        prj_val_x = 1
    else:
        prj_val_x = prj_val
    
    up = part+1 if part < n_nodes-1 else part
    down = part-1 if part > 0 else part
    
    dist_up = np.linalg.norm(node_p[up,:] - Xi)
    dist_down = np.linalg.norm(node_p[down,:] - Xi)
    
    if dist_up < dist_down:
        return node_times[part] + prj_val_x
    else:
        return node_times[part] - prj_val_x