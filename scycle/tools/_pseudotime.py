#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import elpigraph
from anndata import AnnData
from sklearn.mixture import GaussianMixture
def pseudotime(adata: AnnData, scale: bool=True, remap_border = False, border_threshold: float=0.05,
               verbose: bool=True):
    """ Calculate pseudotime from the principal circle
    
    Parameters
    ------------
    adata: AnnData
        The analysis object to be evaluated. Must be previously evaluated by
        `tl.principal_circle` and `tl.remap_nodes`.
    scale: bool
        If True, pseudotime values are given in a scale of 0 to 1.
    remap_border: bool
        If True, the points on the border of cell division will be re-evaluated
        based on their total counts to better predict if they have already
        divided or are about to divide
    border_threshold: float
        A number between 0 and 1 representing the percentages at each end that
        can be considered on the 'border' of cell division (i.e. right before 
        or just after cell division)
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
    X_counts = adata.obs['total_counts']
    
    if verbose: print('Calculating pseudotimes for each cell...')
    
    #-- Project onto graph
    ProjStruct = elpigraph.src.reporting.project_point_onto_graph(X_emb,
                                     NodePositions = node_p,
                                     Edges = edges,
                                     Partition = partition)
    edgeid = ProjStruct['EdgeID'].astype(int)
    prj_vals = ProjStruct['ProjectionValues']
    
    #-- Fix border effect
    if remap_border:
        edgeid = _remap_border_points(edgeid, X_counts, border_threshold, n_nodes)

    #-- Get time points
    node_times = np.array(range(n_nodes))
    
    time_point = np.array([_point_time(X_emb[i,:], prj_vals[i], edgeid[i], node_p, node_times, n_nodes) 
                           for i in range(len(edgeid))])
    # time_point = np.array([_point_time(X_emb[i,:], prj_vals[i], partition[i], node_p, node_times, n_nodes) for i in range(len(partition))])

    time_point[time_point < 0] = n_nodes - time_point[time_point < 0]
    time_point[time_point > n_nodes] = time_point[time_point > n_nodes] - n_nodes
    
    #-- Scale
    pseudotimes = np.array(time_point)/(n_nodes-1) if scale else np.array(time_point)

    #-- Add to adata
    adata.obs['pseudotime'] = pseudotimes
    adata.uns['scycle']['pseudotime'] = {'scale': scale}

def _remap_border_points(edgeid, X_counts, border_threshold, n_nodes):
    border_nnode = np.round(n_nodes * border_threshold)

    #-- Find points on the borders
    pts_lowerborder = edgeid < border_nnode
    pts_upperborder = edgeid >= n_nodes-border_nnode
    border_points = pts_lowerborder | pts_upperborder
    border_edgeids = edgeid[border_points]    
    border_counts = np.array(X_counts[border_points]).reshape(-1,1)
    
    #-- Define reference low counts after division / high counts before division
    after_div_refcnt = np.array(np.min(X_counts[pts_lowerborder])).reshape(-1,1)
    before_div_refcnt = np.array(np.max(X_counts[pts_upperborder])).reshape(-1,1)
    
    #-- Compute Gaussian Mixture model
    gmm = GaussianMixture(n_components=2).fit(border_counts)
    
    before_class = gmm.predict(before_div_refcnt)[0]
    after_class =  gmm.predict(after_div_refcnt)[0]
    
    celldiv_pred = gmm.predict(border_counts)
    remap1 = (border_edgeids >= n_nodes-border_nnode) & (celldiv_pred == after_class)
    border_edgeids[remap1] = 0
    
    remap2 = (border_edgeids < border_nnode) & (celldiv_pred == before_class)
    border_edgeids[remap2] = n_nodes-1
    
    edgeid[border_points] = border_edgeids
    
    return edgeid

    
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
        

        
