#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional

import numpy as np
import elpigraph
from anndata import AnnData

def remap_nodes(adata: AnnData, start_node: Optional[list]= None,
                cycle_direction: Optional[int]=None,
                verbose = True):
    """ Remap principal circle nodes so that it starts at the moment of
    cell division

    Parameters
    -------------
    adata: AnnData
        The analysis object to be evaluated. Must be previously evaluated by
        `tl.principal_circle` and `tl.celldiv_moment` if `celldiv_edge` and
        `cycle_direction` are not provided.
    start_node: optional int
        The node that corresponds to start of cell cycle. If not
        provided, the node suggested by `tl.cell_division` will be used.
    cycle_direction: optional integer 1 or -1
        The direction of increase of the nodes according to the cell cycle
        (1 for clockwise, -1 for counter-clockwise). If not provided, the
        direction calculated by `tl.celldiv_moment` will be used.
    verbose: bool
        If True, the function will print messages.

    Returns
    ------------
    The principal circle coordinates will be updated so that the first node
    is the point immediately after cell division and the last node is the
    moment immediately preceeding it.
    """
    pc_gr = adata.uns['princirc_gr']
    edges = pc_gr['edges']
    node_coords = pc_gr['node_coords']
    edge_coords = pc_gr['edge_coords']
    X_emb = adata.obsm['X_cc'] if 'X_cc' in adata.obsm.keys() else adata.obsm['X_dimRed']
    n_dims = X_emb.shape[1]
    n_nodes = len(node_coords)

    #-- Get cell div edge/cycle direction
    if start_node == None:
        div_node = adata.uns['scycle']['cell_div_moment']['start_node']
    else:
        div_node = start_node

    if cycle_direction == None:
        cdir = adata.uns['scycle']['cell_div_moment']['cell_cycle_direction']
    else:
        cdir = cycle_direction

    #-- Remapping edges and nodes
    if verbose: print('Remapping edges using', div_edge, '...')

    #-- Get remapping positions
    start = start_node
    if cdir > 0:
        remap_vec = np.array(range(start, n_nodes+start))
        idx = remap_vec >= n_nodes
        remap_vec[idx] = remap_vec[idx] - n_nodes
    else:
        remap_vec = [start - el if el <= start else (n_nodes + start) - el for el in range(n_nodes)]
        # remap_edgec = [start - el if el <= start else (n_nodes + start) - el for el in range(n_nodes+1)]

    node_coords = _remap_node_coords(node_coords, remap_vec, n_nodes).reset_index(drop=True)
    node_p = np.array(node_coords.iloc[:,0:n_dims])
    # edge_coords = _remap_edge_coords(edge_coords, remap_edgec)
    edges = _remap_edges(edges, remap_vec).reset_index(drop=True)
    edge_tree = np.array(edges.iloc[1:(n_nodes),0:2])
    partition, dists = elpigraph.src.core.PartitionData(X = X_emb,
                                                        NodePositions = node_p,
                                                        MaxBlockSize = 100000000,
                                                        TrimmingRadius = np.inf,
                                                        SquaredX = np.sum(X_emb**2,axis=1,keepdims=1))

    adata.uns['princirc_gr'] = {'node_coords': node_coords, 'edge_coords': edge_coords,
                                'edges': edges, 'edge_tree': edge_tree}
    adata.uns['scycle']['cell_div_moment']['start_node'] = 0
    adata.uns['scycle']['cell_div_moment']['cell_cycle_direction'] = 1
    adata.obs['partition'] = partition

def _remap_node_coords(node_coords, remap_vec, n_nodes):
    node_coords = node_coords.iloc[remap_vec,:]
    node_coords.sort_values('npos')
    node_coords['npos'] = range(n_nodes)
    return node_coords


def _remap_edges(edges, remap_vec):
    #-- Reorder
    edges_new = edges.iloc[remap_vec,:]
    #-- Update e1
    edges_new['e1'] = range(len(remap_vec))
    #-- Update e2
    e2 = [i for i in range(1, len(remap_vec))]; e2.append(0)
    edges_new['e2'] = e2

    return edges_new
