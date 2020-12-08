#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import elpigraph
from anndata import AnnData
from sklearn.decomposition import PCA


def principal_circle(adata: AnnData, n_nodes: int = 30, verbose: bool = True):
    """Calculates the principal circle nodes and edges for estimation of
    cell-cycle pseudotime

    Parameters
    ------------
    adata: AnnData
        The analysis object to be evaluated. Must first be evaluated by
        `tl.dimensionality_reduction`.
    n_nodes: int
        Number of nodes to be used to calculate the principal circle.
    verbose: bool
        If True, the function will print messages.

    Returns
    ------------
    `adata` will be updated with the coordinates of the nodes and edges of the
    principal circle.
    """
    X_emb = (
        adata.obsm["X_cc"] if "X_cc" in adata.obsm.keys() else adata.obsm["X_dimRed"]
    )
    n_dims = X_emb.shape[1]

    egr = elpigraph.computeElasticPrincipalCircle(
        X_emb, NumNodes=n_nodes, verbose=verbose
    )
    adata.uns["egr"] = egr[0]

    # -- Get graph data
    node_coords, edge_coords = _get_gr_coords(adata)
    node_coords = node_coords.sort_values("npos")
    node_p = np.array(node_coords.iloc[:, 0:n_dims])

    # -- Make edgelist (for npos)
    e1 = [i for i in range(n_nodes)]
    e2 = [i + 1 for i in range(n_nodes)]
    e2[n_nodes - 1] = 0

    # -- Map cells to nodes
    partition, dists = elpigraph.src.core.PartitionData(
        X=X_emb,
        NodePositions=node_p,
        MaxBlockSize=100000000,
        TrimmingRadius=np.inf,
        SquaredX=np.sum(X_emb ** 2, axis=1, keepdims=1),
    )
    adata.obs["partition"] = partition

    # -- Count per node
    total_counts = adata.obs["total_counts"].to_numpy()

    node_read_counts = [
        np.mean(total_counts[np.where(partition == i)[0]]) for i in range(len(node_p))
    ]
    node_coords["total_counts"] = node_read_counts

    edge_data = pd.DataFrame({"e1": e1, "e2": e2, "mean_counts": node_read_counts})

    # -- Project in 3d
    pca = PCA(n_components=3).fit(X_emb)
    node3d = pd.DataFrame(pca.transform(node_coords.iloc[:, 0:n_dims]))
    node3d.columns = ["x", "y", "z"]
    node3d["npos"] = range(n_nodes)
    edge3d = pd.DataFrame(pca.transform(edge_coords.iloc[:, 2 : (n_dims + 2)]))
    edge3d.columns = ["x", "y", "z"]
    edge3d["eid"] = range(n_nodes + 1)

    node_coords = node_coords.merge(node3d, how="left", on="npos")
    edge_coords = edge_coords.merge(edge3d, how="left", on="eid")

    # -- Add to adata
    adata.uns["princirc_gr"] = {
        "node_coords": node_coords,
        "edge_coords": edge_coords,
        "edges": edge_data,
    }
    adata.uns["scycle"]["principal_circle"] = {"n_nodes": n_nodes}


def _get_gr_coords(adata):
    # Get the node coordinates
    node_p = adata.uns["egr"]["NodePositions"]
    node_coords = pd.DataFrame(node_p)
    node_coords.columns = ["dim" + str(i) for i in node_coords.columns]
    node_coords["nid"] = range(len(node_coords.index))
    # Get edges
    edges = pd.DataFrame(adata.uns["egr"]["Edges"][0], columns=["e1", "e2"])
    # Get edge coordinates
    edge_coords, node_order = _get_edge_coords(node_coords, edges)
    edges = adata.uns["egr"]["Edges"][0]
    node_order = node_order[0 : len(node_order) - 1]
    node_ord_df = pd.DataFrame({"nid": node_order, "npos": range(len(node_order))})
    # -- Add node positions
    node_coords = node_coords.merge(node_ord_df, how="left", on="nid")
    node_coords = node_coords.sort_values("npos")

    return node_coords, edge_coords


def _get_edge_coords(point_coords, edges):
    edges = edges.sort_values("e1")
    start_pid = edges["e1"][0]
    s_order = [start_pid]
    pid = None

    while pid != start_pid:
        pid = s_order[len(s_order) - 1]
        query_res = (
            edges.query("e1 == " + str(pid) + "| e2 ==" + str(pid)).to_numpy().flatten()
        )
        connects = query_res[query_res != pid]
        if len(s_order) == 1:
            s_order.append(connects[0])
            pid = connects[0]
        else:
            idx = [not (c in s_order) for c in connects]
            if sum(idx) == 0:
                s_order.append(start_pid)
                pid = start_pid
            else:
                s_order.append(connects[idx][0])

    edge_coords = pd.DataFrame({"eid": range(len(s_order)), "nid": s_order})
    edge_coords = edge_coords.merge(
        point_coords, how="left", left_on="nid", right_on="nid"
    ).sort_values("eid")
    node_order = s_order

    return edge_coords, node_order
