#!/usr/bin/env python3

import numpy as np
from ._trajectory import trajectory

def subtract_cc(adata):
    """Regresses out cell cycle influence with respect to trajectory.

    Parameters
    ------------
    adata: AnnData
    """
    if "partition" not in adata.obs:
        print(
            "-- Trajectory not not computed: computing it with default parameters."
        )
        trajectory(adata)

    X = adata.X
    partition = adata.obs.partition
    
    residue_matrix, r2_scores = _cc_residue_matrix(X, partition)

    adata.obsm["residue_matrix"] = residue_matrix
    adata.X = X - residue_matrix
    adata.var["r2_scores"] = r2_scores

def _cc_residue_matrix(X, partition):
    # Building mapping partition -> set of points
    inds = [[] for _ in range(1 + np.max(partition))]
    for i in range(X.shape[0]):
        inds[partition[i]].append(i)
    if any(len(ind) == 0 for ind in inds):
        print("Warning: empty partitions.")

    # Computing barycenter of each cluster
    means = np.array(
        [
            (np.mean(X[ind, :], axis=0) if len(ind) > 0 else np.zeros((X.shape[1],)))
            for ind in inds
        ]
    )

    # Computing residues, subtracting cell cycle influence
    residue_matrix = means[partition, :]
    r2_scores = np.var(residue_matrix, axis=0) / np.var(X, axis=0)
    
    return((residue_matrix, r2_scores))
    