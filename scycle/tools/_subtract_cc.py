#!/usr/bin/env python3

import numpy as np
from ._principal_circle import principal_circle


def subtract_cc(adata):
    """
    Regresses out cell cycle influence with respect to trajectory.
    """
    if "partition" not in adata.obs:
        print(
            "-- Principal circle not not computed: computing it with default parameters."
        )
        principal_circle(adata)

    X = adata.X

    # Building mapping partition -> set of points
    partition = adata.obs.partition
    inds = [[] for _ in range(1 + np.max(partition))]
    for i in range(adata.n_obs):
        inds[partition[i]].append(i)

    # Computing barycenter of each cluster
    means = np.array([np.mean(X[ind, :], axis=0) for ind in inds])

    # Computing residues, subtracting cell cycle influence
    residue_matrix = means[partition, :]
    adata.obsm["residue_matrix"] = residue_matrix
    adata.X = X - residue_matrix
    adata.var["r2_scores"] = np.var(residue_matrix, axis=0) / np.var(X, axis=0)
