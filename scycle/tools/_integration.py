#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from anndata import AnnData
from ot import dist, sinkhorn
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA


def integration(
    adata_source: AnnData,
    adata_reference: AnnData,
    components: list = [],
    metric: str = "sqeuclidean",
    eps=1e-3,
    new_obsm: str = "X_integrated",
    verbose: bool = False,
):
    """
    Optimal transport-based data integration. Data dimensionality must have been reduced
    before applying the method. Both source and reference must share the same cloud topology,
    otherwise the technique is very prone to overfitting.

    Parameters
    ----------
    adata_source: AnnData
        The dataset to align.
    adata_target: AnnData
        The dataset to use as a reference.
    obsm_key: str
        String identifier for accessing the low-dimension matrix.
        This is used as a AnnData.obsm[$obsm_key]
    metric: str, optional
        Metric to use for computing the cost matrix, passed as
        an argument to scipy.spatial.distance.cdist. Default is
        euclidean distance.
    eps: float, optional
        Error threshold for the sinkhorn algorithm. Tune it if
        encountering overflow/convergence issues.
    new_obsm: str, optional
        Key for source_adata.obsm into which store the integrated matrix.

    """
    assert (
        "dimRed" in adata_reference.uns
    ), "ICA projection matrix missing in source AnnData. \
        Dimensionality reduction must be computed first."

    if verbose:
        print("-- Integrating datasets")

    if "X_dimRed" not in adata_reference.obsm:
        if verbose:
            print("> Projecting reference dataset...")
        adata_reference.obsm["X_dimRed"] = adata_reference.uns["dimRed"].transform(
            adata_reference.X
        )

    print("> Projecting source dataset...")
    adata_source.obsm["X_dimRed"] = adata_reference.uns["dimRed"].transform(
        adata_source.X
    )

    Xs: np.ndarray = adata_source.obsm["X_dimRed"]
    Xt: np.ndarray = adata_reference.obsm["X_dimRed"]

    if len(components) > 0:
        Xs = Xs[:, components]
        Xt = Xt[:, components]

    adata_source.obsm["X_dimRed"] = _raw_ot_integration(
        Xs, Xt, metric=metric, eps=eps, verbose=verbose
    )

    pca = PCA(n_components=2, svd_solver="arpack")
    adata_source.obsm["X_dimRed2d"] = pca.fit_transform(adata_source.obsm["X_dimRed"])


def _raw_ot_integration(
    Xs: np.ndarray, Xt: np.ndarray, metric: str = "sqeuclidean", eps=1e-3, verbose=False
):
    """
    Unregularized optimal transport-based cell cycle integration. Maps
    Xs (source) points onto Xt (target).
    """
    n, d1 = Xs.shape
    m, d2 = Xt.shape

    assert (
        d1 == d2
    ), "Dimensions do not coincide %i (source) vs %i (target). Terminating..." % (
        d1,
        d2,
    )

    w_x, w_y = np.ones((n,)) / n, np.ones((m,)) / m
    print("> Computing distance matrix...")
    M = dist(Xs, Xt, metric=metric)
    M /= M.max()
    print("> Computing optimal transport plan...")
    Gs = sinkhorn(w_x, w_y, M, eps)
    print("> Projecting source dataset...")
    return np.diag(1 / w_x) @ Gs @ Xt
