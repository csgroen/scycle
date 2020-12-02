#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from anndata import AnnData
from ot import dist, sinkhorn
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA

from ._enrich_components import enrich_components


def integration(
    _adata_src: AnnData,
    _adata_ref: AnnData,
    components: list = [],
    metric: str = "sqeuclidean",
    eps=1e-3,
    verbose: bool = False,
):
    """
    Optimal transport-based data integration. Data dimensionality must have been reduced
    before applying the method. Both source and reference must share the same cloud topology,
    otherwise the technique is very prone to overfitting.

    Parameters
    ----------
    adata_src: AnnData
        The dataset to align.
    adata_target: AnnData
        The dataset to use as a reference.
    components: list
        ICs indices caring cell-cycle related information. Let
        it empty for automatic detection.
    metric: str, optional
        Metric to use for computing the cost matrix, passed as
        an argument to scipy.spatial.distance.cdist. Default is
        euclidean distance.
    eps: float, optional
        Error threshold for the sinkhorn algorithm. Tune it if
        encountering overflow/convergence issues.
    """
    assert (
        "dimRed" in _adata_ref.uns
    ), "ICA projection matrix missing in source AnnData. \
        Dimensionality reduction must be computed first."

    # Selecting the cell-cycle ICs components
    if len(components) == 0:
        if verbose:
            print("-- Automatically detecting cell-cycle components...")
        if 'enrich_components' not in _adata_ref.uns['scycle'].keys():
            enrich_components(_adata_ref)
        components = list(_adata_ref.uns["scycle"]["enrich_components"].values())

    if verbose:
        print("-- Integrating datasets...")

    # Filtering only common genes
    if verbose:
        print("> Selecting common genes...")
    src_genes = _adata_src.var_names
    ref_genes = _adata_ref.var_names
    common_genes = [g for g in ref_genes if g in src_genes]
    if verbose:
        print("> %i genes selected." % len(common_genes))
        print("> Slicing matrices...")
    genes_idx, idx = [], 0
    for i, g in enumerate(ref_genes):
        if common_genes[idx] == g:
            genes_idx.append(i)
            idx += 1
    adata_ref = _adata_ref[:, common_genes]
    adata_src = _adata_src[:, common_genes]

    # Projecting both datasets in the same ICs subspace
    if "X_dimRed" not in adata_ref.obsm:
        if verbose:
            print("> Projecting reference dataset...")
        Y_c = adata_ref.X.copy()
        Y_c -= np.mean(Y_c, axis=0)
        adata_ref.obsm["X_dimRed"] = Y_c @ (adata_ref.uns["P_dimRed"].T)[genes_idx, :]

    X_c = adata_src.X.copy()
    X_c -= np.mean(X_c, axis=0)
    adata_src.obsm["X_dimRed"] = X_c @ (adata_ref.uns["P_dimRed"].T)[genes_idx, :]

    Xs: np.ndarray = adata_src.obsm["X_dimRed"]
    Xt: np.ndarray = adata_ref.obsm["X_4ICs"]

    Xs = Xs[:, components]

    _adata_src.obsm["X_4ICs"] = _raw_ot_integration(
        Xs, Xt, metric=metric, eps=eps, verbose=verbose
    )
    _adata_ref.obsm["X_4ICs"] = Xt

    if verbose:
        print("-- Done")

    # -- Improve X_dimRed3d so that all datasets are expressed in the same {PC} subspace
    pca = PCA(n_components=3, svd_solver="arpack")
    _adata_src.obsm["X_dimRed3d"] = pca.fit_transform(_adata_src.obsm["X_4ICs"])
    
    #-- Add integration arguments
    ref_dimred = _adata_ref.uns['scycle']['dimRed']
    ref_dimred['run_on_reference'] = True
    
    _adata_src.uns['scycle']['dimRed'] = ref_dimred
    _adata_src.uns['P_dimRed'] = _adata_ref.uns['P_dimRed']
    
    _adata_src.uns['scycle']['integration'] = dict(
            components = components,
            metric = metric,
            eps=eps
            )


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
    if verbose:
        print("> Computing distance matrix...")
    M = dist(Xs, Xt, metric=metric)
    M /= M.max()
    if verbose:
        print("> Computing optimal transport plan...")
    Gs = sinkhorn(w_x, w_y, M, eps)
    if verbose:
        print("> Projecting source dataset...")
    return np.diag(1 / w_x) @ Gs @ Xt
