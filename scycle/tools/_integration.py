#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import transmorph as tr

from anndata import AnnData
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA

from ._find_cc_components import find_cc_components


def integration(
    _adata_src: AnnData,
    _adata_ref: AnnData,
    components: list = [],
    method: str = "ot",
    entropy: bool = True,
    hreg: float = 1e-3,
    weighted: bool = False,
    verbose: bool = True,
    max_iter: int = 1e7,
    scale: float = 0.1,
    alpha_qp: float = 1.0,
):
    """
    Optimal transport-based data integration. Data dimensionality must have been reduced
    before applying the method. Both source and reference must share the same cloud topology,
    otherwise the technique is very prone to overfitting.

    Parameters
    ----------
    _adata_src: AnnData
        The dataset to align.
    _adata_ref: AnnData
        The dataset to use as a reference.
    components: list
        ICs indices caring cell-cycle related information. Let
        it empty for automatic detection.
    method: str
        In "ot", "gromov". Histogram distance to use.
    entropy: bool
        Use the entropy regularized solver, turn on for larger problems
    hreg: float
        Entropy regularization constant
    weighted:
        False to use uniform weights, True to estimate them.
    verbose: bool
        Outputs information in standard output stream.
    max_iter: int
        Maximum number of iterations for the optimal transport plan computation.
    scale_src: float
        For transmorph only.
        Kernel scaling of the source cloud point.
    scale_ref: float
        For Transmorph only.

    alpha_kde: float
        For Transmorph only.
        Alpha parameter for KDE bandwith selection, between 0 and 1.
    alpha_qp:
        For Transmorph only.
        Alpha parameter for quadratic program solver (OSQP), between 0 and 2
    """
    assert method in ("ot", "gromov"), (
        "Method %s not recognized. Options: 'ot', 'gromov'" % method
    )
    assert (
        "dimRed" in _adata_ref.uns
    ), "ICA projection matrix missing in source AnnData. \
        Dimensionality reduction must be computed first."

    # Selecting the cell-cycle ICs components
    if len(components) == 0:
        if verbose:
            print("-- Automatically detecting cell-cycle components...")
        if "find_cc_components" not in _adata_ref.uns["scycle"].keys():
            find_cc_components(_adata_ref)
        components = list(_adata_ref.uns["scycle"]["find_cc_components"]["indices"].values())

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
            if idx >= len(common_genes):
                break
    adata_ref = _adata_ref[:, common_genes]
    adata_src = _adata_src[:, common_genes]

    # Projecting both datasets in the same ICs subspace
    if "X_dimRed" not in adata_ref.obsm:
        if verbose:
            print("> Projecting reference dataset...")
        Y_c = adata_ref.X.copy()
        Y_c -= np.mean(Y_c, axis=0)
        adata_ref.obsm["X_dimRed"] = Y_c @ adata_ref.varm["P_dimRed"]

    X_c = adata_src.X.copy()
    X_c -= np.mean(X_c, axis=0)
    adata_src.obsm["X_dimRed"] = X_c @ adata_ref.varm["P_dimRed"]

    Xs: np.ndarray = adata_src.obsm["X_dimRed"]
    Xt: np.ndarray = adata_ref.obsm["X_cc"]

    Xs = Xs[:, components]

    # -- Update objects
    if verbose:
        print("> Performing optimal transport based integration using Transmorph...")
    w = tr.Transmorph(
        method=method,
        max_iter=max_iter,
        entropy=entropy,
        hreg=hreg,
        weighted=weighted,
        alpha_qp=alpha_qp,
        scale=scale,
        verbose=verbose
    )
    _adata_src.obsm["X_cc"] = w.fit_transform(Xs, Xt)
    _adata_ref.obsm["X_cc"] = Xt
    _adata_src.obsm["X_dimRed"] = adata_src.obsm["X_dimRed"]

    if verbose:
        print("-- Done")

    # -- Update X_pca_scycle so that all datasets are expressed in the same {PC} subspace
    nc = 3 if _adata_ref.obsm['X_cc'].shape[1] > 2 else 2
    pca = PCA(n_components=nc)
    pca.fit(_adata_ref.obsm['X_cc'])
    _adata_src.obsm["X_pca_scycle"] = pca.transform(_adata_src.obsm['X_cc'])

    # -- Add integration arguments
    ref_dimred = _adata_ref.uns["scycle"]["dimRed"]
    ref_dimred["run_on_reference"] = True

    _adata_src.uns["scycle"]["dimRed"] = ref_dimred
    _adata_src.varm["P_dimRed"] = _adata_ref.varm["P_dimRed"]
