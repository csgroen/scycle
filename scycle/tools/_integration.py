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
        entropy: bool = False,
        hreg: float = 1e-3,
        unbalanced: bool = False,
        mreg: float = 1e-2,
        weighting_strategy: str = "uniform",
        downsampling: bool = False,
        normalize: bool = True,
        jitter_std: float = .02,
        verbose: bool = False,
        max_iter: int = 1e7,
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
    components: list<int>
        ICs indices caring cell-cycle related information. Let
        it empty for automatic detection.
    method: str
        In "ot", "gromov". Histogram distance to use.
    entropy: bool
        Use the entropy regularized solver, turn it on for larger problems
    hreg: float
        Entropy regularization constant, to tune when entropy=True
    unbalanced: bool
        Use unbalanced optimal transport formulation. Helps for problems
        with unbalanced cell types between datasets.
    mreg: float
        Marginal penalty, to tune if unbalanced=True
    weighting_strategy: str
        in "uniform" (default), "woti". How to choose weights before optimal
        transport. In case of severely unbalanced dataset in terms of cell
        phases, try "woti".
    downsampling: bool
        Turn it on for large datasets, at a cost of an approximate solution.
    normalize: bool
        Use column-normalized datasets for the integration. Helps when there
        are strong differences between features.
    jitter_std: float
        Quantity of jittering to apply after integration.
    verbose: bool
        Outputs information in standard output stream.
    max_iter: int
        Maximum number of iterations for the optimal transport plan computation.
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
        entropy=entropy,
        hreg=hreg,
        unbalanced=unbalanced,
        mreg=mreg,
        normalize=normalize,
        metric="sqeuclidean",
        n_hops=downsampling,  # slight abuse bool -> int
        weighting_strategy=weighting_strategy,
        max_iter=max_iter,
        verbose=(2 if verbose else 0)
    )
    _adata_src.obsm["X_cc"] = w.fit_transform(Xs, Xt, jitter_std=jitter_std)
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
