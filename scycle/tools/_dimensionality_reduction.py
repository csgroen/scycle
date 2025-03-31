#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import elpigraph
import warnings
from typing import Optional

# from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from sica.base import StabilizedICA
from anndata import AnnData
from ..data import cc_genes
from ._find_cc_components import find_cc_components
from ._subtract_cc import _cc_residue_matrix

import gc


def dimensionality_reduction(
    adata: AnnData,
    method: str = "pcaCCgenes",
    n_comps: int = 30,
    seed: Optional[int] = None,
    max_iter: int = 200,
    self_consistent_kwargs={
        "n_nodes": 30,
        "max_iter": 10,
        "eps_jaccard": 0.02,
        "Mu": 0.01,
        "r2_threshold": 0.5,
    },
    find_cc_comps: bool = True,
    find_cc_comp_thr: int = 3,
    verbose: bool = True,
):
    """Dimensionality reduction for pseudotime computation

    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        pp.prep_simple or pp.prep_pooling.
    method: str
        Method of dimensionality reduction, currently one of: 'cc_signatures', 'pca', 'ica',
        'pcaCCgenes','self_consistent_CC'.
        The 'self_consistent_pcaCC' find the self-consistent cell cycle space based on an
        initial estimation using ICA and finding the cell cycle ICs, and
        then using the trajectory to find the genes that are most informative
        to find cell cycle space.
    n_comps: int
        Number of components to use for the dimensionality reduction.
    seed: int
        Sets up the random state for FastICA and NMF for results reproducibility.
        Used for 'ica' and'icaCCgenes' methods.
    max_iter: int
        Maximum numhistone_markersber of iterations during FastICA and NMF fit.
    self_consistent_kwargs:
        A dictionary of arguments used if method = 'self_consistent_cc'.
    find_cc_comps: bool
        If True and method='ica', components will be scored to find
        cell cycle-related components
    find_cc_comp_thr: int
      Score threshold for IC recognition. Passed to tl.find_cc_components
    verbose: bool
        If True, messages about function progress will be printed.

    Returns
    ------------------
    `adata` will be updated with the dimensionality reduction results to be
    used for the next steps in the pipeline.
    """

    if "scycle" not in adata.uns:
        adata.uns['scycle'] = dict()
        warnings.warn(
            "Data has not been pre-processed using scycle (`pp.prep_pooling`" +
            "or `pp.prep_simple`), before dimensionality reduction")
    if method in ["pcaCCgenes", "icaCCgenes"]:
        adata_cc = _adata_CCgenes(adata)
        genes = np.array(adata_cc.var.index)
    else:
        genes = np.array(adata.var.index)

    if method == "pca":
        dimred_res = _dimRed_pca(adata, n_comps=n_comps, verbose=verbose)
    elif method == "cc_signatures":
        dimred_res = _dimRed_ccsigs(adata, verbose=verbose)
    elif method == "ica":
        dimred_res = _dimRed_ica(
            adata, n_comps=n_comps, max_iter=max_iter, seed=seed, verbose=verbose
        )
    elif method == "pcaCCgenes":
        dimred_res = _dimRed_pca(adata_cc, n_comps=n_comps, verbose=verbose)
    # elif method == "icaCCgenes":
    #     dimred_res = _dimRed_ica(
    #         adata_cc, max_iter=max_iter, seed=seed, n_comps=n_comps, verbose=verbose
    #     )
    elif method == "self_consistent_CC":
        dimred_res = _dimRed_self_consistent_cc(
            adata,
            n_comps=n_comps,
            max_iter=self_consistent_kwargs["max_iter"],
            n_nodes=self_consistent_kwargs["n_nodes"],
            eps_jaccard=self_consistent_kwargs["eps_jaccard"],
            r2_threshold=self_consistent_kwargs["r2_threshold"],
            verbose=verbose,
        )

    else:
        raise Exception(
            (
                "Not one of the supported methods.\n"
                + "Must be one of: pca, ica, cc_signatures, pcaCCgenes, icaCCgenes or self_consistent_CC"
            )
        )

    X_dimRed = dimred_res["dimred"]
    adata.obsm["X_dimRed"] = X_dimRed
    adata.uns["dimRed"] = dimred_res["obj"]

    # -- 3D
    if method == 'cc_signatures':
        adata.obsm["X_pca_scycle"] = np.array(adata.obs[['G1-S', 'G2-M']])
        adata.uns["scycle"]["dimRed"] = {"method": method}
    else:
        if method in ["ica", "pca"]:
            adata.varm["P_dimRed"] = dimred_res["pMatrix"]
        else:
            adata.uns["P_dimRed"] = dimred_res["pMatrix"]

        pca_dimRed = PCA(n_components=3)
        pca_dimRed.fit(X_dimRed)
        adata.obsm["X_pca_scycle"] = pca_dimRed.transform(X_dimRed)

        adata.uns["scycle"]["dimRed"] = {
            "method": method,
            "features": genes,
            "n_comps": n_comps,
            "seed": seed,
        }
    if method == "self_consistent_CC":
        adata.uns["scycle"].pop("find_cc_components")
        adata.varm.pop('P_dimRed')
        del adata.obsm["X_cc"]
        adata.uns["scycle"]["dimRed"]["cc_genes"] = dimred_res["cc_genes"]

    if method == "ica" and find_cc_comps:
        find_cc_components(adata, thr = find_cc_comp_thr, verbose=verbose)

    gc.collect()


def _adata_CCgenes(adata):
    # -- Select in matrix
    adata_cc = adata.copy()
    idx = [gene in adata_cc.var_names.tolist() for gene in cc_genes]
    genes2keep = np.array(cc_genes)[idx]
    adata_cc = adata_cc[:, genes2keep]
    return adata_cc


def _dimRed_pca(adata, n_comps, verbose=False):
    if verbose:
        print("Dimensionality reduction using PCA...")
    pca = PCA(n_components=n_comps)
    pca.fit(adata.X)
    X_dimRed = pca.transform(adata.X)
    return {"obj": pca, "dimred": X_dimRed, "pMatrix": np.linalg.pinv(pca.components_)}


def _dimRed_ica(adata, n_comps, max_iter, seed, verbose=False):
    # -- Run ICA
    if verbose:
        print("-- Dimensionality reduction using ICA...")
    X = adata.X.copy()
    np.subtract(X, np.mean(X, axis=0), out=X)
    sICA = StabilizedICA(n_components=n_comps, max_iter=2000, n_jobs=-1, n_runs=100, plot=False, normalize=True, fun="logcosh")
    sICA.fit(X)
    if verbose:
        print("-- Done")
    S_pi = np.linalg.pinv(sICA.S_)
    dimred = X @ S_pi
    del X
    return {"obj": sICA, "dimred": dimred, "pMatrix": S_pi}


def _dimRed_ccsigs(adata, verbose=False):
    if verbose:
        print("-- Dimensionality reduction using G1-S and G2-M signatures...")
    dimred = np.array(adata.obs[['G1-S', 'G2-M']])
    return {"obj": None, "dimred": dimred}


def _dimRed_decomp(adata, decomp, feats, ccomps):
    return {"obj": decomp, "dimred": decomp.transform(adata.X[:, feats])[:, ccomps]}


def _dimRed_self_consistent_cc(
    adata,
    n_comps=30,
    n_nodes=30,
    max_iter=10,
    Mu=0.01,
    r2_threshold=0.5,
    eps_jaccard=0.02,
    verbose=True,
):

    if verbose:
        print("Dimensionality reduction using self consistent CC genes..")
        print("Initial estimation using ICA...")

    dimensionality_reduction(adata, method="ica", find_cc_comps=True, verbose=False)
    Xr = adata.obsm["X_cc"]

    # -- Run initial iteration with cc_genes
    r2scores, ind = _genes_ccspace(adata.X, Xr, n_comps, n_nodes, Mu, r2_threshold)
    top_thr = np.min([np.quantile(r2scores[ind], 0.9)])

    if verbose:
        print("Cell cycle genes initially found:", len(ind))
        print("Top ones:", list(adata.var_names[r2scores > top_thr]))

    ind_old = ind
    for counter in range(max_iter):
        X_cc = adata.X[:, ind]
        Xr = PCA(n_components=np.min([n_comps, X_cc.shape[1]])).fit_transform(X_cc)
        r2scores, ind = _genes_ccspace(adata.X, Xr, n_comps, n_nodes, Mu, r2_threshold)

        idx_common_genes = list(set(ind) & set(ind_old))
        union_genes = list(set(ind) | set(ind_old))
        perc = len(idx_common_genes) / len(union_genes)

        if verbose:
            print(
                "\n\nIteration",
                counter + 1,
                "==================\nJaccard coeff:",
                perc,
                "Old:",
                len(ind_old),
                "New:",
                len(ind),
                "\n==================\n",
            )
        if perc > 1 - eps_jaccard:
            break
        ind_old = ind.copy()

    print("Reducing dimensionality using self-consistent CC genes...")
    genes_cc_space = adata.var_names[ind]
    pca = PCA(n_components=n_comps).fit(adata.X[:, ind])
    X_dimRed = pca.transform(adata.X[:, ind])

    return {
        "obj": pca,
        "dimred": X_dimRed,
        "pMatrix": np.linalg.pinv(pca.components_),
        "cc_genes": genes_cc_space.tolist(),
    }


def _genes_ccspace(X, Xr, n_comps, n_nodes, Mu, r2_threshold):
    # -- Get initial trajectory
    egr = elpigraph.computeElasticPrincipalCircle(Xr, n_nodes, Mu=Mu, verbose=False)
    nodep = egr[0]["NodePositions"]
    partition, dists = elpigraph.src.core.PartitionData(
        X=Xr,
        NodePositions=nodep,
        MaxBlockSize=100000000,
        TrimmingRadius=np.inf,
        SquaredX=np.sum(Xr ** 2, axis=1, keepdims=1),
    )

    residue_matrix, r2scores = _cc_residue_matrix(X, partition.flatten())
    ind = np.where(r2scores > r2_threshold)[0]
    return (r2scores, ind)
