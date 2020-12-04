#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from typing import Optional

# from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from algorithms.stabilized_ICA import StabilizedICA
from anndata import AnnData
from ..annot import cellcycle_signatures


def dimensionality_reduction(
    adata: AnnData,
    method: str = "pcaCCgenes",
    n_comps: int = 30,
    sig_names: list = ["S-phase", "G2-M"],
    seed: Optional[int] = None,
    max_iter: int = 200,
    verbose: bool = True,
):
    """Dimensionality reduction for pseudotime computation

    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        pp.prep_simple or pp.prep_pooling.
    method: str
        Method of dimensionality reduction, currently one of: 'pca', 'ica',
        'pcaCCgenes','icaCCgenes','CCgenes'.
        The 'CCgenes' variant use the methods only on the set of cell-cycle genes.
        'CCgenes' uses the G1-S and G2-M signature scores as the reduced dimensions.
    n_comps: int
        Number of components to use for the dimensionality reduction.
    sig_names: list
        Used only for 'CCgenes'. List of signature names to use as reference
        cell cycle dimensions. Must be present as scores per sample in adata.obs
    seed: int
        Sets up the random state for FastICA and NMF for results reproducibility.
        Used for 'nmf', 'ica', 'nmfCCgenes' and 'icaCCgenes' methods.
    max_iter: int
        Maximum number of iterations during FastICA and NMF fit.
    verbose: bool
        If True, messages about function progress will be printed.

    Returns
    ------------------
    `adata` will be updated with the dimensionality reduction results to be
    used for the next steps in the pipeline.
    """

    if "scycle" not in adata.uns:
        raise Exception(
            (
                "Data needs to be pre-processed by `pp.prep_pooling`"
                + "or `pp.prep_simple`, before dimensionality reduction"
            )
        ) 
        
    if method in ['pcaCCgenes', 'icaCCgenes']:
        adata_cc = _adata_CCgenes(adata, sig_names)
        genes = np.array(adata_cc.var.index)
    else:
        genes = np.array(adata.var.index)

    if method == "pca":
        dimred_res = _dimRed_pca(adata, n_comps=n_comps, verbose=verbose)
        # elif method == "pcaCCgenes":
        #    dimred_res = _dimRed_pca(adata_cc, n_comps=n_comps, verbose=verbose)
    elif method == "ica":
        dimred_res = _dimRed_ica(
            adata, n_comps=n_comps, max_iter=max_iter, seed=seed, verbose=verbose
        )
    elif method == 'pcaCCgenes':
        dimred_res = _dimRed_pca(
            adata_cc, n_comps=n_comps, verbose=verbose
        )
    elif method == "icaCCgenes":
        dimred_res = _dimRed_ica(adata_cc, max_iter=max_iter, seed=seed,  n_comps=n_comps, verbose=verbose)

        # elif method == 'nmf': dimred_res = _dimRed_nmf(adata, n_comps = n_comps, max_iter = max_iter, seed = seed, verbose = verbose)
        # elif method == 'nmfCCgenes': dimred_res = _dimRed_nmf(adata_cc, n_comps = n_comps, max_iter = max_iter, seed = seed, verbose = verbose)
        # elif method == "CCgenes":
        #    dimred_res = _dimRed_CCgenes(adata, verbose=verbose)
    else:
        raise Exception(
            (
                "Not one of the supported methods.\n"
                + "Must be one of: pca, ica, nmf, pcaCCgenes, icaCCgenes"
            )
        )

    X_dimRed = dimred_res["dimred"]
    adata.obsm["X_dimRed"] = X_dimRed
    adata.uns["dimRed"] = dimred_res["obj"]
    adata.uns["P_dimRed"] = dimred_res["pMatrix"]

    # -- 3D
    pca_dimRed = PCA(n_components=3)
    pca_dimRed.fit(X_dimRed)
    
    adata.obsm["X_dimRed3d"] = pca_dimRed.transform(X_dimRed)

    adata.uns["scycle"]["dimRed"] = {
        "method": method,
        "features": genes,
        "n_comps": n_comps,
        "seed": seed,
    }


def _adata_CCgenes(adata, sig_names):
    # -- Get CC genes
    cc_sigs = cellcycle_signatures()
    cc_sigs = {sign: cc_sigs[sign] for sign in sig_names}

    flat_sigs = [
        item for sublist in [v for k, v in cc_sigs.items()] for item in sublist
    ]
    cc_genes = np.unique(np.array(flat_sigs))
    # -- Select in matrix
    adata_cc = adata.copy()
    idx = [gene in adata_cc.var_names.tolist() for gene in cc_genes]
    adata_cc = adata_cc[:, cc_genes[idx]]
    return adata_cc


def _dimRed_pca(adata, n_comps, verbose=False):
    if verbose:
        print("Dimensionality reduction using PCA...")
    pca = PCA(n_components=n_comps)
    pca.fit(adata.X)
    X_dimRed = pca.transform(adata.X)
    return {"obj": pca, "dimred": X_dimRed, "pMatrix": pca.components_}


def _dimRed_ica(adata, n_comps, max_iter, seed, verbose=False, max_trials=10):
    # -- Run ICA
    if verbose:
        print("-- Dimensionality reduction using ICA...")
    X = adata.X.copy()
    X -= np.mean(X, axis=0)
    sICA = StabilizedICA(n_components=n_comps, max_iter=2000, n_jobs=-1)
    for i in range(max_trials):
        try:
            sICA.fit(X.T, n_runs=100, plot=False, normalize=True, fun="logcosh")
        except ValueError:
            if verbose:
                print("*- ICA did not converge. Retrying...")
            continue
        break
    if verbose:
        print("-- Done")
    return {"obj": sICA, "dimred": X @ sICA.S_.T, "pMatrix": sICA.S_}


# def _dimRed_nmf(adata, n_comps, max_iter, seed, verbose = False):
#     if verbose: print('-- Dimensionality reduction using NMF...')
#     nmf = NMF(n_components = n_comps, max_iter = max_iter, random_state=seed)
#     nmf.fit(adata.X)
#     X_dimRed = nmf.transform(adata.X)
#     return {'obj': nmf, 'dimred': X_dimRed}


# def _dimRed_CCgenes(adata, verbose=False):
#     if verbose:
#         print("-- Dimensionality reduction using G1 and G2/M signatures...")
#     X_dimRed = adata.obs.loc[:, ["G1", "G2-M"]]
#     return {"obj": "CC_genes", "dimred": X_dimRed}


# def _dimRed_CCA(adata, n_comps, cl_var, verbose = False):
#     if cl_var == None:
#         raise(Exception("cl_var can't be None if using Canonical Correlation Analysis"))
#     if verbose: print('-- Dimensionality reduction using CCA...')
#     #-- Make adata_list
#     cl = adata.obs[cl_var].values
#     idxs = [cl == cat for cat in cl.categories]
#     adata_list = [adata.X[idx,:] for idx in idxs]
#     n = np.min([mat.shape[0] for mat in adata_list])
#     adata_list_min = [mat[0:n,:] for mat in adata_list]
#     cca = rcca.CCA(numCC = n_comps, kernelcca = False)
#     cca.train(adata_list_min)


def _dimRed_decomp(adata, decomp, feats, ccomps):
    return {"obj": decomp, "dimred": decomp.transform(adata.X[:, feats])[:, ccomps]}
