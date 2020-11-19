#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import Optional

from sklearn.decomposition import PCA, FastICA
import numpy as np
from anndata import AnnData
from ..annot import cellcycle_signatures


def dimensionality_reduction(
    adata: AnnData,
    method: str = "pca",
    n_comps: int = 30,
    sig_names: list = ["G1", "S-phase", "Histones", "G2-M"],
    seed: Optional[int] = None,
    max_iter: int = 200,
    verbose: bool = True,
    pp_by_scycle: bool = True,
):
    """Dimensionality reduction for pseudotime computation

    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        pp.prep_simple or pp.prep_pooling.
    method: str
        Method of dimensionality reduction, currently one of: 'pca', 'ica',
        'pcaCCgenes' or 'icaCCgenes'. The 'CCgenes' variant use the methods
        only on the set of cell-cycle genes.
    n_comps: int
        Number of components to use for the dimensionality reduction.
    sig_names: list
        Used only for 'pcaCCgenes' or 'icaCCgenes'. List of signature names
        to reference as 'cell-cycle genes'.
    seed: int
        Sets up the random state for FastICA for results reproducibility.
        Used only for 'ica' or 'icaCCgenes' methods.
    max_iter: int
        Maximum number of iterations during FastICA fit.
    verbose: bool
        If True, messages about function progress will be printed.
    pp_by_scycle: bool
        Set it to true if adata was preprocessed outside of scycle

    Returns
    ------------------
    `adata` will be updated with the dimensionality reduction results to be
    used for the next steps in the pipeline.
    """

    if pp_by_scycle:
        if "scycle" not in adata.uns:
            raise Exception(
                (
                    "Data needs to be pre-processed by `pp.prep_pooling`"
                    + "or `pp.prep_simple`, before dimensionality reduction"
                )
            )
    else:
        adata.uns["scycle"] = {}

    if method == "pca":
        if verbose:
            print("Dimensionality reduction using PCA...")
        pca = PCA(n_components=n_comps)
        pca.fit(adata.X)
        X_dimRed = pca.transform(adata.X)
        adata.obsm["X_dimRed"] = X_dimRed
        adata.uns["dimRed"] = pca

    elif method == "ica":
        # -- Run ICA
        if verbose:
            print("-- Dimensionality reduction using ICA...")
        # ica = StabilizedICA(n_components = n_comps, max_iter = max_iter)
        # ica_fit = ica.fit(adata.X, n_runs = n_runs)
        ica = FastICA(n_components=n_comps, max_iter=max_iter, random_state=seed)
        ica.fit(adata.X)
        X_dimRed = ica.transform(adata.X)

        adata.obsm["X_dimRed"] = X_dimRed
        adata.uns["dimRed"] = ica

    elif method in ["pcaCCgenes", "icaCCgenes"]:
        # -- Get CC genes
        cc_sigs = cellcycle_signatures()
        flat_sigs = [
            item for sublist in [v for k, v in cc_sigs.items()] for item in sublist
        ]
        cc_genes = np.unique(np.array(flat_sigs))
        # -- Select in matrix
        adata_cc = adata.copy()
        idx = [gene in adata_cc.var_names.tolist() for gene in cc_genes]
        adata_cc = adata_cc[:, cc_genes[idx]]

        if method == "icaCCgenes":
            ica = FastICA(n_components=n_comps, random_state=seed)
            ica.fit(adata_cc.X)
            X_dimRed = ica.transform(adata_cc.X)

            adata.obsm["X_dimRed"] = X_dimRed
            adata.uns["dimRed"] = ica

        else:
            pca = PCA(n_components=n_comps)
            pca.fit(adata_cc.X)
            X_dimRed = pca.transform(adata_cc.X)

            adata.obsm["X_dimRed"] = X_dimRed
            adata.uns["dimRed"] = pca

    else:
        raise Exception(
            (
                "Not one of the supported methods.\n"
                + "Must be one of: pca, ica, pcaCCgenes, icaCCgenes"
            )
        )

    pca_dimRed = PCA(n_components=2)
    pca_dimRed.fit(X_dimRed)

    adata.obsm["X_dimRed2d"] = pca_dimRed.transform(X_dimRed)
    adata.uns["scycle"]["dimRed"] = {"method": method, "n_comps": n_comps, "seed": seed}
