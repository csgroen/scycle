#!/usr/bin/env python3

import numpy as np
from ._enrich_components import enrich_components
from ._subtract_cc import subtract_cc


def classify_genes_by_ic(adata, min_r2=0.5, verbose=False):

    if "r2_scores" not in adata.var:
        print(
            "Cell cycle signal was not subtracted. Doing it now with default values..."
        )
        subtract_cc(adata)

    if "enrich_components" not in adata.uns["scycle"]:
        print(
            "Components enrichment was not performed. Doing it now with default values..."
        )
        enrich_components(adata, verbose=verbose)

    # Getting metagenes matrix
    X_4ICs = adata.uns["P_dimRed"][
        list(adata.uns["scycle"]["enrich_components"].values())
    ].T

    # Computing most important IC for each gene
    maxes = np.argmax(X_4ICs, axis=1)

    # Correcting for histone genes
    var_startswith = adata.var_names.str.startswith
    maxes[
        (
            var_startswith("HIST")
            + var_startswith("H1")
            + var_startswith("H2")
            + var_startswith("H3")
            + var_startswith("H4")
        ).astype(bool)
    ] = 3

    # trandforming ids -> str
    maxes = np.where(maxes == 0, "G1/S", maxes)
    maxes = np.where(maxes == "1", "G2/M-", maxes)
    maxes = np.where(maxes == "2", "G2/M+", maxes)
    maxes = np.where(maxes == "3", "HIST", maxes)
    maxes[adata.var["r2_scores"] < min_r2] = "unrelated"

    # writing in adata
    adata.var["cc_class"] = maxes
