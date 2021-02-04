#!/usr/bin/env python3

import numpy as np
from ._find_cc_components import find_cc_components
from ._subtract_cc import subtract_cc


def classify_genes_by_ic(adata, min_r2=0.5, verbose=False):

    if "r2_scores" not in adata.var:
        print(
            "Cell cycle signal was not subtracted. Doing it now with default values..."
        )
        subtract_cc(adata)

    if "find_cc_components" not in adata.uns["scycle"]:
        print(
            "Cell cycle components were not found. Doing it now with default values..."
        )
        find_cc_components(adata, verbose=verbose)

    # Computing most important IC for each gene
    components = list(adata.uns["scycle"]["find_cc_components"].values())
    is_4d = len(components) == 4
    maxes = np.argmax(adata.varm["P_dimRed"][components, :], axis=0)

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
    ] = (3 if is_4d else 2)

    # trandforming ids -> str
    maxes = np.where(maxes == 0, "G1/S", maxes)
    maxes = np.where(maxes == "1", "G2/M+", maxes)
    if is_4d:
        maxes = np.where(maxes == "2", "G2/M-", maxes)
    maxes = np.where(maxes == ("3" if is_4d else "2"), "HIST", maxes)
    maxes[adata.var["r2_scores"] < min_r2] = "unrelated"

    # writing in adata
    adata.var["cc_class"] = maxes
