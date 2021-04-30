#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    geom_point,
    geom_line,
    geom_segment,
    geom_vline,
    geom_text,
    labs,
)
from plotnine.scales import scale_color_manual
from ._pseudotime_scatter import pseudotime_scatter
from scipy.stats import zscore
import warnings

pd.set_option("chained_assignment", None)


def cell_cycle_scores(
    adata,
    scores=["signatures", "components"][0],
    size=1.5,
    alpha=1,
    curvature_shrink=1,
    lab_ypos=2,
):
    """Plots cell cycle signatures vs pseudotime

    Parameters
    ----------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.cell_cycle_phase`.
    scores: str
        A string indicating what to plot as cell cycle scores against pseudotime.
        If 'signatures', standard S-phase, G2-M and Histones signatures are used;
        if 'components', the 4 cell cycle related components are used.
    size: float
        Controls the point size of the plot.
    alpha: float
        A value between 0 and 1. Controls point transparency.
    lab_ypos: float
        Controls the y-axis position of the cell cycle phase annotation.

    Returns
    --------------
    A plotnine scatter plot of pseudotime vs 3 cell cycle signatures.

    """
    if scores == "signatures":
        y = ["G1-S", "G2-M", "Histones"]
        colors = ['#8ca0c9', '#ff8d68', '#5cc2a6', "black"]
    elif scores == "components":
        _add_compScores(adata)
        y = ["G1-S comp", "G2-M comp", "G2-M- comp", "Histone comp"]
        colors = ['#8ca0c9', '#ff8d68', "#e5c494", '#5cc2a6', "black"]

    time_scatter = pseudotime_scatter(adata, y=y, size=size, alpha=alpha) + labs(
        x="Pseudotime", y="Signature scores", color="Signature"
    )

    # -- Add cell cycle annotations
    if "cell_cycle_division" in adata.uns["scycle"]:
        cc_divs = adata.uns["scycle"]["cell_cycle_division"]

        # -- Curvature data
        curv_data = cc_divs["curvature"]
        curv = curv_data["curvature"].values
        cvz = zscore(curv) / curvature_shrink
        cvz = cvz - np.max(cvz)
        curv_data.loc[:, "curvature"] = cvz
        curv_data.loc[:, "signature"] = "Curvature"

        # -- Peak data (for segments)
        gr_min = np.min(curv_data["curvature"])
        pk_data = curv_data[curv_data["ispeak"] == "peak"]
        pk_data.loc[:, "ymin"] = gr_min

        # -- Cell cycle annotation
        cc_phase = pd.DataFrame(
            dict(
                starts=[
                    None,
                    cc_divs["s_start"],
                    cc_divs["g2_start"],
                    cc_divs["m_start"],
                ],
                labels=["G1", "S", "G2", "M"],
                labpos=[
                    np.mean([0, cc_divs["s_start"]]),
                    np.mean([cc_divs["s_start"], cc_divs["g2_start"]]),
                    np.mean([cc_divs["g2_start"], cc_divs["m_start"]]),
                    np.mean([cc_divs["m_start"], 1]),
                ],
                y=lab_ypos,
            )
        )

        cell_cycle_plt = (
            time_scatter
            + geom_point(
                aes("pseudotime", "curvature", color="signature"), data=curv_data
            )
            + geom_line(aes("pseudotime", "curvature"), data=curv_data)
            + scale_color_manual(values=colors)
            + geom_segment(
                aes(x="pseudotime", xend="pseudotime", y="ymin", yend="curvature"),
                linetype="dotted",
                data=pk_data,
            )
            + geom_vline(aes(xintercept="starts"), linetype="dashed", data=cc_phase)
            + geom_text(aes(x="labpos", y="y", label="labels"), data=cc_phase)
        )

        return cell_cycle_plt
    else:
        return time_scatter + scale_color_manual(values=colors[0:-1])


def _add_compScores(adata):
    ecom_idx = adata.uns['scycle']['find_cc_components']['indices']
    for comp, i in ecom_idx.items():
        cname = comp + ' comp'
        adata.obs[cname] = zscore(adata.obsm['X_dimRed'][:,i])

def scatter_cell_cycle(
    adata,
    scores=["signatures", "components"][0],
    size=1.5,
    alpha=1,
    curvature_shrink=1,
    lab_ypos=2,):
    """DEPRECATED: cell_cycle_scores
    """
    warnings.warn("scatter_cell_cycle is deprecated; use cell_cycle_scores", DeprecationWarning)
    return (cell_cycle_scores(adata, scores, size, alpha, curvature_shrink, lab_ypos))


