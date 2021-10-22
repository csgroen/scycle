#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from plotnine import ggplot, aes, geom_point, labs, facet_wrap, geom_vline, geom_text
from ._themes import theme_std
import warnings

def pseudotime_scatter(adata, y, facet = True, size = 1.5, alpha = 1,
                       color = 'black', ncol = 2, lab_ypos = 2):
    """Plots a scatter plot of pseudotime vs one or multiple variables

    Parameters
    --------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.pseudotime`.
    y: str or list
        If type(y) == str, y must be a variable annotated in adata.obs and
        will be used as the y-axis. If type(y) == list, then multiple variables
        will be plotted using a shared y-axis but different point colors.
    facet: bool
        Whether to return a facetted plot or all signatures in a single plot.
        Only used if y is a list.
    size: float
        Controls the point size of the plot.
    alpha: float
        A value between 0 and 1. Controls point transparency.
    color: str
        A supported color name. Controls the point color if type(y)==str.
        Ignored otherwise.
    ncol: int
        Number of columns in the facetting, if facet=True. Ignored otherwise.
    lab_ypos: float
        Controls the y-axis position of the cell cycle phase annotation, if present.

    Returns
    -------------
    A plotnine scatter plot of pseudotime.
    """
    if type(y) == str:
        #-- Get data
        if y in adata.obs.columns:
            plot_df = pd.DataFrame({'x': adata.obs['pseudotime'], 'y': adata.obs[y]})
        elif y in adata.var_names:
            plot_df = pd.DataFrame({'x': adata.obs['pseudotime'], 'y': adata[:,y].X.flatten()})

        #-- Make plot
        if color in adata.obs.columns:
            time_scatter = (ggplot(plot_df, aes(x = 'x', y = 'y'))
              + geom_point(aes(color = color), size = size, alpha = alpha)
              + labs(x = 'Pseudotime', y = y)
              + theme_std)
        else:
            time_scatter = (ggplot(plot_df, aes(x = 'x', y = 'y'))
              + geom_point(size = size, alpha = alpha, color = color)
              + labs(x = 'Pseudotime', y = y)
              + theme_std)

    else:
        #-- Make multiple color plot
        sannot = pd.DataFrame({'pseudotime': adata.obs['pseudotime']})
        sannot['id'] = range(sannot.shape[0])

        #-- Get y from obs or matrix:
        for var in y:
            if var in adata.obs.columns:
                sannot[var] = adata.obs[var]
            elif var in adata.var_names:
                sannot[var] = adata[:,var].X.flatten()
        plot_df = pd.melt(sannot, id_vars = ['id', 'pseudotime'],
                          var_name = 'signature', value_name = 'score')
        plot_df['signature'] = plot_df['signature'].astype('category')
        plot_df['signature'].cat.reorder_categories(y, inplace=True)

        if facet:
            time_scatter = (ggplot(plot_df, aes('pseudotime', 'score'))
             + facet_wrap('signature', scales = 'free_y', ncol = ncol)
             + geom_point(aes(color = 'signature'), alpha = alpha, size = size)
             + theme_std)
        else:
            time_scatter = (ggplot(plot_df, aes('pseudotime', 'score'))
             + geom_point(aes(color = 'signature'), alpha = alpha, size = size)
             + theme_std)

    if "cell_cycle_division" in adata.uns["scycle"]:
        cc_divs = adata.uns["scycle"]["cell_cycle_division"]
        # -- Cell cycle annotation
        cc_phase = pd.DataFrame(
            dict(
                starts=[
                    None,
                    cc_divs["pr_start"],
                    cc_divs["rep_start"],
                    # cc_divs["m_start"],
                ],
                labels=["G1 PM", "G1 PR", "S/G2/M"],
                labpos=[
                    np.mean([0, cc_divs["pr_start"]]),
                    np.mean([cc_divs["pr_start"], cc_divs["rep_start"]]),
                    np.mean([cc_divs["rep_start"], 1]),
                    # np.mean([cc_divs["m_start"], 1]),
                ],
                y=lab_ypos,
            )
        )
        time_scatter = (time_scatter
        + geom_vline(aes(xintercept="starts"), linetype="dashed", data=cc_phase)
        + geom_text(aes(x="labpos", y="y", label="labels"), data=cc_phase))

    return time_scatter

def scatter_pseudotime(adata, y, facet = False, size = 1.5, alpha = 1, color = 'black'):
    """DEPRECATED: use pseudotime_scatter
    """
    warnings.warn("scatter_pseudotime is deprecated; use pseudotime_scatter", DeprecationWarning)
    return (pseudotime_scatter(adata))
