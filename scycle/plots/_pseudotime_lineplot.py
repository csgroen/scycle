#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from plotnine import ggplot, aes, geom_point, labs, facet_wrap, geom_vline, geom_text, geom_smooth
from ._themes import theme_std
import warnings

# import scycle as cc
# adata = cc.tl.read('/home/clarice/Desktop/test.zip')
# ncol = 2; alpha = 1; size = 1; color = 'black'; smoothness = 0.3
y = ['TOP2A', 'CCNB1', 'CCNA2']

def pseudotime_lineplot(adata, y, smoothness = 0.3, facet = True, alpha = 1,
                       size = 1, color = 'black', ncol = 2, lab_ypos = 2):
    """Plots a line plot of pseudotime vs one or multiple variables

    Parameters
    --------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.pseudotime`.
    y: str or list
        If type(y) == str, y must be a variable annotated in adata.obs and
        will be used as the y-axis. If type(y) == list, then multiple variables
        will be plotted using a shared y-axis but different point colors.
    span: float
        Controls the smoothness of the regression. See: span in `plotnine.geom_smooth()`
    facet: bool
        Whether to return a facetted plot or all signatures in a single plot.
        Only used if y is a list.
    alpha: float
        A value between 0 and 1. Controls point transparency.
    size: float
        Controls the line width of the smooth line.
    color: str
        A supported color name. Controls the point color if type(y)==str.
        Ignored otherwise.
    ncol: int
        Number of columns in the facetting, if facet=True. Ignored otherwise.
    lab_ypos: float
        Controls the y-axis position of the cell cycle phase annotation, if present.

    Returns
    -------------
    A plotnine line plot of pseudotime.
    """
    if type(y) == str:
        #-- Get data
        if y in adata.obs.columns:
            plot_df = pd.DataFrame({'x': adata.obs['pseudotime'], 'y': adata.obs[y]})
        elif y in adata.var_names:
            plot_df = pd.DataFrame({'x': adata.obs['pseudotime'], 'y': adata[:,y].X.flatten()})
        else:
            raise Exception('`y` variable not found')

        #-- Make plot
        if color in adata.obs.columns:
            time_line = (ggplot(plot_df, aes(x = 'x', y = 'y'))
              + geom_smooth(aes(color = color), method = 'mavg', size = size, alpha = alpha, span = smoothness, se = False)
              + labs(x = 'Pseudotime', y = y)
              + theme_std)
        else:
            time_line = (ggplot(plot_df, aes(x = 'x', y = 'y'))
              + geom_smooth(method = 'mavg', size = size, alpha = alpha, color = color, span = smoothness, se = False)
              + labs(x = 'Pseudotime', y = y)
              + theme_std)

    else:
        #-- Make multiple color plot
        sannot = pd.DataFrame({'pseudotime': adata.obs['pseudotime']})
        sannot['id'] = range(sannot.shape[0])
        #-- Checks
        check1 = [var in adata.var_names for var in y]
        check2 = [var in adata.obs.columns.values for var in y]
        idx = np.array(check1) | np.array(check2)
        y_arr = np.array(y)
        if not np.any(idx):
            raise Exception('No variables in `y` found.')
        if not np.all(idx):
            warnings.warn('Variable not found! Dropping: ' + ', '.join((y_arr[~idx])))
            y = y_arr[idx]
        #-- Get y from obs or matrix:
        for var in y:
            if var in adata.obs.columns:
                sannot[var] = adata.obs[var]
            elif var in adata.var_names:
                sannot[var] = adata[:,var].X.flatten()
        plot_df = pd.melt(sannot, id_vars = ['id', 'pseudotime'],
                          var_name = 'signature', value_name = 'score')
        plot_df['signature'] = plot_df['signature'].astype('category')
        plot_df['signature'].cat.categories
        plot_df['signature'].cat.reorder_categories(y, inplace=True)

        if facet:
            time_line = (ggplot(plot_df, aes('pseudotime', 'score'))
            + facet_wrap('signature', scales = 'free_y', ncol = ncol)
            + geom_smooth(aes(color = 'signature'), method='mavg', size = size, span = smoothness, se = False)
            + theme_std)
        else:
            time_line = (ggplot(plot_df, aes('pseudotime', 'score'))
            + geom_smooth(aes(color = 'signature'), method='mavg', size = size, span = smoothness, se = False)
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
        time_line = (time_line
        + geom_vline(aes(xintercept="starts"), linetype="dashed", data=cc_phase)
        + geom_text(aes(x="labpos", y="y", label="labels"), data=cc_phase))

    return time_line
