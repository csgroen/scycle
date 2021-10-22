#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotnine import ggplot, aes, geom_bar, coord_flip, guides, labs, theme_light, theme, element_blank, element_line, scale_fill_brewer
import matplotlib.pyplot as plt
import warnings

def cell_cycle_phase_barplot(adata, palette='Set2'):
    """Plots the proportion of cells in each phase of the cell cycle

    See also: cell_cycle_phase_pieplot for the matplotlib pie chart


    Parameters
    -----------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.annotate_cell_cycle`.

    Returns
    -----------
    A plotnine barplot with the total counts of cell in each phase of the
    cell cycle.

    """
    plt_data = adata.obs.copy()
    plt_data['cell_cycle_phase'] = pd.Categorical(plt_data['cell_cycle_phase'], categories = ['G1 post-mitotic', 'G1 pre-replication', 'S/G2/M'])

    cycle_plot = (ggplot(plt_data, aes('cell_cycle_phase', fill = 'cell_cycle_phase'))
    + geom_bar()
    + coord_flip()
    + guides(fill = False)
    + labs(y = '', x = 'Cell cycle phase')
    + theme_light()
    + theme(panel_grid_major_y = element_blank(),
            panel_grid_minor_y = element_blank(),
            panel_grid_major_x = element_line(size = 1.5),
            panel_grid_minor_x = element_line(size = 1.5))
    + scale_fill_brewer(type='qual', palette=palette)
    )

    return cycle_plot

def cell_cycle_phase_pieplot(adata, colors = [ '#66c2a5', '#fc8d62', '#8da0cb']):
    """Plots the proportion of cells in each phase of the cell cycle

    See also: cell_cycle_phase_barplot for the plotnine barplot

    Parameters
    -----------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.annotate_cell_cycle`.
    colors: list of length 4
        A list of 4 colors to be used for the pie-chart for G1 post-mitotic,
        G1 pre-replication and S/G2/M, respectively.

    Returns
    -----------
    A matplotlib pieplot of the proportion of cells in each phase of the
    cell cycle.
    """
    colors.reverse()
    cell_counts = adata.obs.groupby('cell_cycle_phase').size()
    cell_counts = cell_counts.iloc[::-1]
    pie, ax = plt.subplots(figsize=[10,6])
    labels = cell_counts.keys()
    return (plt.pie(x=cell_counts, autopct="%.1f%%", explode=[0.05]*3,
                    labels=labels, colors = colors, pctdistance=0.5))

def barplot_cycle_phase(adata):
    """DEPRECATED: use cell_cycle_phase_barplot
    """
    warnings.warn("barplot_cycle_phase is deprecated; use cell_cycle_phase_barplot", DeprecationWarning)
    return cell_cycle_phase_barplot(adata)
