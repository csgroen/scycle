#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotnine import ggplot, aes, geom_bar, coord_flip, guides, labs, theme_light, theme, element_blank, element_line

def barplot_cycle_phase(adata):
    """Plots the proportion of cells in each phase of the cell cycle
    
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
    plt_data['cell_cycle_phase'] = pd.Categorical(plt_data['cell_cycle_phase'], categories = ['M', 'G2', 'S', 'G1'])
    
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
    )
    
    return cycle_plot