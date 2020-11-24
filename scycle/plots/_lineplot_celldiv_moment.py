#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from plotnine import ggplot, aes, geom_point, geom_path, geom_smooth, annotate, labs, geom_col, scale_fill_manual
from ._themes import theme_std

def lineplot_celldiv_moment(adata): 
    """ Plots total_counts as a function of the principal circle nodes to
    visualize the moment of cell division.
    
    Parameters
    ----------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.celldiv_moment`.
        
    Returns
    ------------
    A plotnine line-plot to help visualize the moment of cell division and
    direction of the cell cycle. The moment of cell division is defined by the 
    largest drop in total_counts. The changes in counts are represented by the
    bars at the bottom, and the suggested moment of cell division is marked in
    red. The cell cycle should follow an incremental increase in total counts
    until around the moment of cell division.
    """
    edge_to_0 = adata.uns['scycle']['cell_div_moment']['cell_div_edge'][0]
    edges = adata.uns['princirc_gr']['edges']
    edges['cell_div'] = edges['e1'] == edge_to_0
    ref_var = adata.uns['scycle']['cell_div_moment']['ref_var']
    
    cell_div_count = edges[edges['e1'] == edge_to_0]['mean_var']
        
    cell_div_plot = (ggplot(edges, aes('e1', 'mean_var'))
     + geom_point(aes(y = 'mean_var'), size = 2)
     + geom_path(aes(y = 'mean_var'))
     + geom_smooth(aes(y = 'mean_var'), method = 'lm', linetype = 'dashed')
     + annotate("point", x = edge_to_0, y = cell_div_count, color = 'red', size = 2)
     + labs(x = 'Edge position', y = ref_var)
     + geom_col(aes(y = 'diff_var', fill = 'cell_div'))
     + scale_fill_manual(values = ['darkgrey', 'red'], guide = False)
     + theme_std)
    
    return cell_div_plot