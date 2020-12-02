#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotnine import ggplot, aes, geom_point, labs
from ._themes import theme_std

def scatter_pseudotime(adata, y, size = 1.5, alpha = 1, color = 'black'):
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
    size: float
        Controls the point size of the plot.
    alpha: float
        A value between 0 and 1. Controls point transparency.
    color: str
        A supported color name. Controls the point color if type(y)==str.
        Ignored otherwise.
        
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
        plt_vars = y
        plt_vars.append('pseudotime')
        sannot = adata.obs.copy()[plt_vars]
        sannot['id'] = range(sannot.shape[0])
            
        plot_df = pd.melt(sannot, id_vars = ['id', 'pseudotime'], 
                            var_name = 'signature', value_name = 'score')
    
        time_scatter = (ggplot(plot_df, aes('pseudotime', 'score'))
         + geom_point(aes(color = 'signature'), alpha = alpha, size = size)
         + theme_std)

    return time_scatter
