#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotnine import ggplot, aes, geom_point, scale_color_cmap, labs
# from ._themes import theme_std
from plotnine import theme_light, theme, element_blank
theme_std = theme_light() + theme(panel_grid=element_blank())

def scatter_projection (adata, project_main = True, comps = [0,1], 
                col_var = 'total_counts', palette = 'viridis', size = 1.5,
                alpha = 0.7):
    """Plots the 2-D projection of the cells in the dimensions found by
    `tl.dimensionality_reduction`.
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    project_main: bool
        If True, the first 2 components or the components defined using
        `tl.choose_components` will be plotted and `comps` will be ignored.
    comps: list
        If `project_main` is False, the components to be used for the x- and y-
        axes will be the components listed in `comps`.
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
    palette: str
        A `cmap` palette to be used for coloring the scatterplot.
    size: float
        Controls the size of the points of the scatterplot.
    alpha: float
        Controls the transparency of the points of the scatterplot.
        [0 = completly transparent, 1 = completly opaque]
        
    Returns
    --------------
    A plotnine scatter plot of the 2-D projection of all cells.
    """
    # Get coordinates
    if project_main:
        X_dimRed = adata.obsm['X_dimRed2d']
        dimRed_comps = [0,1]
    else:
        X_dimRed = adata.obsm['X_dimRed'][:, comps]
        dimRed_comps = comps
        
    # Set axes names    
    dimRed_method = adata.uns['scycle']['dimRed']['method']
    
    if dimRed_method == 'CCgenes':
        ax_names = ['G1-S', 'G2-M']
    else:
        ax_names = ['PC' + str(i+1) for i in  dimRed_comps]
        
    plot_data = pd.DataFrame(X_dimRed, columns = ['PCA', 'PCB'])
    plot_data[col_var] = adata.obs[col_var].values

    # Plot
    proj_plot = (ggplot(plot_data, aes('PCA', 'PCB'))
     + geom_point(aes(color = col_var), size = size, alpha = alpha)
     + theme_std
     + labs(x = ax_names[0], y = ax_names[1]))
        
    # Palettes
    if palette != 'viridis':
        proj_plot = (proj_plot + scale_color_cmap(cmap_name = palette))
    
    return proj_plot