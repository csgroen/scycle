#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from plotnine import aes, geom_path, geom_point, geom_text
from ._scatter_projection import scatter_projection

def scatter_projection_circle (adata, col_var = 'total_counts', palette = 'viridis', 
                       size = 1.5, node_size = 7, node_color = 'lightgrey', show_nid = True):
    """Plots the 2-D projection of cells with the principal circle nodes and
    edges overlayed.
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.principal_circle`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
    palette: str
        A `cmap` palette to be used for coloring the scatterplot.
    size: float
        Controls the size of the points of the scatterplot.
    node_size: floar
        Controls the size of the node points of the principal circle.
    node_color: str
        If 'total_counts', shows the node average total_counts. Otherwise, must
        be a supported color name. e.g. node_color = 'black'.
    show_nid: bool
        If True, shows the node identified by the node points.
        
    Returns
    --------------
    A plotnine scatter plot of the 2-D projection of all cells with the
    principal circle nodes and edges overlayed.
    """
    # Make the projection plot
    proj_plot = scatter_projection(adata, col_var = col_var, palette = palette, size = size)
    
    # Get the circle coordinates
    node_coords = adata.uns['princirc_gr']['node_coords']
    edge_coords = adata.uns['princirc_gr']['edge_coords']
    
    # Add to plot
    circle_plot = (proj_plot
     + geom_path(aes(x = 'x', y = 'y'), data = edge_coords)
     + geom_point(aes(x = 'x', y = 'y', color = 'total_counts'), 
                  size = node_size,  data = node_coords)
     )
    
    if node_color == 'total_counts':
            circle_plot = (proj_plot
             + geom_path(aes(x = 'x', y = 'y'), data = edge_coords)
             + geom_point(aes(x = 'x', y = 'y', color = 'total_counts'), 
                          size = node_size,  data = node_coords)
             )
    else:
            circle_plot = (proj_plot
             + geom_path(aes(x = 'x', y = 'y'), data = edge_coords)
             + geom_point(aes(x = 'x', y = 'y'), color = node_color,
                          size = node_size,  data = node_coords)
             )
    
    if show_nid:
        circle_plot = (circle_plot 
                       + geom_text(aes('x', 'y', label = 'npos'), 
                                   data = node_coords, size = 10))
    
    return circle_plot