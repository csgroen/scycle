#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotnine import (
    ggplot,
    aes,
    geom_point,
    geom_path,
    geom_text,
    scale_color_cmap,
    labs,
)
from ._themes import theme_std

def cell_cycle_projection(
    adata,
    col_var='total_counts',
    shape_var=None,
    palette='viridis',
    size=1.5,
    alpha=0.7,
    trajectory=False,
    node_size=7,
    node_color='lightgrey',
    show_nid=True
):
    """Plots the 2-D projection of the cells in the G1-S and G2-M dimensions.

    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
    shape_var: str
        The variable to be mapped to the shape of the points. Must be present in adata.obs
    palette: str
        A `cmap` palette to be used for coloring the scatterplot.
    size: float
        Controls the size of the points of the scatterplot.
    alpha: float
        Controls the transparency of the points of the scatterplot.
        [0 = completly transparent, 1 = completly opaque]
    trajectory: bool
        If True and the principal circle has already been calculated, the trajectory
        is added.
    node_size: float
        Controls the size of the node points of the principal circle.
    node_color: str
        If 'total_counts', shows the node average total_counts. Otherwise, must
        be a supported color name. e.g. node_color = 'black'.
    show_nid: bool
        If True, shows the node identified by the node points.

    Returns
    --------------
    A plotnine scatter plot of the 2-D projection on G1-S vs G2-M of all cells.
    """
    ax_names = ['G1/S', 'G2/M']

    plot_data = adata.obs[['G1-S', 'G2-M']]
    plot_data[col_var] = adata.obs[col_var].values

    if shape_var is None:
        proj_plot = (
            ggplot(plot_data, aes("G1-S", "G2-M"))
            + geom_point(aes(color=col_var), size=size, alpha=alpha)
            + theme_std
            + labs(x=ax_names[0], y=ax_names[1])
        )
    else:
        plot_data[shape_var] = adata.obs[shape_var].values
        proj_plot = (
            ggplot(plot_data, aes("G1-S", "G2-M"))
            + geom_point(aes(color=col_var, shape=shape_var), size=size, alpha=alpha)
            + theme_std
            + labs(x=ax_names[0], y=ax_names[1])
        )

    if trajectory and "principal_circle" in adata.uns["scycle"].keys():
        # Get the circle coordinates
        node_coords = adata.uns["princirc_gr"]["node_coords"]
        edge_coords = adata.uns["princirc_gr"]["edge_coords"]

        # Add to plot
        proj_plot = (
            proj_plot
            + geom_path(aes(x="G1-S", y="G2-M"), data=edge_coords)
            + geom_point(
                aes(x="G1-S", y="G2-M"),
                color=node_color,
                size=node_size,
                data=node_coords,
            )
        )

        if show_nid:
            proj_plot = proj_plot + geom_text(
                aes("G1-S", "G2-M", label="npos"), data=node_coords, size=10
            )
    return proj_plot

def cell_cycle_pca(
    adata,
    col_var="total_counts",
    shape_var=None,
    palette="viridis",
    size=1.5,
    alpha=0.7,
    trajectory=False,
    node_size=7,
    node_color="lightgrey",
    show_nid=True,
):
    """Plots the 2-D projection of the cells in the dimensions found by
    `tl.dimensionality_reduction`.

    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
    shape_var: str
        The variable to be mapped to the shape of the points. Must be present in adata.obs
    palette: str
        A `cmap` palette to be used for coloring the scatterplot.
    size: float
        Controls the size of the points of the scatterplot.
    alpha: float
        Controls the transparency of the points of the scatterplot.
        [0 = completly transparent, 1 = completly opaque]
    trajectory: bool
        If True and the principal circle has already been calculated, the trajectory
        is added.
    node_size: float
        Controls the size of the node points of the principal circle.
    node_color: str
        If 'total_counts', shows the node average total_counts. Otherwise, must
        be a supported color name. e.g. node_color = 'black'.
    show_nid: bool
        If True, shows the node identified by the node points.

    Returns
    --------------
    A plotnine scatter plot of the 2-D projection of all cells.
    """
    # Get coordinates
    X_dimRed = adata.obsm["X_pca_scycle"][:, 0:2]

    # Set axes names
    dimRed_method = adata.uns["scycle"]["dimRed"]["method"]

    if dimRed_method == "cc_signatures":
        ax_names = ["G1-S", "G2-M"]
    else:
        ax_names = ["PC 1", "PC 2"]

    plot_data = pd.DataFrame(X_dimRed, columns=["PCA", "PCB"])
    plot_data[col_var] = adata.obs[col_var].values

    # Plot
    if shape_var is None:
        proj_plot = (
            ggplot(plot_data, aes("PCA", "PCB"))
            + geom_point(aes(color=col_var), size=size, alpha=alpha)
            + theme_std
            + labs(x=ax_names[0], y=ax_names[1])
        )
    else:
        plot_data[shape_var] = adata.obs[shape_var].values
        proj_plot = (
            ggplot(plot_data, aes("PCA", "PCB"))
            + geom_point(aes(color=col_var, shape=shape_var), size=size, alpha=alpha)
            + theme_std
            + labs(x=ax_names[0], y=ax_names[1])
        )

    # Palettes
    if palette != "viridis":
        proj_plot = proj_plot + scale_color_cmap(cmap_name=palette)

    if trajectory and "principal_circle" in adata.uns["scycle"].keys():
        # Get the circle coordinates
        node_coords = adata.uns["princirc_gr"]["node_coords"]
        edge_coords = adata.uns["princirc_gr"]["edge_coords"]

        # Add to plot
        if node_color == "total_counts":
            proj_plot = (
                proj_plot
                + geom_path(aes(x="x", y="y"), data=edge_coords)
                + geom_point(
                    aes(x="x", y="y", color="total_counts"),
                    size=node_size,
                    data=node_coords,
                )
            )
        else:
            proj_plot = (
                proj_plot
                + geom_path(aes(x="x", y="y"), data=edge_coords)
                + geom_point(
                    aes(x="x", y="y"),
                    color=node_color,
                    size=node_size,
                    data=node_coords,
                )
            )

        if show_nid:
            proj_plot = proj_plot + geom_text(
                aes("x", "y", label="npos"), data=node_coords, size=10
            )

    return proj_plot

import warnings
def scatter_projection(
    adata,
    col_var="total_counts",
    shape_var=None,
    palette="viridis",
    size=1.5,
    alpha=0.7,
    trajectory=False,
    node_size=7,
    node_color="lightgrey",
    show_nid=True,
):
    """DEPRECATED. see: cell_cycle_pca

    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
    shape_var: str
        The variable to be mapped to the shape of the points. Must be present in adata.obs
    palette: str
        A `cmap` palette to be used for coloring the scatterplot.
    size: float
        Controls the size of the points of the scatterplot.
    alpha: float
        Controls the transparency of the points of the scatterplot.
        [0 = completly transparent, 1 = completly opaque]
    trajectory: bool
        If True and the principal circle has already been calculated, the trajectory
        is added.
    node_size: float
        Controls the size of the node points of the principal circle.
    node_color: str
        If 'total_counts', shows the node average total_counts. Otherwise, must
        be a supported color name. e.g. node_color = 'black'.
    show_nid: bool
        If True, shows the node identified by the node points.

    Returns
    --------------
    A plotnine scatter plot of the 2-D projection of all cells.
    """
    warnings.warn("scatter_projection is deprecated; use cell_cycle_pca", DeprecationWarning)
    return (cell_cycle_pca(adata, col_var, shape_var, palette, size, alpha, trajectory,
    node_size, node_color, show_nid))
