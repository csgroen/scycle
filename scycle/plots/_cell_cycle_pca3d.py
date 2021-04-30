import plotly.graph_objs as go
import plotly.express as px
import numpy as np
import pandas as pd
import warnings

def cell_cycle_pca3d(adata,
                     col_var = 'total_counts', 
                     size = 5,
                     palette = 'viridis', alpha = 0.7,
                     trajectory = False,  node_size = 10, 
                     node_color = 'lightgrey', show_nid = True): 
    """Plots the 3-D projection of the cells in the dimensions found by
    `tl.dimensionality_reduction`.
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
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
    A plotly.graph_objs scatter plot of the 3-D projection of all cells in the 
    3-D projection of the cell cycle dimensions.
    
    Note
    -------------
    We use plotly.graph_objs to render the 3-D plots. You may have to adjust
    the renderer. For using the browser as renderer you can use the following
    snippet:
        import plotly.io as pio
        pio.renderers.default = 'browser'
    """
    #-- Get cell projection coordinates
    assert 'X_pca_scycle' in adata.obsm, "3D PCA must be computed first."
    
    x = adata.obsm['X_pca_scycle'][:,0]
    y = adata.obsm['X_pca_scycle'][:,1]
    z = adata.obsm['X_pca_scycle'][:,2]
    col = adata.obs[col_var].values
    
    plotdf = pd.DataFrame(dict(x = x, y = y, z = z, col = col))
    
    fig = px.scatter_3d(plotdf, x='x',y='y',z='z', color='col', labels = {'col': col_var}, 
                            color_continuous_scale = palette)

    if type(col[0]) == str:
        fig.update_traces(marker = dict(size = size,
                                        opacity = alpha))
    else:
        fig.update_traces(marker = dict(size = size,
                                        opacity = alpha))
        
    if trajectory and 'principal_circle' in adata.uns['scycle'].keys():
        xn = adata.uns['princirc_gr']['node_coords']['x'].values
        xn = np.append(xn, xn[0])
        yn = adata.uns['princirc_gr']['node_coords']['y'].values
        yn = np.append(yn, yn[0])
        zn =  adata.uns['princirc_gr']['node_coords']['z'].values
        zn = np.append(zn, zn[0])
        npos = adata.uns['princirc_gr']['node_coords']['npos'].values
        npos = np.append(npos, 0)


        fig.add_traces(data = go.Scatter3d(x=xn, y=yn, z=zn,
                                           mode = 'lines+markers+text',
                                           marker = dict(
                                               color = node_color,
                                               size = node_size, 
                                               opacity = 1),
                                           line = dict(
                                               color = 'black'
                                               ),
                                           text = npos,
                                           name = 'trajectory'))
    if show_nid: 
        fig.update_layout(scene = dict(xaxis_title = 'PC1', 
                                       yaxis_title = 'PC2', 
                                       zaxis_title = 'PC3'))
    return fig


def scatter_projection3d(adata,
                         col_var = 'total_counts', 
                         size = 5,
                         palette = 'viridis', alpha = 0.7,
                         trajectory = False,  node_size = 10, 
                         node_color = 'lightgrey', show_nid = True): 
    """DEPRECATED. see: cell_cycle_pca3d
    
    Parameters
    -------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.dimensionality_reduction`
    col_var: str
        The variable to be used to color the points. Must be present in adata.obs
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
    A plotly.graph_objs scatter plot of the 3-D projection of all cells in the 
    3-D projection of the cell cycle dimensions.
    
    Note
    -------------
    We use plotly.graph_objs to render the 3-D plots. You may have to adjust
    the renderer. For using the browser as renderer you can use the following
    snippet:
        import plotly.io as pio
        pio.renderers.default = 'browser'
    """
    warnings.warn("scatter_projection3d is deprecated; use cell_cycle_pca3d", DeprecationWarning)
    return (cell_cycle_pca3d(adata, col_var, palette, size, alpha, trajectory,
    node_size, node_color, show_nid))
