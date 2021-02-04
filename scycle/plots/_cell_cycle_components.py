import pandas as pd
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine.scales import scale_alpha_manual, scale_size_manual, scale_color_manual
from plotnine.facets import facet_wrap
import matplotlib as mpl
import numpy as np
from ._themes import theme_std
import warnings


def cell_cycle_components(adata, plot_type = 'panel', palette = 'Set1'):
    """Plots a scatter plot of trajectory vs component scores for each component
    from the dimensionality reduction, highlighting the selected cell cycle components

    Parameters
    --------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.find_cc_components` and `tl.principal_circle`.
    plot_type: str
        One of 'all' or 'panel'
    palette: 'str'
        A palette supported by matplotlib.cm.get_cmap
    
    Returns
    -------------
    A plotnine scatter plot of IC scores vs trajectory. It can be used to
    diagnose whether the cell cycle ICs vary through the trajectory, and if
    others do not.
    """
    #-- Get projection data
    proj = adata.obsm['X_dimRed']
    n_ics = proj.shape[1]
    spart = adata.obs['partition'].values
    comps = adata.uns['scycle']['find_cc_components']['indices']
    comps_found = list(comps.keys())
    
    #-- Make IC dataframe
    ic_df = pd.DataFrame(proj)
    ic_names = ['IC'+str(i) for i in range(proj.shape[1])]
    ic_df.columns = ic_names
    ic_df['partition'] = spart
    
    #-- Melt for plotting
    ic_traj = ic_df.groupby('partition').sum()
    # ic_traj = pd.DataFrame(zscore(ic_traj))
    ic_traj['partition'] = [i for i in range(np.max(spart)+1)]
    ic_trajm = pd.melt(ic_traj, id_vars = 'partition', var_name = 'IC')
    
    if plot_type == 'all':
        #-- Add variables for mapping plotting
        for cp in comps_found:
            ic_trajm = _update_ictrajm(ic_trajm, comps, cp)
        
        idx = [i not in comps_found for i in ic_trajm['IC']]
        ic_trajm['ccIC'] = 'cell cycle IC'
        ic_trajm.loc[idx,'ccIC'] = 'other'
        
        #-- Get colors
        cmap = mpl.cm.get_cmap(palette, n_ics)
        colors = np.array([ mpl.colors.rgb2hex(cmap(i)) for i in range(n_ics) ])
        jmp = int(np.round(n_ics / 3))
        cidx = np.array([0, 0 + jmp, 0 + 2*jmp, n_ics-1])
        oidx = np.array([i not in cidx for i in range(n_ics)])
        cc_cols = np.append(colors[cidx], colors[oidx])
        
        #-- Plot
        splot = (ggplot(ic_trajm, aes(x='partition', y='value', color ='IC', alpha = 'ccIC', size = 'ccIC'))
         + geom_point(size = 3)
         + geom_line()
         + scale_alpha_manual(values = [1,0.2], name = 'IC type')
         + scale_size_manual(values = [1.5,1], name = 'IC type')
         + scale_color_manual(values = cc_cols)
         + labs(x = 'Trajectory', y = 'IC score')
         + theme_std)
    
    elif plot_type == 'panel':
        #-- Add variables for mapping plotting
        if 'Histone' and 'G2/M-' in comps_found:
            ic_trajm1 = _multi_ictrajm(ic_trajm, comps, 'G1/S')
            ic_trajm2 = _multi_ictrajm(ic_trajm, comps, 'G2/M')
            ic_trajm3 = _multi_ictrajm(ic_trajm, comps, 'G2/M-')
            ic_trajm4 = _multi_ictrajm(ic_trajm, comps, 'Histone')
            ic_trajm4plot = pd.concat([ic_trajm1, ic_trajm2, ic_trajm3, ic_trajm4])
        elif 'Histone' in comps_found:
            ic_trajm1 = _multi_ictrajm(ic_trajm, comps, 'G1/S')
            ic_trajm2 = _multi_ictrajm(ic_trajm, comps, 'G2/M')
            ic_trajm4 = _multi_ictrajm(ic_trajm, comps, 'Histone')
            ic_trajm4plot = pd.concat([ic_trajm1, ic_trajm2, ic_trajm4])
        elif 'G2/M-' in comps_found:
            ic_trajm1 = _multi_ictrajm(ic_trajm, comps, 'G1/S')
            ic_trajm2 = _multi_ictrajm(ic_trajm, comps, 'G2/M')
            ic_trajm3 = _multi_ictrajm(ic_trajm, comps, 'G2/M-')
            ic_trajm4plot = pd.concat([ic_trajm1, ic_trajm2, ic_trajm3])
        else:
            ic_trajm1 = _multi_ictrajm(ic_trajm, comps, 'G1/S')
            ic_trajm2 = _multi_ictrajm(ic_trajm, comps, 'G2/M')
            ic_trajm4plot = pd.concat([ic_trajm1, ic_trajm2])

        #-- Get mapping colors
        cmap = mpl.cm.get_cmap(palette, 5)
        cc_cols = np.append(
            np.array([ mpl.colors.rgb2hex(cmap(i)) for i in range(4) ]), 
            'grey')
        
        #-- Plot
        splot = (ggplot(ic_trajm4plot, aes(x='partition', y='value', color ='IC', alpha = 'IC', size = 'IC'))
             + facet_wrap(facets='facet')
             + geom_point(size = 3)
             + geom_line()
             + scale_size_manual(values = [1.5,1.5,1.5,1.5,1])
             + scale_alpha_manual(values = [1,1,1,1,0.2])
             + scale_color_manual(values = cc_cols)
             + theme_std
             + labs(x = 'Trajectory', y = 'IC score'))
        
    return splot

def _update_ictrajm(ic_trajm, comps, var):
    ic_trajm.loc[ic_trajm['IC'] == 'IC' + str(comps[var]),'IC'] = var
    return ic_trajm

def _multi_ictrajm(icdf, comps, var):
    icdf = icdf.copy()
    icdf['facet'] = var + ' signature'
    idx = icdf['IC'] == 'IC' + str(comps[var])
    icdf.loc[idx, 'IC'] = var
    icdf.loc[np.invert(idx),'IC'] = 'other'
    return icdf

def scatter_enrich_components(adata, plot_type = 'panel', palette = 'Set1'):
    """DEPRECATED: use cell_cycle_components
    """
    warnings.warn("scatter_enrich_components is deprecated; use cell_cycle_components", DeprecationWarning)
    return (cell_cycle_components(adata, plot_type, palette))
    

