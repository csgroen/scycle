#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from plotnine import ggplot, aes, geom_histogram, labs
from ._themes import theme_std
import warnings

def pseudotime_histogram(adata, fill = '#595959', alpha = 1, bins = 30):
    """Plots a histogram of pseudotime
    
    Parameters
    --------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.pseudotime`.
    fill: str
        Controls the color of the histogram bars. Must be a supported color
        name or hex-code.
    alpha: float
        A float between 0 and 1. Controls the transparency of the bars.
        
    Returns
    -----------
    A plotnine histogram of pseudotime.
    """
    if fill in adata.obs.columns:
        hist_plt = (ggplot(adata.obs, aes('pseudotime', fill = fill))
                    + geom_histogram(alpha = alpha, bins = bins))
    
    else:
        hist_plt = (ggplot(adata.obs, aes('pseudotime'))
                    + geom_histogram(fill = fill, alpha = alpha, bins = bins))
        
    hist_plt = (hist_plt         
                + labs(x = 'Pseudotime', y = 'Count')
                + theme_std)
    
    return hist_plt

def hist_pseudotime(adata, fill = '#595959', alpha = 1, bins = 30):
    """DEPRECATED: use pseudotime_histogram
    
    Parameters
    --------------
    adata: AnnData
        The AnnData object being used for the analysis. Must be previously
        evaluated by `tl.pseudotime`.
    fill: str
        Controls the color of the histogram bars. Must be a supported color
        name or hex-code.
    alpha: float
        A float between 0 and 1. Controls the transparency of the bars.
        
    Returns
    -----------
    A plotnine histogram of pseudotime.
    """
    warnings.warn("hist_pseudotime is deprecated; use pseudotime_histogram", DeprecationWarning)
    return pseudotime_histogram(adata, fill = fill, alpha = alpha, bins = bins)