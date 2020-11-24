import re
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as mt
from collections import OrderedDict
from ..annot import pathway_annotation


def enrich_components(adata, collection = 'REACTOME', 
                      terms2select = ['^G1_S', '^G2', '^G2_M', 'HISTONE', '^DNA_REPLICATION', '^M_PHASE', 'CELL_CYCLE'], 
                      ntop = 5, nsd = 3, pcutoff = 0.05, verbose = True):
    """ Pathway enrichment for components from the dimensionality reduction
    
    Parameters
    --------------
    adata: AnnData
        AnnData object for the analysis. Must be previously evaluated by
        tl.dimensionality_reduction.
    collection: str
        One of 'REACTOME','KEGG' or 'BIOCARTA'. A pathway collection to be
        used for pathway enrichment.
    terms2select: list
        A list of regular expressions to be matched with enriched terms and
        select the components.
    ntop: int
        Number of top enriched terms to be used to search for cell cycle
        components.
    nsd: int
        Number of standard deviations for selection of genes associated with
        each component.
    pcutoff: float
        Cut-off of adjusted p-value to be used for the list of enriched terms
        for each component.
    verbose: bool
        If True, messages about function progress will be printed.
    
    Returns
    ------------------
    `adata` will be updated with the results of the component enrichment and
    automatic selection attempt using the function parameters.
    """
    
    #-- Get components
    comps = adata.uns['dimRed'].components_
    comps = pd.DataFrame(np.transpose(comps))
    n_comps = comps.shape[1]
    comps.index = adata.var.index
    
    #-- Do enrichment (hypergeometric)
    comp_enrich = {'comp'+str(i): _enrich_component(comps.iloc[:,i], collection, nsd, pcutoff) for i in range(n_comps)}
    comp_names = ['comp'+str(i) for i in range(n_comps)]
    
    top_enrich = {comp: comp_enrich[comp]['up_pathways'].head(ntop)['pathway'].values for comp in comp_names}
    bottom_enrich = {comp: comp_enrich[comp]['down_pathways'].head(ntop)['pathway'].values for comp in comp_names}
    
    #-- Search enrichment terms
    r = re.compile('|'.join(terms2select))
    us_comps = np.where([np.any([bool(r.search(term)) for term in top_enrich[comp]]) for comp in comp_names])[0]
    ds_comps = np.where([np.any([bool(r.search(term)) for term in bottom_enrich[comp]]) for comp in comp_names])[0]

    #-- Get suggested
    suggest_comps = list(us_comps) + list(ds_comps)
    suggest_comps_names = [str(i) for i in suggest_comps]
    ndim = len(suggest_comps)
    
    #-- Inform
    if verbose:
        print('Suggesting', ndim, 'components:\n', ', '.join(suggest_comps_names))
        for i in us_comps:
            cn='comp'+str(i)
            print('--', cn+':')
            print(top_enrich[cn])
        for i in ds_comps:
            cn='comp'+str(i)
            print('--', cn+':')
            print(bottom_enrich[cn])

    #-- Return
    adata.uns['scycle']['enrich_components'] = {'comp_enrich_results': comp_enrich,
                                               'up_enrich': top_enrich,
                                               'down_enrich': bottom_enrich,
                                               'suggested_comps': suggest_comps,
                                               'ndims': ndim}

def _enrich_component(component, collection, nsd, pcutoff):
    loadings = component.values
    up = np.mean(loadings) + nsd*np.std(loadings)
    down = np.mean(loadings) - nsd*np.std(loadings)
    
    #-- Get top genes contributing to component
    up_genes = component[component > up].sort_values(0, ascending = False).index.values
    down_genes = component[component < down].sort_values(0).index.values
    
    #-- Enrich
    up_res = __enrich_pathways(up_genes, source = collection, pcutoff = pcutoff)
    down_res = __enrich_pathways(down_genes, source = collection, pcutoff = pcutoff)
    
    return {'up_pathways': up_res, 'down_pathways': down_res}


def __enrich_pathways(genes, source, pcutoff = 0.05):
    paths = pathway_annotation(source)
    paths = OrderedDict(paths)

    #-- Get unique genes in pathway annotation (universe)
    unique_genes = np.unique([gene for path in paths.values() for gene in path])
    M = len(unique_genes)
    Np = np.sum([gene in unique_genes for gene in genes])

    #-- Get results
    enrich_res = [___enrich_path(genes, path, M, Np) for path in paths.values()]
    pvals = [res[0] for res in enrich_res]
    padj = mt.fdrcorrection(pvals)[1]
    path_ratios = [res[1] for res in enrich_res]
    universe_ratios = [str(Np)+'/'+str(M) for i in range(len(paths))]
    
    enrich_table = pd.DataFrame({'pathway': list(paths.keys()), 'pval': pvals, 
                                 'padj': padj, 'path_ratio': path_ratios, 
                                 'universe_ratio': universe_ratios})
    enrich_table = enrich_table[enrich_table['padj'] < pcutoff]
    enrich_table.sort_values('padj')
    return enrich_table

def ___enrich_path(genes, path, M, Np):
    #-- Calculate hypergeometric test
    n = len(path)
    # N = len(genes)
    k = np.sum([gene in path for gene in genes])
    pval = stats.hypergeom(M, n, Np).pmf(k)
    return (pval, str(k)+'/'+str(n))

