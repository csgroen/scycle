#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import elpigraph
from ..data import g1s_markers, g2m_markers

def find_nonproliferating_cells(adata,
                                nonprolif_frac=0.3, 
                                n_nodes=50,
                                max_iter = 20,
                                n_sd = 3.0,
                                Mu = 1.0,
                                verbose = True):
    """ Find non-proliferating cells from the trajectory on cell cycle space
    
    Parameters
    ----------------------
    adata: AnnData
        AnnData to be annotated. Should be previously pre-processed.
    nonprolif_frac: float
        Initial estimation of fraction of non-proliferating cells
    n_nodes: int
        Number of nodes used for finding the trajectory
    max_iter: float
        Maximum number of iterations of new trajectories to find non-proliferating
        cells
    n_sd: float
        Number of standard deviations of distance that non-proliferating cells
        deviate from trajectory
    Mu: float
        Parameter passed to elpigraph.computeElasticPrincipalCircle
    verbose: bool
        If True, the function will print messages
    """
    all_markers = g1s_markers + g2m_markers
    Xccm = adata[:,list(set(all_markers)&set(adata.var_names))].X
    cc_score = np.array(list(np.mean(Xccm,axis=1)))

    ind_sorted_prolif = np.argsort(cc_score)
    ind_nonprolif = ind_sorted_prolif[0:int(len(adata)*nonprolif_frac)]
    adata.obs['proliferating'] = np.empty(len(adata)).astype(np.bool)
    adata.obs['proliferating'][:] = True
    adata.obs['proliferating'][ind_nonprolif] = False
    
    # sc.pl.scatter(adata,x='G1-S',y='G2-M',color='proliferating')

    fraction_nonprolif_old = nonprolif_frac

    for i in range(max_iter):
        X_elpigraph_training = adata.obs[['G1-S','G2-M']].to_numpy().astype(np.float64)
        u = X_elpigraph_training.copy()
        X_elpigraph_training = X_elpigraph_training[adata.obs['proliferating'],:]

        egr = elpigraph.computeElasticPrincipalCircle(X_elpigraph_training,n_nodes,
                                                      Mu=Mu,verbose=False)
        partition, dists = elpigraph.src.core.PartitionData(X = u, NodePositions = egr[0]['NodePositions'], 
                                                                MaxBlockSize = 100000000, TrimmingRadius = np.inf,
                                                                SquaredX = np.sum(u**2,axis=1,keepdims=1))
        ind_prolif = adata.obs['proliferating']
        mndist = np.mean(dists[ind_prolif])
        intervaldist = np.std(dists[ind_prolif])*n_sd
        tt1 = [not b for b in adata.obs['proliferating']]
        tt2 = [(d>mndist+intervaldist)[0] for d in dists]
        nonprolif_new = np.array(tt1) & np.array(tt2)
        adata.obs['proliferating'] = [not b for b in nonprolif_new]
        fraction_nonprolif = 1-np.sum(adata.obs['proliferating'])/len(adata)
        if verbose:
            print('\n\n===========\nIteration',i,'Fraction of non-proliferating cells:',fraction_nonprolif,'\n==============\n\n\n')    
        if np.abs(fraction_nonprolif-fraction_nonprolif_old)<0.01:
            break
        fraction_nonprolif_old = fraction_nonprolif
