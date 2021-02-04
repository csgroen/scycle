# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# import numpy as np
# import elpigraph
# from sklearn.decomposition import PCA

# from ._subtract_cc import subtract_cc

# def self_consistent_trajectory(adata,
#                                max_iter=10,
#                                eps_jaccard = 0.01,
#                                n_pc = 30,
#                                n_nodes = 30,
#                                r2_threshold = 0.5, 
#                                Mu=0.1,
#                                verbose = True):
#     """Uses the trajectory-determined cell cycle subspace to determine the trajectory
    
#     Used for `pp.normalize_by_partition`

#     Parameters
#     ------------
#     adata: AnnData
#         The analysis object to be evaluated. Must first be evaluated by
#         `tl.dimensionality_reduction`.
#     max_iter: int
#         Maximum number of iterations to find the cell cycle subspace from the
#         trajectory.
#     eps_jaccard: float
#         Parameter used for the convergence of the cell cycle subspace, compared
#         to previous iteration.
#     n_pc: int
#         Number of principal components of the cell cycle genes used as the
#         cell cycle subspace.
#     n_nodes: int
#         Number of nodes used for the calculation of the trajectory.
#     r2_threshold: float
#         R^2 threshold for cell-cycle related genes when regressed against the
#         trajectory.
#     Mu: float
#         Passed to elpigraph.computeElasticPrincipalCircle: the lambda 
#         parameter used the compute the elastic energy
#     verbose: bool
#         If True, the function will print messages.

#     Returns
#     ------------
#     `adata` will be updated with the partitions calculated from the self-consistent
#     trajectory
#     """
    
#     Xtemp = adata.X
    
#     subtract_cc(adata)
#     adata.X = Xtemp
#     r2scores = adata.var['r2_scores']
#     ind = np.where(r2scores>r2_threshold)[0]
#     if verbose:
#         print('Cell cycle genes initially found:',len(ind))
#         print('Top ones:',list(adata.var_names[np.where(r2scores>0.9)[0]]))    
    
#     ind_old = ind
#     for counter in range(max_iter):
#         Xcc = adata.X[:,ind]
#         pca = PCA(n_components=np.min((n_pc,len(ind)-1)))
#         Xcc_reduced = pca.fit_transform(Xcc)
#         Xr = Xcc_reduced.astype(np.float64)
        
#         egr = elpigraph.computeElasticPrincipalCircle(Xr,n_nodes,Mu=Mu,verbose=False)
        
#         nodep = egr[0]['NodePositions']
#         partition, dists = elpigraph.src.core.PartitionData(X = Xr, NodePositions = nodep, 
#                                                         MaxBlockSize = 100000000, TrimmingRadius = np.inf,
#                                                         SquaredX = np.sum(Xr**2,axis=1,keepdims=1))        
#         adata.obs['partition'] = partition
#         Xtemp = adata.X
#         subtract_cc(adata)
#         adata.X = Xtemp
#         r2scores_sc = adata.var['r2_scores']
        
#         ind = np.where(r2scores_sc>r2_threshold)[0]
#         idx_common_genes = list(set(ind)&set(ind_old))
#         union_genes = list(set(ind)|set(ind_old))
#         perc = len(idx_common_genes)/len(union_genes)
        
#         if verbose:
#             print('\n\nIteration',counter+1,'==================\nJaccard coeff:',perc,'Old:',len(ind_old),'New:',len(ind),'\n==================\n')
#         if perc>1-eps_jaccard:
#             break
#         ind_old = ind.copy()
    
#     if verbose:
#         print('Cell cycle genes found from the self-consistent trajectory:',len(ind))
#         print('Top ones:',list(adata.var_names[np.where(r2scores>0.9)[0]]))
#     adata.uns['scycle']['self-consistent_trajectory'] = {
#         'cc_genes': list(adata.var_names[ind]),
#         'Xr': Xr,
#         }