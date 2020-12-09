#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import UnivariateSpline


def curvature(adata, smoothing_factor=10, s=0, k=3):
    """
    Computes principal circle curvature.
    adata: AnnData
        Annotated data of interest, with principal circle already computed.
    smoothing_factor: float
        Final curvature curve additional smoothing.
    s: float
        Positive smoothing factor used in scipy.UnivariateSpline
    k: float
        Degree of the spline, increasing it improves accuracy but
        can cause overfitting.
    """
    if "egr" not in adata.uns:
        print("Error: Principal circle must be computed before computing curvature.")
        return

    node_pos: np.ndarray = adata.uns["egr"]["NodePositions"]
    n_nodes: int = node_pos.shape[0]
    d_nodes: int = node_pos.shape[1]
    x: np.ndarray = np.linspace(0, n_nodes - 1, n_nodes)

    # Computing splines and derivatives
    splines = [UnivariateSpline(x, node_pos[:, i], s=s, k=k) for i in range(d_nodes)]
    dd_splines = [spline.derivative(n=2) for spline in splines]

    temp = np.array([deriv(x) for deriv in dd_splines])

    curv_raw = np.sqrt(np.sum(temp ** 2, axis=0))
    curv_spl = UnivariateSpline(x, curv_raw, s=np.var(curv_raw) * smoothing_factor, k=k)
    adata.uns["scycle"]["principal_circle"]["x_curv"] = x
    adata.uns["scycle"]["principal_circle"]["curv"] = curv_spl(x)
