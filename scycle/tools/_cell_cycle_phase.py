#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from ._curvature import curvature
from ._annotate_cell_cycle import annotate_cell_cycle


def cell_cycle_phase(adata,
                     smoothing_factor=10,
                     transition_refs = [0.35, 0.65],
                     max_refdist = 0.2,
                     annotate = True,
                     verbose = True):
    """Estimates the phase of cell cycle for each cell based on the curvature
    of the trajectory embedded in the cell cycle space

    Parameters
    -----------
    adata: AnnData
        The analysis object to be evaluated. Must first be evaluated by
        `tl.pseudotime`.
    smoothing_factor: int
        A smoothing factor used for calculating the Univariate Spline used to
        estimate the curvature of the trajectory
    transition_refs: list
        A list of length 3 with reference time-points for transitions between
        G1 and S, S and G2 and G2 and M, respectively.
    max_refdist: float
        Maximum distance between reference transition time and transition point.
        If no peak is within max_refdist of a transition reference point, the
        reference point is used instead.
    verbose: bool
        If True, the function will print messages.

    Returns
    ------------
    The `adata` object will be updated with the points of suggested start of
    cell cycle phases. To assign cell cycle phase to each cell, you need to
    run `tl.annotate_cell_cycle`.
    """

    # -- Get node positions
    node_coords = adata.uns["princirc_gr"]["node_coords"]
    idx = np.array(["dim" in cname for cname in node_coords.columns])

    nodep = np.array(node_coords.iloc[:, idx])
    nnodes = nodep.shape[0]

    # -- Calculate curvature and find peaks
    curvature(adata)
    curv_data = adata.uns["scycle"]["principal_circle"]
    x, curv = curv_data["x_curv"], curv_data["curv"]
    peaks = find_peaks(curv)[0]
    peak_times = peaks / nnodes

    # -- Check ref_times
    if len(transition_refs) != 2:
        raise Exception("`transition_refs` must be a list of length 2")

    #-- Get peaks closest to transition_refs
    sref_time, g2ref_time = transition_refs
    # sref_time, g2ref_time, mref_time = transition_refs
    pr_start = _transition_time(sref_time, peak_times, max_refdist)
    rep_start = _transition_time(g2ref_time, peak_times, max_refdist)
    # m_start = _transition_time(mref_time, peak_times, max_refdist)

    # -- Save curvature info
    curv_data = pd.DataFrame(dict(x=x, pseudotime=x / nnodes, curvature=curv)).merge(
        pd.DataFrame(dict(x=peaks, ispeak="peak")), on="x", how="left"
    )

    # -- Inform
    if verbose:
        print("-- Suggested cell cycle division:")
        print("G1 post-mitotic:", " 0  ", "-", pr_start)
        print("G1 pre-replication: ", pr_start, "-", rep_start)
        print("S/G2/M:", rep_start, "-", 1)

    # -- Return
    adata.uns["scycle"]["cell_cycle_division"] = {
        "pr_start": pr_start,
        "rep_start": rep_start,
        "curvature": curv_data,
    }

    if annotate:
        annotate_cell_cycle(adata)


def _transition_time(ref_time, peak_times, max_refdist):
    sdists = np.abs(peak_times - ref_time)
    idx = np.argmin(sdists)
    if sdists[idx] < max_refdist: return peak_times[idx]
    else: return ref_time
