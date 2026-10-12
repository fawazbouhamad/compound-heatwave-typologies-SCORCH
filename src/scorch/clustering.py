# SPDX-License-Identifier: GPL-3.0-only
"""Stage 3: DBSCAN clustering of the heatwave cells on each selected day.

Every selected day is clustered with every (eps, min_samples) pair of the
settings. The most frequent non-zero cluster count decides which pairs are
kept; the day's settings are the mean eps and the mean min_samples (rounded
half up) of the kept pairs. Within a sequence of consecutive selected days
the largest daily eps and min_samples are applied to every day.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

NOISE = -1   # DBSCAN label for cells that belong to no cluster


def longitude_major_order(lon, lat):
    """Sort cells by longitude, then latitude (west to east, south to north).

    The order matters for DBSCAN: a border cell within reach of two clusters
    joins the cluster that reaches it first. Labels follow the order in which
    clusters are first found.
    """
    return np.lexsort((np.asarray(lat), np.asarray(lon)))


def grid_distance_matrix(lon, lat, grid_step=1.0):
    """Euclidean distances between cells in grid-step units.

    One unit is one grid step (1 degree).
    """
    from scipy.spatial.distance import cdist
    pts = np.column_stack((np.asarray(lon, dtype=float) / float(grid_step),
                           np.asarray(lat, dtype=float) / float(grid_step)))
    return cdist(pts, pts, metric="euclidean").astype(float)


def run_dbscan(distances, eps, min_samples):
    """DBSCAN on a precomputed distance matrix (scikit-learn).

    A cell is a core cell when at least min_samples cells, itself included,
    lie within eps. Returns (labels, is_core): labels are 0, 1, 2, ... in
    order of first appearance and -1 for noise.
    """
    from sklearn.cluster import DBSCAN
    model = DBSCAN(eps=float(eps), min_samples=int(min_samples),
                   metric="precomputed").fit(distances)
    is_core = np.zeros(len(model.labels_), dtype=bool)
    is_core[model.core_sample_indices_] = True
    return model.labels_, is_core


def count_clusters(labels):
    return len(set(np.asarray(labels).tolist()) - {NOISE})


def parameter_search(distances, eps_values, min_samples_values):
    """Number of clusters for every (eps, min_samples) pair on one day."""
    rows = []
    for eps in eps_values:
        for min_samples in min_samples_values:
            labels, _ = run_dbscan(distances, eps, min_samples)
            rows.append({"eps": float(eps), "min_samples": int(min_samples),
                         "n_clusters": count_clusters(labels),
                         "n_noise_cells": int((labels == NOISE).sum())})
    return rows


def round_half_up(x):
    """Nearest integer with halves rounded up: floor(x + 0.5).

    This is the rounding used for the daily minimum sample size. It is not
    the same as Python's round(), which rounds halves to the nearest even
    number.
    """
    return int(np.floor(float(x) + 0.5))


def select_daily_parameters(search_rows):
    """Choose one day's DBSCAN settings from its parameter search.

    Find the most frequent nonzero cluster count (all counts if every pair
    gives zero clusters). If several counts tie, keep the pairs of all of them;
    daily eps is the mean retained eps and daily min_samples the mean
    retained min_samples rounded half up (never below 1).
    """
    if not search_rows:
        raise ValueError("the parameter search has no rows")
    counts = [r["n_clusters"] for r in search_rows]
    nonzero = [c for c in counts if c > 0]
    pool = nonzero if nonzero else counts
    frequency = pd.Series(pool).value_counts()
    top = int(frequency.max())
    modal_counts = sorted(int(v) for v, n in frequency.items() if n == top)
    kept = [r for r in search_rows if r["n_clusters"] in modal_counts]
    mean_min_samples = float(np.mean([r["min_samples"] for r in kept]))
    return {
        "modal_cluster_counts": modal_counts,
        "modal_frequency": top,
        "n_pairs_retained": len(kept),
        "daily_eps": float(np.mean([r["eps"] for r in kept])),
        "daily_min_samples_mean": mean_min_samples,
        "daily_min_samples": max(1, round_half_up(mean_min_samples)),
        "retained_pairs": {(r["eps"], r["min_samples"]) for r in kept},
    }


# DBSCAN summaries for Appendix A
DBSCAN_GROUPS = ["All Types", "Type 1", "Type 2", "Type 3", "Type 4"]


def values_by_group(frame, column):
    groups = {"All Types": frame[column].to_numpy(float)}
    for group in DBSCAN_GROUPS[1:]:
        groups[group] = frame.loc[frame["type"] == group, column].to_numpy(float)
    return groups


def box_statistics(values):
    """Box-plot statistics; quartiles by numpy's default linear interpolation.

    An empty group (a type without events) gives a zero count and NaN values.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan, "min": np.nan,
                "max": np.nan, "q1": np.nan, "q3": np.nan}
    return {
        "n": len(values),
        "mean": values.mean(),
        "median": float(np.median(values)),
        "min": values.min(),
        "max": values.max(),
        "q1": np.percentile(values, 25),
        "q3": np.percentile(values, 75),
    }


def figure_a2_panels(per_day, per_event):
    """(values by group, y-axis label, population) for panels (a)-(d)."""
    min_samples_label = "Selected DBSCAN min_samples"
    eps_label = "Selected DBSCAN eps (grid-cell units)"
    days_note = f"{len(per_day)} event-days"
    events_note = f"{len(per_event)} events"
    return {
        "a": (values_by_group(per_day, "min_samples"), min_samples_label, days_note),
        "b": (values_by_group(per_event, "min_samples"), min_samples_label, events_note),
        "c": (values_by_group(per_day, "eps"), eps_label, days_note),
        "d": (values_by_group(per_event, "eps"), eps_label, events_note),
    }


def statistics_table(panels):
    rows = []
    for letter, (groups, ylabel, population) in panels.items():
        for group in DBSCAN_GROUPS:
            rows.append({
                "panel": letter,
                "metric": ylabel,
                "population": population,
                "group": group,
                **box_statistics(groups[group]),
            })
    return pd.DataFrame(rows)
