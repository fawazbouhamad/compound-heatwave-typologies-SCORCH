# SPDX-License-Identifier: GPL-3.0-only
"""Stage 5: compound events, their four types and the manuscript summaries.

A compound event is a run of consecutive selected days. Its type depends on
how many clusters (ellipses) there were on each of its days:

* Type 1, Widespread (Isolated): one day with one cluster.
* Type 2, Spatially Clustered: one day with several clusters.
* Type 3, Temporally Clustered: several days, with one cluster on every day
  or with several clusters on every day.
* Type 4, Compound Clustering (Multi-Type): several days that mix
  one-cluster days and several-cluster days.

Every selected day has at least one cluster, so an event whose day has no
cluster is a data error and is rejected rather than classified.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TYPE_NAMES = {1: "Widespread (Isolated)", 2: "Spatially Clustered",
              3: "Temporally Clustered", 4: "Compound Clustering (Multi-Type)"}

TYPE_DESCRIPTIONS = {
    1: "a single selected day with one heatwave cluster",
    2: "a single selected day with two or more separate clusters",
    3: "two or more consecutive days that all have one cluster, or all have several",
    4: "two or more consecutive days that mix one-cluster and several-cluster days",
}


def classify_event(daily_cluster_counts):
    """Type 1-4 from the number of clusters on each day of an event.

    Type 1: one day, one cluster.  Type 2: one day, several clusters.
    Type 3: several days, either one cluster every day or several clusters
    every day.  Type 4: several days mixing one-cluster and several-cluster
    days.

    Raises ValueError for an empty event or a day without any cluster: such
    an event cannot occur in a consistent catalog and is not classified.
    """
    counts = [int(c) for c in daily_cluster_counts]
    if not counts:
        raise ValueError("an event needs at least one day; empty events are not classified")
    if any(c < 1 for c in counts):
        raise ValueError(f"every event day needs at least one cluster; got daily counts {counts}")
    if len(counts) == 1:
        return 1 if counts[0] == 1 else 2
    all_one = all(c == 1 for c in counts)
    all_multi = all(c > 1 for c in counts)
    return 3 if (all_one or all_multi) else 4


# Axial statistics and trend tests
def axial_mean_and_resultant(azimuths_deg):
    """Doubled-angle (axial) mean of orientations and its resultant length.

    mu = 0.5 * atan2(mean sin 2theta, mean cos 2theta), wrapped to (-90, 90];
    R = hypot(mean cos 2theta, mean sin 2theta). Axes 180 degrees apart are
    the same axis, so an ordinary arithmetic mean is not used.
    """
    arr = np.asarray(azimuths_deg, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        raise ValueError("axial mean needs at least one orientation")
    two = np.deg2rad(2.0 * arr)
    c = float(np.sum(np.cos(two)) / arr.size)
    s = float(np.sum(np.sin(two)) / arr.size)
    r = float(np.hypot(c, s))
    if r <= 1e-12:
        raise ValueError("axial mean is undefined for these orientations")
    mu = 0.5 * float(np.degrees(np.arctan2(s, c)))
    mu = float(np.mod(mu + 90.0, 180.0) - 90.0)
    return (mu + 180.0 if mu <= -90.0 else mu), r


def mann_kendall_test(y):
    """Two-sided Mann-Kendall test with tie-corrected variance and continuity
    correction."""
    from scipy.stats import norm
    y = np.asarray(y, dtype=float)
    n = len(y)
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(y[j] - y[k])
    _, counts = np.unique(y, return_counts=True)
    var_s = n * (n - 1) * (2 * n + 5)
    var_s = (var_s - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18.0
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    if p < 0.05 and z > 0:
        trend = "increasing"
    elif p < 0.05 and z < 0:
        trend = "decreasing"
    else:
        trend = "no significant trend"
    return s, z, p, trend


def sens_slope(x, y):
    """Median of all pairwise slopes."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slopes = []
    for i in range(len(x) - 1):
        for j in range(i + 1, len(x)):
            if x[j] != x[i]:
                slopes.append((y[j] - y[i]) / (x[j] - x[i]))
    return np.median(slopes)


# Manuscript summaries
def table_1(events):
    """Compound heatwave statistics by type (manuscript Table 1).

    A type without events gets a row with zero events and empty statistics
    instead of stopping the run.
    """
    rows = []
    for t in (1, 2, 3, 4):
        sub = events[events["event_type"] == t]
        row = {
            "compound_heatwave_type": f"Type {t}: {TYPE_NAMES[t]}",
            "selected_event_days": int(sub["n_selected_days"].sum()) if len(sub) else 0,
            "number_of_compound_events": int(len(sub)),
        }
        if len(sub):
            row.update({
                "longest_event_days": int(sub["duration_days"].max()),
                "mean_duration_days_per_event": round(float(sub["duration_days"].mean()), 2),
                "mean_ellipses_per_event": round(float(sub["total_clusters"].mean()), 2),
                "minimum_ellipses_per_event": int(sub["total_clusters"].min()),
                "maximum_ellipses_per_event": int(sub["total_clusters"].max()),
            })
        else:
            row.update({
                "longest_event_days": pd.NA,
                "mean_duration_days_per_event": pd.NA,
                "mean_ellipses_per_event": pd.NA,
                "minimum_ellipses_per_event": pd.NA,
                "maximum_ellipses_per_event": pd.NA,
            })
        rows.append(row)
    table = pd.DataFrame(rows)
    # Keep whole-number columns as integers when every type has events, as in
    # the manuscript table; otherwise pandas' nullable integers are used.
    for column in ("longest_event_days", "minimum_ellipses_per_event",
                   "maximum_ellipses_per_event"):
        table[column] = pd.array(table[column], dtype="Int64")
    if not table["longest_event_days"].isna().any():
        for column in ("longest_event_days", "minimum_ellipses_per_event",
                       "maximum_ellipses_per_event"):
            table[column] = table[column].astype(int)
    for column in ("mean_duration_days_per_event", "mean_ellipses_per_event"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return table


MINIMUM_TREND_YEARS = 3


def trend_statistics(events):
    """Trend tests on the annual mean duration of Type 3 and Type 4 events.

    The year of an event is the year of its first day. Only years with at
    least one event of the type enter the series. Annual event counts are
    reported descriptively and are not trend-tested. A type with fewer than
    three such years gets a row of missing values (no test is possible).
    """
    from scipy.stats import linregress
    rows, series = [], []
    for t in (3, 4):
        sub = (events[events["event_type"] == t].groupby("year")["duration_days"]
               .mean().reset_index().sort_values("year"))
        x = sub["year"].to_numpy(dtype=float)
        y = sub["duration_days"].to_numpy(dtype=float)
        if len(sub) < MINIMUM_TREND_YEARS:
            rows.append({
                "event_type": t, "event_type_name": TYPE_NAMES[t], "n_years": len(sub),
                "ols_slope_days_per_year": np.nan, "ols_slope_days_per_decade": np.nan,
                "ols_intercept": np.nan, "ols_p_value": np.nan, "ols_r_value": np.nan,
                "ols_r_squared": np.nan, "sen_slope_days_per_year": np.nan,
                "sen_slope_days_per_decade": np.nan, "mann_kendall_s": np.nan,
                "mann_kendall_z": np.nan, "mann_kendall_p_value": np.nan,
                "mann_kendall_result": f"not computed (fewer than {MINIMUM_TREND_YEARS} years)",
                "residual_mae_days": np.nan, "residual_rmse_days": np.nan,
            })
            series.append(pd.DataFrame({"event_type": t, "year": sub["year"],
                                        "mean_duration_days": y,
                                        "ols_fitted_days": np.full(len(y), np.nan)}))
            continue
        ols = linregress(x, y)
        fitted = ols.intercept + ols.slope * x
        residual = y - fitted
        s, z, p, trend = mann_kendall_test(y)
        sen = sens_slope(x, y)
        rows.append({
            "event_type": t, "event_type_name": TYPE_NAMES[t], "n_years": len(sub),
            "ols_slope_days_per_year": ols.slope,
            "ols_slope_days_per_decade": ols.slope * 10,
            "ols_intercept": ols.intercept, "ols_p_value": ols.pvalue,
            "ols_r_value": ols.rvalue, "ols_r_squared": ols.rvalue ** 2,
            "sen_slope_days_per_year": sen, "sen_slope_days_per_decade": sen * 10,
            "mann_kendall_s": s, "mann_kendall_z": z, "mann_kendall_p_value": p,
            "mann_kendall_result": trend,
            "residual_mae_days": float(np.mean(np.abs(residual))),
            "residual_rmse_days": float(np.sqrt(np.mean(residual ** 2))),
        })
        series.append(pd.DataFrame({"event_type": t, "year": sub["year"],
                                    "mean_duration_days": y,
                                    "ols_fitted_days": fitted}))
    return pd.DataFrame(rows), pd.concat(series, ignore_index=True)
