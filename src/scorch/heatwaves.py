# SPDX-License-Identifier: GPL-3.0-only
"""Stage 1: cell thresholds, exceedance and local heatwaves.

A grid cell is in a local heatwave when its daily maximum temperature (Tmax)
is at or above the cell's own warm-season percentile threshold on at least
``min_exceeding_days`` days within a run of at least ``min_duration_days``
days. Single non-exceeding days inside the run are bridged (they belong to
the heatwave); two consecutive non-exceeding days end it. Seasons are
labelled separately: 30 September and the following 1 April are not
consecutive days.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def cell_thresholds(tmax, quantile):
    """Empirical percentile of each cell's Tmax over all study dates.

    Use linear interpolation between ordered values, with float64 Tmax.
    """
    n_time, n_lat, n_lon = tmax.shape
    values = pd.DataFrame({
        "cell": np.repeat(np.arange(n_lat * n_lon), n_time),
        "val": tmax.reshape(n_time, -1).T.ravel(),
    })
    thresholds = values.groupby("cell", sort=True)["val"].quantile(quantile)
    return thresholds.to_numpy().reshape(n_lat, n_lon)


def label_heatwaves(exceed_series, min_len=3, min_ones=3,
                    allow_unlimited_single_zeros=True):
    """Label local heatwaves in one 0/1 daily exceedance series.

    * split the series wherever two consecutive days do not exceed,
    * trim each piece so it starts and ends on an exceeding day,
    * keep it when it spans at least ``min_len`` days and contains at least
      ``min_ones`` exceeding days; single non-exceeding days inside it are
      part of the heatwave.

    Returns (labels, episodes): the heatwave number of every day (0 = not in
    a heatwave) and a list of (number, start, end, length, n_exceeding,
    n_bridged) tuples.
    """
    s = pd.Series(exceed_series).astype(int)
    idx = pd.DatetimeIndex(s.index)
    vals = s.values
    zero = vals == 0

    # Position of the first zero of every "00" pair.
    double_zero_breaks = np.where(zero[:-1] & zero[1:])[0]
    boundaries = [-1] + double_zero_breaks.tolist() + [len(s) - 1]
    raw_segments = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        i0, i1 = a + 1, b
        if i0 <= i1:
            raw_segments.append((i0, i1))

    hw_id = 0
    labels = np.zeros(len(s), dtype=int)
    episodes = []
    for i0, i1 in raw_segments:
        seg_vals = vals[i0:i1 + 1]
        if seg_vals.sum() == 0:
            continue
        ones_pos = np.where(seg_vals == 1)[0]
        j0 = i0 + int(ones_pos[0])
        j1 = i0 + int(ones_pos[-1])
        event_vals = vals[j0:j1 + 1]
        ones_count = int(event_vals.sum())
        zeros_count = int((event_vals == 0).sum())
        length_days = int(j1 - j0 + 1)
        if allow_unlimited_single_zeros:
            ok = length_days >= min_len and ones_count >= min_ones
        else:
            ok = (length_days >= min_len and ones_count >= min_ones
                  and zeros_count <= 1)
        if ok:
            hw_id += 1
            labels[j0:j1 + 1] = hw_id
            episodes.append((hw_id, idx[j0].date(), idx[j1].date(),
                             length_days, ones_count, zeros_count))
    return labels, episodes


def label_cell_by_season(exceed, dates, min_len, min_ones):
    """Apply label_heatwaves to one cell, one warm season at a time.

    30 September and the following 1 April are not consecutive days, so a
    heatwave may not continue across that gap. Each season is labelled
    separately, and heatwave numbers keep counting upward through the
    record so every heatwave in a cell has its own number (1, 2, 3, ...).
    """
    years = np.array([int(d[:4]) for d in dates])
    labels = np.zeros(len(exceed), dtype=np.int32)
    episodes = []
    offset = 0
    for year in np.unique(years):
        rows = np.where(years == year)[0]
        season = pd.Series(exceed[rows], index=pd.to_datetime(dates[rows]))
        season_labels, season_episodes = label_heatwaves(
            season, min_len=min_len, min_ones=min_ones)
        labels[rows] = np.where(season_labels > 0, season_labels + offset, 0)
        for number, start, end, length, n_exceed, n_bridged in season_episodes:
            episodes.append((number + offset, start, end, length, n_exceed, n_bridged))
        offset += len(season_episodes)
    return labels, episodes


def detect_all_cells(tmax, thresholds, dates, settings):
    """Exceedance, heatwave membership and heatwave summaries for all cells.

    Returns
      exceed     int8 array (date, lat, lon): 1 when Tmax >= threshold
      number     int32 array (date, lat, lon): heatwave number within the
                 cell (0 = not in a heatwave), counted 1, 2, 3, ... through
                 the record in each cell
      episodes   list of per-heatwave rows
                 (lat_index, lon_index, number, start, end, length,
                 n_exceeding, n_bridged)
    """
    n_time, n_lat, n_lon = tmax.shape
    # A value equal to the threshold counts as exceedance.
    exceed = (tmax >= thresholds[None, :, :]).astype(np.int8)
    number = np.zeros((n_time, n_lat, n_lon), dtype=np.int32)
    episodes = []
    for i in range(n_lat):
        for j in range(n_lon):
            labels, cell_episodes = label_cell_by_season(
                exceed[:, i, j], dates,
                min_len=settings["min_duration_days"],
                min_ones=settings["min_exceeding_days"])
            number[:, i, j] = labels
            for hw_number, start, end, length, n_exceed, n_bridged in cell_episodes:
                episodes.append((i, j, hw_number, str(start), str(end),
                                 length, n_exceed, n_bridged))
    return exceed, number, episodes
