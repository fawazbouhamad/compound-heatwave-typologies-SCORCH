# SPDX-License-Identifier: GPL-3.0-only
"""Stage 1 rules: thresholds, exceedance equality, bridging and season boundaries."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import heatwaves
from conftest import SEASON_2000, exceedance_series

SETTINGS = {"min_duration_days": 3, "min_exceeding_days": 3}


def test_cell_thresholds_use_linear_interpolation():
    rng = np.random.default_rng(0)
    tmax = rng.normal(30, 3, size=(50, 2, 3))
    thresholds = heatwaves.cell_thresholds(tmax, 0.95)
    expected = np.percentile(tmax.reshape(50, -1), 95, axis=0, method="linear").reshape(2, 3)
    assert np.allclose(thresholds, expected)


def test_value_equal_to_threshold_counts_as_exceedance():
    dates = np.array([d.strftime("%Y-%m-%d") for d in SEASON_2000])
    tmax = np.full((len(dates), 1, 1), 20.0)
    tmax[10:13, 0, 0] = 30.0                      # exactly the threshold on three days
    thresholds = np.array([[30.0]])
    exceed, number, episodes = heatwaves.detect_all_cells(tmax, thresholds, dates, SETTINGS)
    assert exceed[10:13, 0, 0].tolist() == [1, 1, 1] and exceed.sum() == 3
    assert len(episodes) == 1 and episodes[0][5:] == (3, 3, 0)   # length, exceeding, bridged
    assert number[10:13, 0, 0].tolist() == [1, 1, 1]


def test_isolated_non_exceeding_day_is_bridged():
    series = exceedance_series(SEASON_2000, [10, 11, 13, 14])     # 1 1 0 1 1
    labels, episodes = heatwaves.label_heatwaves(series)
    assert len(episodes) == 1
    number, start, end, length, n_exceeding, n_bridged = episodes[0]
    assert (length, n_exceeding, n_bridged) == (5, 4, 1)
    assert labels[10:15].tolist() == [1, 1, 1, 1, 1] and labels[12] == 1


def test_two_consecutive_non_exceeding_days_end_a_heatwave():
    series = exceedance_series(SEASON_2000, [10, 11, 12, 15, 16, 17])   # 1 1 1 0 0 1 1 1
    labels, episodes = heatwaves.label_heatwaves(series)
    assert len(episodes) == 2
    assert [e[3] for e in episodes] == [3, 3] and labels[13:15].tolist() == [0, 0]


def test_isolated_or_short_exceedances_are_not_heatwaves():
    assert heatwaves.label_heatwaves(exceedance_series(SEASON_2000, [10]))[1] == []
    assert heatwaves.label_heatwaves(exceedance_series(SEASON_2000, [10, 11]))[1] == []
    # Three days long but only two exceeding days: too few exceedances.
    assert heatwaves.label_heatwaves(exceedance_series(SEASON_2000, [10, 12]))[1] == []
    # Four exceeding days spread as 1 0 1 0 1 0 1: one heatwave of 7 days with 3 bridged days.
    _, episodes = heatwaves.label_heatwaves(exceedance_series(SEASON_2000, [10, 12, 14, 16]))
    assert len(episodes) == 1 and episodes[0][3:] == (7, 4, 3)


def test_season_boundary_splits_heatwaves():
    """30 September and 1 April of the next year are not consecutive days."""
    two_seasons = SEASON_2000.append(pd.date_range("2001-04-01", "2001-09-30", freq="D"))
    dates = np.array([d.strftime("%Y-%m-%d") for d in two_seasons])
    exceed = np.zeros(len(dates), dtype=np.int8)
    exceed[180:183] = 1                           # 28-30 September 2000
    exceed[183:186] = 1                           # 1-3 April 2001
    labels, episodes = heatwaves.label_cell_by_season(exceed, dates, 3, 3)
    assert len(episodes) == 2
    assert [str(e[1]) for e in episodes] == ["2000-09-28", "2001-04-01"]
    assert labels[180:186].tolist() == [1, 1, 1, 2, 2, 2]     # numbers continue through the record
    # Two days at the end of one season and one at the start of the next: no heatwave.
    exceed[:] = 0
    exceed[181:184] = 1
    assert heatwaves.label_cell_by_season(exceed, dates, 3, 3)[1] == []


def test_detect_all_cells_peak_and_summary_shapes():
    dates = np.array([d.strftime("%Y-%m-%d") for d in SEASON_2000])
    tmax = np.full((len(dates), 2, 2), 20.0)
    tmax[5:8, 1, 0] = [31.0, 35.0, 35.0]          # one heatwave in one cell
    thresholds = np.full((2, 2), 30.0)
    exceed, number, episodes = heatwaves.detect_all_cells(tmax, thresholds, dates, SETTINGS)
    assert exceed.shape == number.shape == (len(dates), 2, 2)
    assert episodes == [(1, 0, 1, "2000-04-06", "2000-04-08", 3, 3, 0)]


@pytest.mark.parametrize("ones, expected", [
    ([0, 1, 2, 3], 1),                 # one run of four
    ([0, 1, 2, 5, 6, 7], 2),           # two runs separated by two empty days
    ([0, 2, 4], 1),                    # bridged single gaps form one heatwave
])
def test_label_counts(ones, expected):
    assert len(heatwaves.label_heatwaves(exceedance_series(SEASON_2000, ones))[1]) == expected
