# SPDX-License-Identifier: GPL-3.0-only
"""Stage 5 rules: event types, Table 1 with a missing type, axial mean, trends."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import events as ev


@pytest.mark.parametrize("counts, expected", [
    ([1], 1), ([3], 2), ([1, 1], 3), ([2, 2, 3], 3), ([1, 2], 4), ([2, 1, 1], 4),
    (np.array([1, 1, 1, 1]), 3), ((5, 1), 4),
])
def test_classify_event(counts, expected):
    assert ev.classify_event(counts) == expected


@pytest.mark.parametrize("counts", [[], [0], [1, 0], [0, 2], np.array([], dtype=int)])
def test_empty_events_and_zero_cluster_days_are_rejected(counts):
    with pytest.raises(ValueError, match="at least one"):
        ev.classify_event(counts)


def events_frame(types):
    rows = []
    for k, t in enumerate(types, start=1):
        days = 1 if t in (1, 2) else 3
        clusters = days * (1 if t in (1, 3) else 2)
        rows.append({"compound_event_id": k, "event_type": t, "n_selected_days": days,
                     "duration_days": days, "total_clusters": clusters,
                     "start_date": f"{1980 + k}-07-01", "year": 1980 + k})
    return pd.DataFrame(rows)


def test_table_1_with_every_type():
    table = ev.table_1(events_frame([1, 2, 3, 3, 4]))
    assert table["number_of_compound_events"].tolist() == [1, 1, 2, 1]
    assert table["selected_event_days"].tolist() == [1, 1, 6, 3]
    assert table["longest_event_days"].tolist() == [1, 1, 3, 3]
    assert table["mean_ellipses_per_event"].tolist() == [1.0, 2.0, 3.0, 6.0]
    assert table["longest_event_days"].dtype.kind == "i"
    assert list(table.columns) == [
        "compound_heatwave_type", "selected_event_days", "number_of_compound_events",
        "longest_event_days", "mean_duration_days_per_event", "mean_ellipses_per_event",
        "minimum_ellipses_per_event", "maximum_ellipses_per_event"]


def test_table_1_handles_a_type_without_events(tmp_path):
    table = ev.table_1(events_frame([1, 3, 4, 4]))
    row = table.set_index("compound_heatwave_type").loc["Type 2: Spatially Clustered"]
    assert row["number_of_compound_events"] == 0 and row["selected_event_days"] == 0
    assert pd.isna(row["longest_event_days"]) and pd.isna(row["mean_ellipses_per_event"])
    assert table["number_of_compound_events"].tolist() == [1, 0, 1, 2]
    path = tmp_path / "table_1.csv"
    table.to_csv(path, index=False, lineterminator="\n")      # must not raise
    assert pd.read_csv(path)["number_of_compound_events"].tolist() == [1, 0, 1, 2]


def test_axial_mean_treats_opposite_axes_as_equal():
    mean, resultant = ev.axial_mean_and_resultant([45.0, 45.0])
    assert mean == pytest.approx(45.0) and resultant == pytest.approx(1.0)
    mean, resultant = ev.axial_mean_and_resultant([-80.0, 80.0])   # both nearly north-south
    assert mean == pytest.approx(90.0) and 0 < resultant < 1
    with pytest.raises(ValueError, match="undefined"):
        ev.axial_mean_and_resultant([0.0, 90.0])
    with pytest.raises(ValueError, match="at least one"):
        ev.axial_mean_and_resultant([])


def test_trend_statistics_need_three_years():
    frame = events_frame([3, 3, 4])                       # two Type-3 years, one Type-4 year
    trends, series = ev.trend_statistics(frame)
    assert trends["n_years"].tolist() == [2, 1]
    assert trends["mann_kendall_result"].str.startswith("not computed").all()
    assert trends["ols_slope_days_per_year"].isna().all()
    assert len(series) == 3
    frame = events_frame([3, 3, 3, 3, 4, 4, 4])
    frame.loc[frame["event_type"] == 3, "duration_days"] = [2, 3, 4, 5]
    trends, _ = ev.trend_statistics(frame)
    type3 = trends.set_index("event_type").loc[3]
    assert type3["n_years"] == 4 and type3["ols_slope_days_per_year"] == pytest.approx(1.0)
    assert type3["sen_slope_days_per_year"] == pytest.approx(1.0)
    assert np.isfinite(type3["mann_kendall_p_value"])


def test_mann_kendall_and_sens_slope():
    s, z, p, trend = ev.mann_kendall_test([1, 2, 3, 4, 5, 6, 7, 8])
    assert s == 28 and p < 0.05 and trend == "increasing"
    assert ev.sens_slope([1, 2, 3, 4], [2, 4, 6, 8]) == 2.0
