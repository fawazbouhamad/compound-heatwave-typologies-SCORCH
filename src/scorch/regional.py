# SPDX-License-Identifier: GPL-3.0-only
"""Stage 2: regionally extensive heatwave days and consecutive-day sequences.

A day is selected when its number of heatwave cells reaches the regional
threshold, a high percentile of the daily heatwave-cell counts over all
study dates (days with zero heatwave cells included). Consecutive selected
dates form one sequence; in stage 5 every sequence becomes one compound
event.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def regional_threshold(daily_heatwave_cells, quantile, method):
    """Regional threshold from heatwave-cell counts on all study dates.

    The "higher" quantile method returns the observed count at or just above
    the quantile position.
    """
    counts = pd.Series(daily_heatwave_cells)
    return float(counts.quantile(quantile, interpolation=method))


def group_consecutive_days(selected_dates):
    """Identify consecutive-date sequences by their start date.

    Step 3 pools DBSCAN settings within each sequence. Compound event IDs
    are assigned only in step 5, after clustering and ellipse fitting.
    """
    dates = pd.DatetimeIndex(sorted(pd.to_datetime(list(selected_dates))))
    if len(dates) == 0:
        return pd.DataFrame({"date": pd.Series(dtype=str),
                             "sequence_start_date": pd.Series(dtype=str)})
    gaps = np.diff(dates.values).astype("timedelta64[D]").astype(int)
    new_sequence = np.concatenate([[True], gaps > 1])
    starts = np.maximum.accumulate(np.where(new_sequence, np.arange(len(dates)), 0))
    return pd.DataFrame({"date": dates.strftime("%Y-%m-%d"),
                         "sequence_start_date": dates[starts].strftime("%Y-%m-%d")})


def attach_compound_events(records, event_days):
    """Link dated records to the completed step-5 event-day catalog in memory."""
    if "compound_event_id" in records:
        raise ValueError("records already contain compound event IDs")
    lookup = event_days.set_index("date")["compound_event_id"]
    if not lookup.index.is_unique:
        raise ValueError("the event-day catalog contains duplicate dates")
    result = records.copy()
    result["compound_event_id"] = result["date"].map(lookup)
    if result["compound_event_id"].isna().any():
        raise ValueError("a record date is missing from the step-5 event-day catalog")
    result["compound_event_id"] = result["compound_event_id"].astype(int)
    return result
