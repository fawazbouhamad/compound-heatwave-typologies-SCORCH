# SPDX-License-Identifier: GPL-3.0-only
"""Stage 2 rules: regional threshold, consecutive-day sequences, event attachment."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import regional


def test_regional_threshold_higher_method_returns_an_observed_count():
    counts = np.array([0] * 97 + [10, 20, 30])    # zero-count days are included
    threshold = regional.regional_threshold(counts, 0.975, "higher")
    # Position 0.975 * 99 = 96.525 lies between the last zero and the 10: "higher" takes the 10.
    assert threshold == 10.0
    assert regional.regional_threshold(counts, 0.975, "linear") == pytest.approx(5.25)
    assert (counts >= threshold).sum() == 3


def test_group_consecutive_days_identifies_sequences():
    dates = ["2001-08-08", "2001-08-07", "2001-08-09", "2001-08-12", "1983-07-13"]
    sequences = regional.group_consecutive_days(dates)
    assert sequences["date"].tolist() == ["1983-07-13", "2001-08-07", "2001-08-08",
                                          "2001-08-09", "2001-08-12"]
    assert sequences["sequence_start_date"].tolist() == ["1983-07-13", "2001-08-07", "2001-08-07",
                                                         "2001-08-07", "2001-08-12"]


def test_group_consecutive_days_handles_no_dates():
    sequences = regional.group_consecutive_days([])
    assert len(sequences) == 0 and list(sequences.columns) == ["date", "sequence_start_date"]


def test_attach_compound_events():
    event_days = pd.DataFrame({"date": ["2001-08-07", "2001-08-08"], "compound_event_id": [10, 10],
                               "day_in_event": [1, 2]})
    records = pd.DataFrame({"date": ["2001-08-08", "2001-08-07"], "value": [1.0, 2.0]})
    attached = regional.attach_compound_events(records, event_days)
    assert attached["compound_event_id"].tolist() == [10, 10]
    with pytest.raises(ValueError, match="already contain"):
        regional.attach_compound_events(attached, event_days)
    with pytest.raises(ValueError, match="missing from the step-5"):
        regional.attach_compound_events(pd.DataFrame({"date": ["1999-01-01"]}), event_days)
    duplicated = pd.concat([event_days, event_days])
    with pytest.raises(ValueError, match="duplicate dates"):
        regional.attach_compound_events(records, duplicated)
