# SPDX-License-Identifier: GPL-3.0-only
"""End-to-end test: the synthetic demo finds exactly the designed events."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scorch import heatwaves
from scorch.demo import run_demo, synthetic

REPO = Path(__file__).resolve().parents[1]
HAS_MATPLOTLIB = importlib.util.find_spec("matplotlib") is not None


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    directory = tmp_path_factory.mktemp("scorch_demo")
    return run_demo(directory, make_figure=HAS_MATPLOTLIB, log=lambda *_: None)


def test_designed_events_are_recovered(demo):
    comparison = demo["comparison"]
    assert comparison["matches_design"].all(), comparison.to_string()
    assert len(demo["extra_events"]) == 0
    assert demo["events_by_type"] == {"1": 1, "2": 1, "3": 2, "4": 2}


def test_demo_matches_the_recorded_expectations(demo):
    expected = json.loads((REPO / "demo" / "expected_results.json").read_text(encoding="utf-8"))
    events = pd.read_csv(demo["directory"] / "output" / "catalogs" / "05_compound_heatwave_events"
                         / "compound_events.csv")
    found = [{"start_date": r.start_date, "end_date": r.end_date, "event_type": int(r.event_type),
              "daily_cluster_counts": str(r.daily_cluster_counts)} for r in events.itertuples()]
    assert found == expected["events"]
    selected = pd.read_csv(demo["directory"] / "output" / "catalogs"
                           / "02_regionally_extensive_heatwaves" / "selected_days.csv")
    assert len(selected) == expected["selected_days"]
    assert expected["seed"] == synthetic.DEMO_SEED


def test_demo_writes_the_promised_files(demo):
    directory = demo["directory"]
    for relative in ["input/tmax_processed.nc", "input/maps/grid_cells.csv", "settings.yml",
                     "output/tables/table_1.csv", "summary.md",
                     "output/run_records/05_summarize_compound_events.json"]:
        assert (directory / relative).is_file(), relative
    table = pd.read_csv(directory / "output" / "tables" / "table_1.csv")
    assert table["number_of_compound_events"].tolist() == [1, 1, 2, 2]
    summary = (directory / "summary.md").read_text(encoding="utf-8")
    assert "invented" in summary and "| yes |" in summary and "| NO |" not in summary


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib is not installed")
def test_demo_figure_is_drawn(demo):
    assert demo["figure"] is not None and demo["figure"].stat().st_size > 10_000


def test_background_weather_alone_stays_below_the_designed_events():
    """Without the designed events the synthetic noise never reaches a blob's size."""
    dates, lat, lon, tmax = synthetic.make_synthetic_tmax(include_events=False)
    labels = np.array([d.strftime("%Y-%m-%d") for d in dates])
    thresholds = heatwaves.cell_thresholds(tmax, 0.95)
    _, number, _ = heatwaves.detect_all_cells(tmax, thresholds, labels,
                                              {"min_duration_days": 3, "min_exceeding_days": 3})
    counts = (number > 0).reshape(len(labels), -1).sum(axis=1)
    smallest_designed_day = 49                                # one radius-4 blob
    assert counts.max() < smallest_designed_day
