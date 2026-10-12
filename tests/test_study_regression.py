# SPDX-License-Identifier: GPL-3.0-only
"""Regression test against the published study results.

Needs the archived input: set SCORCH_STUDY_INPUT to the path of
tmax_processed.nc (DOI 10.5281/zenodo.21717874). The five stages then run
into a temporary folder and every catalog, Table 1 and supporting CSV must
equal the committed files byte for byte (after line-ending normalization).
Takes about three minutes.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from scorch import config as cfg, io
from scorch.stages import run_stages

REPO = Path(__file__).resolve().parents[1]
STUDY_INPUT = os.environ.get("SCORCH_STUDY_INPUT", "")
EXPECTED_SHA256 = "fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90"
CATALOGS = [
    "catalogs/01_heatwave_events/heatwave_events.csv",
    "catalogs/02_regionally_extensive_heatwaves/selected_days.csv",
    "catalogs/02_regionally_extensive_heatwaves/selected_day_cells.csv",
    "catalogs/03_extreme_heatwave_clusters/clusters.csv",
    "catalogs/03_extreme_heatwave_clusters/cluster_cells.csv",
    "catalogs/04_extreme_heatwave_ellipses/ellipses.csv",
    "catalogs/05_compound_heatwave_events/compound_events.csv",
    "catalogs/05_compound_heatwave_events/event_days.csv",
    "tables/table_1.csv",
]
SUPPORTING = [
    "heatwave_detection/cell_thresholds.csv",
    "regional_selection/daily_summary.csv",
    "clustering/parameter_selection.csv",
    "clustering/daily_settings.csv",
    "clustering/sequence_settings.csv",
    "annual_event_counts.csv",
    "annual_mean_duration.csv",
    "duration_trend_statistics.csv",
    "in_text_statistics.csv",
]
COMPARED_FILES = ([(f"output/{name}", name) for name in CATALOGS] +
                  [(f"supporting_files/additional_results/main_workflow/{name}",
                    f"supporting_results/main_workflow/{name}") for name in SUPPORTING])

pytestmark = [pytest.mark.study, pytest.mark.slow,
              pytest.mark.skipif(not STUDY_INPUT, reason="set SCORCH_STUDY_INPUT to tmax_processed.nc")]


@pytest.fixture(scope="module")
def study_run(tmp_path_factory):
    source = Path(STUDY_INPUT)
    assert source.is_file(), source
    assert io.sha256_of(source) == EXPECTED_SHA256, "not the version 1.0.0 study input"
    out = tmp_path_factory.mktemp("study")
    config = cfg.load_config(overrides={"paths": {
        "base_dir": str(REPO), "input_dir": str(source.parent), "input_file": source.name,
        "output_dir": str(out), "supporting_results_dir": str(out / "supporting_results"),
        "run_records_dir": str(out / "run_records")}}, environment={})
    results = dict(run_stages(config, log=lambda *_: None))
    return out, results


def test_headline_numbers(study_run):
    _, results = study_run
    assert results[1]["local_events"] == 176061
    assert results[2]["selected_days"] == 395 and results[2]["regional_threshold_cells"] == 371
    assert results[3]["clusters"] == 753 and results[4]["ellipses"] == 753
    assert results[5]["events_by_type"] == {"1": 3, "2": 3, "3": 20, "4": 25}


@pytest.mark.parametrize("tracked, produced", COMPARED_FILES)
def test_files_match_the_published_results(study_run, tracked, produced):
    out, _ = study_run
    expected = (REPO / tracked).read_bytes().replace(b"\r\n", b"\n")
    actual = (out / produced).read_bytes().replace(b"\r\n", b"\n")
    assert actual == expected, f"{produced} differs from {tracked}"


def test_flags_file_matches(study_run):
    out, _ = study_run
    tracked = REPO / "supporting_files/additional_results/main_workflow/heatwave_detection/heatwave_flags.nc"
    produced = out / "supporting_results/main_workflow/heatwave_detection/heatwave_flags.nc"
    assert produced.read_bytes() == tracked.read_bytes()
