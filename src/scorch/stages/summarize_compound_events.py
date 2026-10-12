# SPDX-License-Identifier: GPL-3.0-only
"""Stage 5: classify compound events and calculate the manuscript summaries.

Reads the preceding catalogs; writes Catalog 5 (``compound_events.csv``,
``event_days.csv``), ``tables/table_1.csv`` and supporting calculations
(annual counts, annual mean durations, trend statistics and, for the study
region, the in-text statistics of Section 4).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import config as cfg, events as ev, io, regional
from . import StageError

STAGE_NAME = "05_summarize_compound_events"


def summarize_events(days, settings, ellipses):
    rows = []
    for event_id, event_days in days.groupby("compound_event_id", sort=True):
        dates = sorted(event_days["date"])
        members = ellipses[ellipses["compound_event_id"] == event_id]
        # Clusters on each day, in date order (a day without clusters is 0 and
        # is rejected by classify_event).
        per_day = members.groupby("date").size().reindex(dates, fill_value=0).to_numpy()
        try:
            event_type = ev.classify_event(per_day)
        except ValueError as error:
            raise StageError(f"compound event {int(event_id)} ({dates[0]} to {dates[-1]}) "
                             f"cannot be classified: {error}") from error
        azimuth, resultant = ev.axial_mean_and_resultant(members["azimuth_north_cw_deg"])
        # Largest ellipse of the event; if two share the largest area the
        # earlier one (by date, then cluster number) is taken.
        largest = members.loc[members["area_km2"].idxmax()]
        s = settings.loc[event_id]
        start, end = pd.Timestamp(dates[0]), pd.Timestamp(dates[-1])
        rows.append({
            "compound_event_id": int(event_id),
            "start_date": dates[0], "end_date": dates[-1],
            "duration_days": (end - start).days + 1,
            "n_selected_days": len(dates),
            "event_type": event_type,
            "event_type_name": ev.TYPE_NAMES[event_type],
            "daily_cluster_counts": "-".join(str(int(n)) for n in per_day),
            "total_clusters": int(per_day.sum()),
            "min_daily_clusters": int(per_day.min()),
            "max_daily_clusters": int(per_day.max()),
            "n_single_cluster_days": int((per_day == 1).sum()),
            "n_multi_cluster_days": int((per_day > 1).sum()),
            "event_eps": s["event_eps"],
            "event_min_samples": int(s["event_min_samples"]),
            "mean_area_km2": float(members["area_km2"].mean()),
            "mean_axis_ratio": float(members["axis_ratio"].mean()),
            "axial_mean_azimuth_deg": azimuth,
            "axial_resultant_length": resultant,
            "largest_ellipse_id": int(largest["ellipse_id"]),
            "largest_ellipse_date": largest["date"],
            "largest_ellipse_area_km2": float(largest["area_km2"]),
        })
    return pd.DataFrame(rows)


def in_text_statistics(events, ellipses, daily, n_local_events, n_cells, boundaries):
    """Numbers quoted in the text of Sections 3.1 and 4 that are not in a
    table or figure, each with its definition.

    The northern Arabian group (Section 4) is the set of ellipses whose
    temperature-weighted centroid lies in Syria or Iraq, or in Saudi Arabia
    at or north of 25 degrees N (Natural Earth 1:10m country polygons).
    Northwest-southeast axes have an azimuth between -67.5 and -22.5 degrees.
    """
    from .. import geography

    e = ellipses.merge(events[["compound_event_id", "event_type"]], on="compound_event_id")
    e["country"] = geography.country_of_points(e["centroid_lon_weighted"],
                                               e["centroid_lat_weighted"], boundaries)
    group = e["country"].isin(["Syria", "Iraq"]) | (
        (e["country"] == "Saudi Arabia") & (e["centroid_lat_weighted"] >= 25.0))
    nwse = e["azimuth_north_cw_deg"].between(-67.5, -22.5)
    group_type4 = group & (e["event_type"] == 4)
    events_with_group_type4 = e.loc[group_type4, "compound_event_id"].unique()
    in_those_events = e["compound_event_id"].isin(events_with_group_type4)
    n_years = daily["date"].str[:4].nunique()

    rows = [
        ("local heatwave events", n_local_events, "count", "Section 3.1"),
        ("local events per grid cell", n_local_events / n_cells, "events per cell", "Section 3.1"),
        ("local events per grid cell per year", n_local_events / n_cells / n_years,
         "events per cell per year", "Section 3.1"),
        ("local events per study day", n_local_events / len(daily),
         "events per day", "Section 3.1 (the text calls this simultaneous cells)"),
        ("mean simultaneous heatwave cells per study day",
         daily["n_heatwave_cells"].mean(), "cells per day", "Section 3.1 check"),
    ]

    def describe(mask, label):
        return [
            (f"{label}: ellipses", int(mask.sum()), "count", "Section 4"),
            (f"{label}: share of all ellipses", mask.mean() * 100, "percent", "Section 4"),
            (f"{label}: mean area", e.loc[mask, "area_km2"].mean() / 1e6,
             "million km2", "Section 4"),
            (f"{label}: mean axis ratio", e.loc[mask, "axis_ratio"].mean(), "L2/L1", "Section 4"),
            (f"{label}: northwest-southeast share", nwse[mask].mean() * 100, "percent",
             "Section 4"),
        ]

    rows += describe(group, "northern Arabian group")
    rows += describe(~group, "all other ellipses")
    rows += describe(group_type4, "northern Arabian group, Type 4 ellipses")
    rows += [("events containing northern Arabian Type 4 ellipses",
              len(events_with_group_type4), "count", "Section 4.2.1"),
             ("ellipses in those events", int(in_those_events.sum()), "count",
              "Section 4.2.1")]
    rows += describe(in_those_events & ~group_type4,
                     "other ellipses of those events")
    return pd.DataFrame(rows, columns=["statistic", "value", "unit", "where_reported"])


def main(config=None, log=print):
    config = config if config is not None else cfg.load_config()
    stage2 = cfg.catalog_dir(config, "02_regionally_extensive_heatwaves")
    stage4 = cfg.catalog_dir(config, "04_extreme_heatwave_ellipses")
    out = cfg.catalog_dir(config, "05_compound_heatwave_events")
    tables = cfg.output_subdir(config, "tables")
    results = cfg.supporting_results_dir(config) / "main_workflow"
    results.mkdir(parents=True, exist_ok=True)
    inputs = [stage2 / "selected_days.csv", results / "clustering" / "sequence_settings.csv",
              stage4 / "ellipses.csv", results / "regional_selection" / "daily_summary.csv",
              results / "heatwave_detection" / "cell_thresholds.csv"]
    for needed in inputs:
        if not needed.is_file():
            raise io.InputError(f"result of an earlier stage not found: {needed}\n"
                                "Run stages 1 to 4 first")
    days = pd.read_csv(inputs[0])
    sequences = regional.group_consecutive_days(days["date"])
    event_days = sequences[["date"]].copy()
    event_days["compound_event_id"] = pd.factorize(sequences["sequence_start_date"], sort=True)[0] + 1
    event_days["day_in_event"] = event_days.groupby("compound_event_id").cumcount() + 1
    days = days.merge(event_days, on="date", validate="one_to_one")
    settings = pd.read_csv(inputs[1], float_precision="round_trip")
    sequence_events = sequences.assign(compound_event_id=event_days["compound_event_id"])
    settings = settings.merge(sequence_events[["sequence_start_date", "compound_event_id"]]
                              .drop_duplicates(), on="sequence_start_date", validate="one_to_one")
    settings = settings.rename(columns={"sequence_eps": "event_eps",
                                        "sequence_min_samples": "event_min_samples"})
    settings = settings.set_index("compound_event_id")
    ellipses = regional.attach_compound_events(
        pd.read_csv(inputs[2], float_precision="round_trip"), event_days)

    daily = pd.read_csv(inputs[3])
    cell_table = pd.read_csv(inputs[4])
    events = summarize_events(days, settings, ellipses)
    if events.empty:
        raise StageError("no compound events were formed from the selected days")
    events["year"] = events["start_date"].str[:4].astype(int)
    table = ev.table_1(events)
    trends, duration_series = ev.trend_statistics(events)

    # Annual event counts by type, every year from the first event to the last study year.
    years = np.arange(events["year"].min(), int(config["study_area"]["last_year"]) + 1)
    counts = (events.groupby(["year", "event_type"]).size().unstack(fill_value=0)
              .reindex(index=years, columns=[1, 2, 3, 4], fill_value=0))
    counts.columns = [f"type_{t}_events" for t in counts.columns]
    counts["all_events"] = counts.sum(axis=1)
    counts.index.name = "year"

    paths = {
        "events": out / "compound_events.csv",
        "event_days": out / "event_days.csv",
        "table": tables / "table_1.csv",
        "counts": results / "annual_event_counts.csv",
        "durations": results / "annual_mean_duration.csv",
        "trends": results / "duration_trend_statistics.csv",
    }
    events.drop(columns="year").to_csv(paths["events"], index=False, lineterminator="\n")
    event_days.to_csv(paths["event_days"], index=False, lineterminator="\n")
    table.to_csv(paths["table"], index=False, lineterminator="\n")
    counts.to_csv(paths["counts"], lineterminator="\n")
    duration_series.to_csv(paths["durations"], index=False, lineterminator="\n")
    trends.to_csv(paths["trends"], index=False, lineterminator="\n")

    # In-text statistics describe the study region (Section 4 of the manuscript).
    boundaries = cfg.country_boundaries_path(config)
    if config["summaries"]["in_text_statistics"]:
        if not boundaries.is_file():
            raise io.InputError(f"country boundaries not found: {boundaries}\n"
                                "Set summaries.in_text_statistics to false for a study "
                                "outside the Eastern Mediterranean and Middle East")
        in_text = in_text_statistics(events, ellipses, daily,
                                     int(cell_table["n_local_events"].sum()),
                                     len(cell_table), boundaries.with_suffix(""))
        paths["in_text"] = results / "in_text_statistics.csv"
        in_text.to_csv(paths["in_text"], index=False, lineterminator="\n")
        inputs.append(boundaries)
    else:
        log("in-text statistics of the study region are switched off (summaries.in_text_statistics)")

    # Consistency checks within this stage.
    if events["n_selected_days"].sum() != len(days) or \
            events["total_clusters"].sum() != len(ellipses):
        raise StageError("event totals do not add up to the selected days and ellipses")
    if (events["duration_days"] != events["n_selected_days"]).any():
        raise StageError("an event has a gap between its selected days")

    type_counts = events["event_type"].value_counts().sort_index().to_dict()
    log(f"compound events: {len(events)}; by type {type_counts}; "
        f"longest {events['duration_days'].max()} days")
    notes = {"events": int(len(events)),
             "events_by_type": {str(k): int(v) for k, v in type_counts.items()}}
    io.write_run_record(config, STAGE_NAME, inputs, list(paths.values()), notes=notes)
    return notes
