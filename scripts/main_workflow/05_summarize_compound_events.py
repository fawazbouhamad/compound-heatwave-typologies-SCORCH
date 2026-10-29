"""Classify compound events and calculate the manuscript summaries.

Reads the preceding catalogs; writes Catalog 5, output/tables/table_1.csv,
and supporting calculations in supporting_files/additional_results/main_workflow/.
Run: python scripts/main_workflow/05_summarize_compound_events.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions as f  # noqa: E402


def summarize_events(days, settings, ellipses):
    rows = []
    for event_id, event_days in days.groupby("compound_event_id", sort=True):
        dates = sorted(event_days["date"])
        members = ellipses[ellipses["compound_event_id"] == event_id]
        # Clusters on each day, in date order.
        per_day = members.groupby("date").size().reindex(dates).to_numpy()
        event_type = f.classify_event(per_day)
        azimuth, resultant = f.axial_mean_and_resultant(members["azimuth_north_cw_deg"])
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
            "event_type_name": f.TYPE_NAMES[event_type],
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


def table_1(events):
    """Compound heatwave statistics by type (manuscript Table 1)."""
    rows = []
    for t in (1, 2, 3, 4):
        sub = events[events["event_type"] == t]
        rows.append({
            "compound_heatwave_type": f"Type {t}: {f.TYPE_NAMES[t]}",
            "selected_event_days": int(sub["n_selected_days"].sum()),
            "number_of_compound_events": len(sub),
            "longest_event_days": int(sub["duration_days"].max()),
            "mean_duration_days_per_event": round(float(sub["duration_days"].mean()), 2),
            "mean_ellipses_per_event": round(float(sub["total_clusters"].mean()), 2),
            "minimum_ellipses_per_event": int(sub["total_clusters"].min()),
            "maximum_ellipses_per_event": int(sub["total_clusters"].max()),
        })
    return pd.DataFrame(rows)


def trend_statistics(events):
    """Trend tests on the annual mean duration of Type 3 and Type 4 events.

    The year of an event is the year of its first day. Only years with at
    least one event of the type enter the series. Annual event counts are
    reported descriptively and are not trend-tested.
    """
    rows, series = [], []
    for t in (3, 4):
        sub = (events[events["event_type"] == t].groupby("year")["duration_days"]
               .mean().reset_index().sort_values("year"))
        x = sub["year"].to_numpy(dtype=float)
        y = sub["duration_days"].to_numpy(dtype=float)
        ols = linregress(x, y)
        fitted = ols.intercept + ols.slope * x
        residual = y - fitted
        s, z, p, trend = f.mann_kendall_test(y)
        sen = f.sens_slope(x, y)
        rows.append({
            "event_type": t, "event_type_name": f.TYPE_NAMES[t], "n_years": len(sub),
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


def in_text_statistics(events, ellipses, daily, n_local_events, n_cells, boundaries):
    """Numbers quoted in the text of Sections 3.1 and 4 that are not in a
    table or figure, each with its definition.

    The northern Arabian group (Section 4) is the set of ellipses whose
    temperature-weighted centroid lies in Syria or Iraq, or in Saudi Arabia
    at or north of 25 degrees N (Natural Earth 1:10m country polygons).
    Northwest-southeast axes have an azimuth between -67.5 and -22.5 degrees.
    """
    e = ellipses.merge(events[["compound_event_id", "event_type"]], on="compound_event_id")
    e["country"] = f.country_of_points(e["centroid_lon_weighted"],
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


def main():
    config = f.load_config()
    stage2 = f.catalog_dir(config, "02_regionally_extensive_heatwaves")
    stage4 = f.catalog_dir(config, "04_extreme_heatwave_ellipses")
    out = f.catalog_dir(config, "05_compound_heatwave_events")
    tables = f.output_subdir(config, "tables")
    results = f.supporting_results_dir(config) / "main_workflow"
    results.mkdir(parents=True, exist_ok=True)
    boundaries = f.supporting_input_dir(config) / "maps" / "countries_10m" / "ne_10m_admin_0_countries.shp"
    inputs = [stage2 / "selected_days.csv", results / "clustering/sequence_settings.csv",
              stage4 / "ellipses.csv", results / "regional_selection/daily_summary.csv",
              results / "heatwave_detection/cell_thresholds.csv", boundaries]
    days = pd.read_csv(inputs[0])
    sequences = f.group_consecutive_days(days["date"])
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
    ellipses = f.attach_compound_events(pd.read_csv(inputs[2], float_precision="round_trip"),
                                       event_days)

    daily = pd.read_csv(inputs[3])
    cell_table = pd.read_csv(inputs[4])
    events = summarize_events(days, settings, ellipses)
    events["year"] = events["start_date"].str[:4].astype(int)
    table = table_1(events)
    trends, duration_series = trend_statistics(events)
    in_text = in_text_statistics(events, ellipses, daily,
                                 int(cell_table["n_local_events"].sum()),
                                 len(cell_table), boundaries.with_suffix(""))

    # Annual event counts by type, every year from the first event to 2025.
    years = np.arange(events["year"].min(), config["study_area"]["last_year"] + 1)
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
        "in_text": results / "in_text_statistics.csv",
    }
    events.drop(columns="year").to_csv(paths["events"], index=False, lineterminator="\n")
    event_days.to_csv(paths["event_days"], index=False, lineterminator="\n")
    table.to_csv(paths["table"], index=False, lineterminator="\n")
    counts.to_csv(paths["counts"], lineterminator="\n")
    duration_series.to_csv(paths["durations"], index=False, lineterminator="\n")
    trends.to_csv(paths["trends"], index=False, lineterminator="\n")
    in_text.to_csv(paths["in_text"], index=False, lineterminator="\n")

    # Consistency checks within this stage.
    if events["n_selected_days"].sum() != len(days) or \
            events["total_clusters"].sum() != len(ellipses):
        raise SystemExit("event totals do not add up to the selected days and ellipses")
    if (events["duration_days"] != events["n_selected_days"]).any():
        raise SystemExit("an event has a gap between its selected days")

    type_counts = events["event_type"].value_counts().sort_index().to_dict()
    print(f"compound events: {len(events)}; by type {type_counts}; "
          f"longest {events['duration_days'].max()} days")
    f.write_run_record(config, __file__, inputs, list(paths.values()),
                       notes={"events_by_type": {str(k): int(v) for k, v in type_counts.items()}})


if __name__ == "__main__":
    main()
