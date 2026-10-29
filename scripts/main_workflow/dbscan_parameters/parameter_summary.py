"""Calculate the DBSCAN summaries used by Figures A1 and A2.

Reads: Cluster and compound-event catalogs.
Writes: supporting_files/additional_results/main_workflow/dbscan_parameters/
Run: python scripts/main_workflow/dbscan_parameters/parameter_summary.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts" / "main_workflow"))
sys.path.insert(0, str(REPO / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402
from functions import figure_a2_panels, statistics_table

OUTPUT_DIR = REPO / "supporting_files" / "additional_results" / "main_workflow" / "dbscan_parameters"

# Figure A1 uses Event 10, the Type 4 example in Figure 7.
FIGURE_A1_EVENT_ID = 10
def write_figure_a1_tables(selection, daily, event_id, output_dir):
    """Matrices, retained pairs and per-day settings of one event."""
    event_days = daily[daily["compound_event_id"] == event_id]
    days = event_days["date"].tolist()
    eps_values = sorted(selection["eps"].unique())
    min_samples_values = sorted(selection["min_samples"].unique())
    paths, retained_rows, settings_rows = [], [], []
    for day in days:
        rows = selection[selection["date"] == day]
        if len(rows) != len(eps_values) * len(min_samples_values):
            raise SystemExit(f"{day}: incomplete parameter search")
        matrix = rows.pivot(index="min_samples", columns="eps", values="n_clusters")
        matrix = matrix.loc[min_samples_values, eps_values]
        path = output_dir / f"figure_A1_matrix_{day}.csv"
        matrix.to_csv(path, lineterminator="\n")
        paths.append(path)
        retained = rows[rows["retained"] == 1]
        retained_rows.append(retained[["date", "eps", "min_samples", "n_clusters"]])
        settings = event_days[event_days["date"] == day].iloc[0]
        # The mean of the retained pairs is the day's selected setting; check
        # it against the value the clustering step recorded. The tolerance
        # allows for pandas' default number reader, which can change the last
        # binary digit of a value read from CSV (about 4e-16 here).
        if (abs(retained["eps"].mean() - settings["daily_eps"]) > 1e-12
                or abs(retained["min_samples"].mean()
                       - settings["daily_min_samples_mean"]) > 1e-12):
            raise SystemExit(f"{day}: retained pairs do not give the recorded daily setting")
        settings_rows.append({
            "date": day, "compound_event_id": event_id,
            "n_pairs_retained": len(retained),
            "modal_cluster_counts": settings["modal_cluster_counts"],
            "mean_eps_of_retained_pairs": retained["eps"].mean(),
            "mean_min_samples_of_retained_pairs": retained["min_samples"].mean(),
            "daily_min_samples": settings["daily_min_samples"],
            "event_eps": settings["event_eps"],
            "event_min_samples": settings["event_min_samples"],
            "n_clusters_with_event_settings": settings["n_clusters"]})
    pd.concat(retained_rows).to_csv(output_dir / "figure_A1_retained_pairs.csv", index=False,
                                    lineterminator="\n")
    pd.DataFrame(settings_rows).to_csv(output_dir / "figure_A1_daily_settings.csv", index=False,
                                       lineterminator="\n")
    return days, paths + [output_dir / "figure_A1_retained_pairs.csv",
                          output_dir / "figure_A1_daily_settings.csv"]


# Figure A.2


def main():
    started = time.time()
    config = f.load_config()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clustering_dir = f.supporting_results_dir(config) / "main_workflow/clustering"
    events_dir = f.catalog_dir(config, "05_compound_heatwave_events")
    inputs = [clustering_dir / "parameter_selection.csv", clustering_dir / "daily_settings.csv",
              clustering_dir / "sequence_settings.csv", events_dir / "compound_events.csv",
              events_dir / "event_days.csv"]
    # Use pandas' default number reader for the manuscript statistics.
    # Reading with an exact decimal parser changes some values in the 16th digit.
    selection, daily, sequence_settings, events, event_days = (pd.read_csv(p) for p in inputs)
    daily = f.attach_compound_events(daily, event_days)
    daily = daily.rename(columns={"sequence_eps": "event_eps",
                                  "sequence_min_samples": "event_min_samples"})
    event_settings = sequence_settings.merge(
        events[["compound_event_id", "start_date", "end_date"]],
        left_on=["sequence_start_date", "end_date"],
        right_on=["start_date", "end_date"], validate="one_to_one")
    event_settings = event_settings.rename(columns={"sequence_eps": "event_eps",
                                                     "sequence_min_samples": "event_min_samples"})

    # Figure A.1
    days, outputs = write_figure_a1_tables(selection, daily, FIGURE_A1_EVENT_ID, OUTPUT_DIR)
    if days != config["representative_events"]["figure_7_type_4"]:
        raise SystemExit(f"Event {FIGURE_A1_EVENT_ID} days {days} are not the Figure 7 days")
    print(pd.read_csv(OUTPUT_DIR / "figure_A1_daily_settings.csv").to_string(index=False))

    # Figure A.2: settings applied on each day (the event's) and per event.
    type_names = events.set_index("compound_event_id")["event_type"].map(lambda t: f"Type {t}")
    per_day = pd.DataFrame({"date": daily["date"],
                            "compound_event_id": daily["compound_event_id"],
                            "type": daily["compound_event_id"].map(type_names),
                            "min_samples": daily["event_min_samples"],
                            "eps": daily["event_eps"]})
    per_event = pd.DataFrame({"compound_event_id": event_settings["compound_event_id"],
                              "type": event_settings["compound_event_id"].map(type_names),
                              "min_samples": event_settings["event_min_samples"],
                              "eps": event_settings["event_eps"]})
    per_day.to_csv(OUTPUT_DIR / "figure_A2_values_per_day.csv", index=False, lineterminator="\n")
    per_event.to_csv(OUTPUT_DIR / "figure_A2_values_per_event.csv", index=False,
                     lineterminator="\n")
    panels = figure_a2_panels(per_day, per_event)
    statistics = statistics_table(panels)
    statistics.to_csv(OUTPUT_DIR / "figure_A2_statistics.csv", index=False, lineterminator="\n")
    outputs += [OUTPUT_DIR / "figure_A2_values_per_day.csv",
                OUTPUT_DIR / "figure_A2_values_per_event.csv",
                OUTPUT_DIR / "figure_A2_statistics.csv"]

    counts = per_event["type"].value_counts().sort_index().to_dict()
    print(f"Figure A.2: {len(per_day)} days, {len(per_event)} events, events by type {counts}")
    print(statistics[["panel", "group", "n", "mean", "median", "min", "max"]].to_string(index=False))

    f.write_json(f.REPO_DIR / ".scorch/logs/run-records/appendix_A.json", {
        "script": f.relative_to_repo(__file__),
        "finished_utc": f.now_utc(),
        "runtime_seconds": round(time.time() - started, 1),
        "software": f.software_versions(),
        "inputs": {f.relative_to_repo(p): f.sha256_of(p) for p in inputs},
        "outputs": {f.relative_to_repo(p): f.sha256_of(p) for p in outputs}})
    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
