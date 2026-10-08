"""Choose DBSCAN settings and group the heatwave cells on selected days.

Reads Catalog 2; writes clusters and cell membership in Catalog 3.
Parameter searches and consecutive-date sequence settings are supporting results.
Run: python scripts/main_workflow/03_cluster_heatwaves.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions as f  # noqa: E402


def main():
    config = f.load_config()
    settings = config["clustering"]
    grid_step = config["study_area"]["grid_step_deg"]
    stage2 = f.catalog_dir(config, "02_regionally_extensive_heatwaves")
    out = f.catalog_dir(config, "03_extreme_heatwave_clusters")
    support = f.supporting_results_dir(config) / "main_workflow/clustering"
    support.mkdir(parents=True, exist_ok=True)
    cells_path = stage2 / "selected_day_cells.csv"
    days_path = stage2 / "selected_days.csv"
    cells = pd.read_csv(cells_path)
    days = pd.read_csv(days_path)
    days = days.merge(f.group_consecutive_days(days["date"]), on="date",
                      validate="one_to_one")

    # Each day's heatwave cells in west-to-east order.
    day_cells = {}
    for date, group in cells.groupby("date", sort=True):
        group = group.iloc[f.longitude_major_order(group["lon"], group["lat"])]
        day_cells[date] = group.reset_index(drop=True)

    # 1-2. Parameter search and daily settings.
    search_rows, daily_rows, distances = [], [], {}
    for day in days.itertuples():
        group = day_cells[day.date]
        d = f.grid_distance_matrix(group["lon"], group["lat"], grid_step)
        distances[day.date] = d
        rows = f.parameter_search(d, settings["eps_values"], settings["min_samples_values"])
        choice = f.select_daily_parameters(rows)
        for r in rows:
            search_rows.append({
                "date": day.date, "sequence_start_date": day.sequence_start_date,
                **r, "retained": int((r["eps"], r["min_samples"]) in choice["retained_pairs"])})
        daily_rows.append({
            "date": day.date, "sequence_start_date": day.sequence_start_date,
            "n_heatwave_cells": len(group),
            "modal_cluster_counts": "|".join(map(str, choice["modal_cluster_counts"])),
            "modal_frequency": choice["modal_frequency"],
            "n_pairs_retained": choice["n_pairs_retained"],
            "daily_eps": choice["daily_eps"],
            "daily_min_samples_mean": choice["daily_min_samples_mean"],
            "daily_min_samples": choice["daily_min_samples"]})
    search = pd.DataFrame(search_rows)
    daily = pd.DataFrame(daily_rows)

    # 3. Pool settings over consecutive selected dates, as in the manuscript.
    # This computational grouping does not assign final compound event IDs.
    sequence = (daily.groupby("sequence_start_date")
             .agg(end_date=("date", "max"),
                  n_days=("date", "size"), sequence_eps=("daily_eps", "max"),
                  sequence_min_samples=("daily_min_samples", "max"))
             .reset_index())
    daily = daily.merge(sequence[["sequence_start_date", "sequence_eps", "sequence_min_samples"]],
                        on="sequence_start_date")

    # 4. Final clustering of each day with its consecutive-sequence settings.
    cluster_rows, cell_rows = [], []
    next_cluster_id = 1
    for day in daily.itertuples():
        group = day_cells[day.date]
        labels, is_core = f.run_dbscan(distances[day.date], day.sequence_eps,
                                       day.sequence_min_samples)
        n_within = (distances[day.date] <= day.sequence_eps).sum(axis=1)
        # Clusters are numbered 0, 1, 2, ... within the day in the order of
        # their DBSCAN labels.
        day_labels = sorted(set(labels.tolist()) - {f.NOISE})
        global_id = {}
        for number, label in enumerate(day_labels):
            global_id[label] = next_cluster_id
            cluster_rows.append({
                "cluster_id": next_cluster_id, "date": day.date,
                "cluster_number_in_day": number, "dbscan_label": label,
                "n_member_cells": int((labels == label).sum()),
                "n_core_cells": int(((labels == label) & is_core).sum()),
                "eps_applied": day.sequence_eps,
                "min_samples_applied": day.sequence_min_samples,
                "n_clusters_on_day": len(day_labels),
                "n_noise_cells_on_day": int((labels == f.NOISE).sum())})
            next_cluster_id += 1
        cell_rows.append(pd.DataFrame({
            "date": day.date,
            "cell_id": group["cell_id"], "lat": group["lat"], "lon": group["lon"],
            "local_event_id": group["local_event_id"],
            "dbscan_label": labels,
            "cluster_id": pd.array([global_id.get(x) for x in labels], dtype="Int64"),
            "is_noise": (labels == f.NOISE).astype(int),
            "is_core_cell": is_core.astype(int),
            "n_cells_within_eps": n_within}))
    clusters = pd.DataFrame(cluster_rows)
    cluster_cells = pd.concat(cell_rows, ignore_index=True)
    daily["n_clusters"] = daily["date"].map(clusters.groupby("date").size()).fillna(0).astype(int)

    paths = {name: support / f"{name}.csv" for name in
             ("parameter_selection", "daily_settings", "sequence_settings")}
    paths.update({name: out / f"{name}.csv" for name in ("clusters", "cluster_cells")})
    search.to_csv(paths["parameter_selection"], index=False, lineterminator="\n")
    daily.to_csv(paths["daily_settings"], index=False, lineterminator="\n")
    sequence.to_csv(paths["sequence_settings"], index=False, lineterminator="\n")
    clusters.to_csv(paths["clusters"], index=False, lineterminator="\n")
    cluster_cells.to_csv(paths["cluster_cells"], index=False, lineterminator="\n")

    # Consistency checks within this stage.
    if len(search) != 63 * len(days):
        raise SystemExit("parameter search is incomplete")
    members = cluster_cells.dropna(subset=["cluster_id"]).groupby("cluster_id").size()
    if not (members.reindex(clusters["cluster_id"]).to_numpy()
            == clusters["n_member_cells"].to_numpy()).all():
        raise SystemExit("cluster sizes and cell memberships disagree")
    if (daily["n_clusters"] == 0).any():
        raise SystemExit("a selected day has no cluster")

    print(f"parameter-search rows: {len(search):,}")
    print(f"clusters: {len(clusters)} on {daily['date'].nunique()} days "
          f"({daily['n_clusters'].min()}-{daily['n_clusters'].max()} per day)")
    f.write_run_record(config, __file__, [cells_path, days_path], list(paths.values()),
                       notes={"clusters": int(len(clusters)),
                              "parameter_search_rows": int(len(search))})


if __name__ == "__main__":
    main()
