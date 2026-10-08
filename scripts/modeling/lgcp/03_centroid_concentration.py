"""Measure centroid concentration within the fitted intensity zones.

Reads: supporting_files/additional_results/modeling/lgcp/
Writes: supporting_files/additional_results/modeling/lgcp/
Run: python scripts/modeling/lgcp/03_centroid_concentration.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import functions as f  # noqa: E402
import zone_scoring as zs  # noqa: E402

OUTPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp"


def main():
    started = time.time()
    config = f.load_config()
    zone_fractions = config["appendix_c"]["zone_area_fractions"]

    cells = pd.read_csv(OUTPUT_DIR / "predicted_intensity_1deg.csv")
    cells["area_km2"] = zs.cell_area_km2(cells["lat"])
    # Convert the stored intensity per million km^2 back to intensity per km^2.
    intensity = cells["mean_intensity_per_million_km2"].to_numpy(float) / 1e6
    if len(np.unique(intensity)) != len(intensity):
        raise ValueError("tied intensities: the zones would depend on the sort order")

    # Use pandas' default number reader for the centroid coordinates.
    # The manuscript distances use these parsed values, which can differ
    # from the stored values in the last binary digit.
    centroids = pd.read_csv(OUTPUT_DIR / "centroids.csv")
    row_of_cell = pd.Series(np.arange(len(cells)), index=cells["cell_id"])
    centroids["cell_row"] = row_of_cell.reindex(centroids["nearest_cell_id"]).to_numpy()

    per_centroid, zones = zs.score_centroids(intensity, cells, centroids, zone_fractions)
    per_centroid.insert(1, "compound_event_id", centroids["compound_event_id"].to_numpy())
    per_centroid.insert(2, "date", centroids["date"].to_numpy())
    per_centroid.insert(3, "centroid_lon", centroids["centroid_lon"].to_numpy())
    per_centroid.insert(4, "centroid_lat", centroids["centroid_lat"].to_numpy())
    per_centroid.insert(5, "is_daily_largest", centroids["is_daily_largest"].to_numpy())
    per_centroid.insert(6, "is_event_largest", centroids["is_event_largest"].to_numpy())
    groups = zs.summarize_groups(per_centroid, zone_fractions)

    summary_rows = []
    for fraction in zone_fractions:
        label = zs.zone_label(fraction)
        for _, group in groups.iterrows():
            n_in = int(group[f"n_in_{label}_zone"])
            summary_rows.append({
                "zone_area_fraction": fraction,
                "group": group["group"],
                "n_centroids_in_group": int(group["n_centroids"]),
                "n_in_zone": n_in,
                "percent_in_zone": 100.0 * n_in / group["n_centroids"],
                "mean_distance_to_zone_km": group[f"mean_distance_to_{label}_zone_km"],
                "mean_area_weighted_rank": group["mean_area_weighted_rank"],
            })
    summary = pd.DataFrame(summary_rows)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_DIR / "centroid_concentration_summary.csv", index=False, lineterminator="\n")
    per_centroid.to_csv(OUTPUT_DIR / "centroid_concentration_per_centroid.csv", index=False, lineterminator="\n")
    groups.to_csv(OUTPUT_DIR / "centroid_concentration_group_summary.csv", index=False, lineterminator="\n")
    zones.to_csv(OUTPUT_DIR / "concentration_zones.csv", index=False, lineterminator="\n")
    print(summary.to_string(index=False))

    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
