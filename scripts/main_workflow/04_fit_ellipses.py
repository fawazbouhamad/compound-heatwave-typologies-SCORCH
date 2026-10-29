"""Fit a PCA ellipse and temperature-weighted centroid to each cluster.

Reads Catalog 3 and input temperatures; writes Catalog 4.
Run: python scripts/main_workflow/04_fit_ellipses.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions as f  # noqa: E402


def main():
    config = f.load_config()
    settings = config["ellipses"]
    stage3 = f.catalog_dir(config, "03_extreme_heatwave_clusters")
    out = f.catalog_dir(config, "04_extreme_heatwave_ellipses")
    clusters_path = stage3 / "clusters.csv"
    cells_path = stage3 / "cluster_cells.csv"
    clusters = pd.read_csv(clusters_path)
    cells = pd.read_csv(cells_path).dropna(subset=["cluster_id"])
    cells["cluster_id"] = cells["cluster_id"].astype(int)

    data = f.read_tmax_input(config)
    date_index = {d: k for k, d in enumerate(data["dates"])}
    lat_index = {v: k for k, v in enumerate(data["lat"])}
    lon_index = {v: k for k, v in enumerate(data["lon"])}

    # Member-cell Tmax at the precision used for the manuscript centroids.
    t = cells["date"].map(date_index).to_numpy()
    i = cells["lat"].map(lat_index).to_numpy()
    j = cells["lon"].map(lon_index).to_numpy()
    member_tmax = data["tmax"][t, i, j]
    if settings["weighted_centroid_tmax_precision"] == "float32":
        member_tmax = f.centroid_temperature_weights(member_tmax).astype(np.float64)
    cells["tmax_c"] = member_tmax

    rows = []
    for cluster in clusters.itertuples():
        # Members in west-to-east, south-to-north order.
        members = cells[cells["cluster_id"] == cluster.cluster_id]
        members = members.iloc[f.longitude_major_order(members["lon"], members["lat"])]
        geometry = f.ellipse_geometry(members["lon"].to_numpy(), members["lat"].to_numpy(),
                                      settings["scale_factor"],
                                      settings["eigenvalue_floor_km2"])
        lon_w, lat_w, weight_sum = f.weighted_centroid(
            members["lon"].to_numpy(), members["lat"].to_numpy(),
            members["tmax_c"].to_numpy())
        row = {"ellipse_id": cluster.cluster_id, "cluster_id": cluster.cluster_id,
               "date": cluster.date,
               "cluster_number_in_day": cluster.cluster_number_in_day}
        row.update({k: v for k, v in geometry.items() if not k.startswith("_")})
        row.update({"centroid_lon_weighted": lon_w, "centroid_lat_weighted": lat_w,
                    "tmax_weight_sum_c": weight_sum,
                    "tmax_min_c": float(members["tmax_c"].min()),
                    "tmax_max_c": float(members["tmax_c"].max()),
                    "scale_factor": settings["scale_factor"],
                    "eigenvalue_floor_km2": settings["eigenvalue_floor_km2"]})
        rows.append(row)
    ellipses = pd.DataFrame(rows)

    path = out / "ellipses.csv"
    ellipses.to_csv(path, index=False, lineterminator="\n")

    # Consistency checks within this stage: one ellipse per cluster and
    # positive geometry.
    if ellipses["cluster_id"].tolist() != clusters["cluster_id"].tolist():
        raise SystemExit("ellipses and clusters do not match one to one")
    if (ellipses[["major_axis_km", "minor_axis_km", "area_km2"]] <= 0).any().any():
        raise SystemExit("non-positive ellipse geometry")
    if ((ellipses["axis_ratio"] <= 0) | (ellipses["axis_ratio"] > 1)).any():
        raise SystemExit("axis ratio outside (0, 1]")

    print(f"ellipses: {len(ellipses)}; area {ellipses['area_km2'].min():,.0f}"
          f" to {ellipses['area_km2'].max():,.0f} km^2")
    f.write_run_record(config, __file__, [clusters_path, cells_path, data["path"]], [path],
                       notes={"ellipses": int(len(ellipses)),
                              "weighted_centroid_tmax_precision":
                                  settings["weighted_centroid_tmax_precision"]})


if __name__ == "__main__":
    main()
