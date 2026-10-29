"""Calculate ellipse purity over the scale-factor grid for Figure B1.

Reads: Cluster-cell catalog and study settings in the shared Python functions.
Writes: supporting_files/additional_results/main_workflow/ellipse_sensitivity/
Run: python scripts/main_workflow/ellipse_sensitivity/pca_sensitivity.py
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402

OUTPUT_DIR = REPO / "supporting_files" / "additional_results" / "main_workflow" / "ellipse_sensitivity"


MATRIX_FILES = {"minimum": "purity_matrix_minimum.csv",
                "maximum": "purity_matrix_maximum.csv",
                "mean": "purity_matrix_mean.csv",
                "sum": "purity_matrix_sum.csv"}


def cell_key(lon, lat):
    """Cell centre as a (lon, lat) pair rounded to 6 decimals, for set matching."""
    return (round(float(lon), 6), round(float(lat), 6))


def ellipse_footprint(anchor_lon, anchor_lat, centre_lon, centre_lat,
                      eigvecs, eigvals, scale_major, scale_minor,
                      eigenvalue_floor_km2):
    """Grid-cell centres inside the ellipse, as a set of (lon, lat) pairs.

    The candidate centres form the 1-degree lattice through the anchor cell
    (one of the cluster's own cells). The box is centred on the lattice
    point nearest the ellipse centre and has a half-width of
    ceil(R / km_per_degree + 1) + 1 lattice steps, R being the longer
    semi-axis in km, so it always extends beyond the ellipse.
    """
    a = scale_major * math.sqrt(max(float(eigvals[0]), eigenvalue_floor_km2))
    b = scale_minor * math.sqrt(max(float(eigvals[1]), eigenvalue_floor_km2))
    longest_semi_axis_km = max(a, b)
    km_per_deg_lon = f.KM_PER_DEG_LON_AT_EQUATOR * math.cos(math.radians(centre_lat))
    half_width_lon = int(math.ceil(longest_semi_axis_km / km_per_deg_lon + 1)) + 1
    half_width_lat = int(math.ceil(longest_semi_axis_km / f.KM_PER_DEG_LAT + 1)) + 1
    lon_steps = round(centre_lon - anchor_lon)
    lat_steps = round(centre_lat - anchor_lat)
    candidate_lon = anchor_lon + np.arange(lon_steps - half_width_lon,
                                           lon_steps + half_width_lon + 1)
    candidate_lat = anchor_lat + np.arange(lat_steps - half_width_lat,
                                           lat_steps + half_width_lat + 1)
    grid_lon, grid_lat = np.meshgrid(candidate_lon, candidate_lat)
    grid_lon, grid_lat = grid_lon.ravel(), grid_lat.ravel()
    # The ellipse is centred on (centre_lon, centre_lat), which is also the
    # origin of the km projection, so its centre in km is exactly (0, 0).
    inside = f.points_inside_ellipse(grid_lon, grid_lat, centre_lon, centre_lat,
                                     np.zeros(2), eigvecs, eigvals,
                                     scale_major, eigenvalue_floor_km2,
                                     scale_factor_minor=scale_minor)
    return {cell_key(x, y) for x, y in zip(grid_lon[inside], grid_lat[inside])}


def read_clusters(cells_path):
    """One entry per cluster with its PCA fit, in cluster_id order (1..753).

    cluster_id runs by date and then by DBSCAN label. Members are kept
    west to east, then south to north, as in the main ellipse step.
    The order affects the last digit of the mean position.
    """
    cells = pd.read_csv(cells_path).dropna(subset=["cluster_id"])
    cells["cluster_id"] = cells["cluster_id"].astype(int)
    clusters = []
    for cluster_id, members in cells.groupby("cluster_id", sort=True):
        members = members.iloc[f.longitude_major_order(members["lon"], members["lat"])]
        lon = members["lon"].to_numpy(dtype=float)
        lat = members["lat"].to_numpy(dtype=float)
        centre_lon, centre_lat = float(lon.mean()), float(lat.mean())
        _, eigvecs, eigvals = f.fit_pca(f.project_lonlat_to_km(lon, lat, centre_lon, centre_lat))
        clusters.append({"cluster_id": int(cluster_id),
                         "date": members["date"].iloc[0],
                         "anchor_lon": float(lon[0]), "anchor_lat": float(lat[0]),
                         "centre_lon": centre_lon, "centre_lat": centre_lat,
                         "eigvecs": eigvecs, "eigvals": eigvals,
                         "own_cells": {cell_key(x, y) for x, y in zip(lon, lat)}})
    return clusters


def main():
    started = time.time()
    config = f.load_config()
    scale_factors = [float(s) for s in config["appendix_b"]["scale_factors"]]
    eigenvalue_floor_km2 = float(config["ellipses"]["eigenvalue_floor_km2"])
    cells_path = f.catalog_dir(config, "03_extreme_heatwave_clusters") / "cluster_cells.csv"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    clusters = read_clusters(cells_path)
    n_scales, n_clusters = len(scale_factors), len(clusters)
    print(f"clusters: {n_clusters}; scale factors per axis: {scale_factors}")

    # purity[i, j, k]: major-axis scale i, minor-axis scale j, cluster k.
    purity = np.zeros((n_scales, n_scales, n_clusters), dtype=float)
    long_rows = []
    for i, scale_major in enumerate(scale_factors):
        for j, scale_minor in enumerate(scale_factors):
            for k, c in enumerate(clusters):
                footprint = ellipse_footprint(c["anchor_lon"], c["anchor_lat"],
                                              c["centre_lon"], c["centre_lat"],
                                              c["eigvecs"], c["eigvals"],
                                              scale_major, scale_minor,
                                              eigenvalue_floor_km2)
                n_intersection = len(c["own_cells"] & footprint)
                value = 100.0 * n_intersection / len(footprint) if footprint else float("nan")
                purity[i, j, k] = value
                long_rows.append((c["cluster_id"], c["date"], scale_major, scale_minor,
                                  len(c["own_cells"]), len(footprint), n_intersection, value))

    # Summaries over clusters. A small ellipse can contain no cell centre at
    # all; its purity is undefined (NaN), so it is left out of the minimum,
    # maximum and mean for that combination.
    minimum = np.nanmin(purity, axis=2)
    maximum = np.nanmax(purity, axis=2)
    mean = np.nanmean(purity, axis=2)
    total = minimum + maximum + mean
    n_empty = np.isnan(purity).sum(axis=2)

    index = pd.Index(scale_factors, name="scale_factor_major_axis")
    columns = pd.Index(scale_factors, name="scale_factor_minor_axis")
    for name, values in (("minimum", minimum), ("maximum", maximum),
                         ("mean", mean), ("sum", total)):
        pd.DataFrame(values, index, columns).to_csv(OUTPUT_DIR / MATRIX_FILES[name],
                                                    lineterminator="\n")
    by_cluster = pd.DataFrame(long_rows, columns=[
        "cluster_id", "date", "scale_major", "scale_minor", "n_cluster_cells",
        "n_footprint_cells", "n_intersection", "purity"])
    by_cluster = by_cluster.sort_values(["cluster_id", "scale_major", "scale_minor"],
                                        kind="stable")
    by_cluster.to_csv(OUTPUT_DIR / "purity_by_cluster.csv", index=False, lineterminator="\n")

    diagonal = pd.DataFrame({
        "scale_factor": scale_factors,
        "minimum_purity": np.diag(minimum), "maximum_purity": np.diag(maximum),
        "mean_purity": np.diag(mean), "sum": np.diag(total),
        "n_clusters_used": n_clusters - np.diag(n_empty)})
    diagonal["largest_sum_on_diagonal"] = diagonal["sum"] == diagonal["sum"].max()
    diagonal.to_csv(OUTPUT_DIR / "diagonal_summary.csv", index=False, lineterminator="\n")


    best = diagonal.loc[diagonal["sum"].idxmax()]
    print(diagonal.to_string(index=False))
    print(f"largest sum on the diagonal: scale factor {best['scale_factor']:g} "
          f"(sum {best['sum']:.6f})")
    outputs = ([OUTPUT_DIR / name for name in MATRIX_FILES.values()]
               + [OUTPUT_DIR / "purity_by_cluster.csv",
                  OUTPUT_DIR / "diagonal_summary.csv"])
    # Record of what was read and written (checksums, software versions).
    f.write_json(f.REPO_DIR / ".scorch/logs/run-records/appendix_B.json", {
        "script": f.relative_to_repo(__file__),
        "finished_utc": f.now_utc(),
        "runtime_seconds": round(time.time() - started, 1),
        "software": f.software_versions(),
        "inputs": {f.relative_to_repo(cells_path): f.sha256_of(cells_path)},
        "outputs": {f.relative_to_repo(p): f.sha256_of(p) for p in outputs},
        "clusters": n_clusters,
        "best_diagonal_scale_factor": float(best["scale_factor"])})
    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
