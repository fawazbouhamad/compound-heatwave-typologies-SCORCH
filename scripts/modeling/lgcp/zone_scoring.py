"""Shared calculations for centroid concentration and validation zones.

Reads: Called by LGCP and Appendix C scripts.
Writes: Returns calculations to the calling script.
Run: Imported by the calculation scripts.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KM_PER_DEGREE_LATITUDE = 110.574
KM_PER_DEGREE_LONGITUDE_AT_EQUATOR = 111.320
EARTH_RADIUS_KM = 6371.0088


def cell_area_km2(lat):
    """Area of a one-degree grid cell centred at latitude `lat` (degrees)."""
    return (KM_PER_DEGREE_LONGITUDE_AT_EQUATOR * np.cos(np.radians(lat)) *
            KM_PER_DEGREE_LATITUDE)


def haversine_km(lon1, lat1, lon2, lat2):
    """Great-circle distance in km."""
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    dlon, dlat = lon2 - lon1, lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return EARTH_RADIUS_KM * 2 * np.arcsin(np.sqrt(a))


def zones_and_area_weighted_rank(intensity, area, zone_fractions):
    """Area-weighted rank of every cell and the cells of each zone.

    Returns (area_weighted_rank, zone_cells, zone_area_fraction), where
    zone_cells[fraction] is a True/False array over the cells and
    zone_area_fraction[fraction] is the area share the zone actually covers.
    """
    total_area = float(area.sum())
    ascending = np.argsort(intensity)
    area_weighted_rank = np.empty_like(intensity)
    area_weighted_rank[ascending] = np.cumsum(area[ascending]) / total_area
    descending = ascending[::-1]
    cumulative_share = np.cumsum(area[descending]) / total_area
    zone_cells, zone_area_fraction = {}, {}
    for fraction in zone_fractions:
        n_cells = int(np.searchsorted(cumulative_share, fraction)) + 1
        in_zone = np.zeros(intensity.size, dtype=bool)
        in_zone[descending[:n_cells]] = True
        zone_cells[fraction] = in_zone
        zone_area_fraction[fraction] = float(area[in_zone].sum() / total_area)
    return area_weighted_rank, zone_cells, zone_area_fraction


def distance_to_zone_km(centroid_lon, centroid_lat, centroid_cell_row, in_zone,
                        cell_lon, cell_lat):
    """0 km inside the zone, else the distance to the nearest zone cell centre."""
    zone_lon, zone_lat = cell_lon[in_zone], cell_lat[in_zone]
    distance = np.zeros(centroid_lon.size)
    for j in range(centroid_lon.size):
        if in_zone[centroid_cell_row[j]]:
            continue
        distance[j] = float(np.min(haversine_km(centroid_lon[j], centroid_lat[j],
                                                zone_lon, zone_lat)))
    return distance


def zone_label(fraction):
    return f"top{int(round(fraction * 100))}"


def score_centroids(intensity, cells, centroids, zone_fractions):
    """Scores of the given centroids on one mean-intensity surface.

    intensity : mean intensity per km^2 of each cell, in the row order of `cells`
    cells     : table with lon, lat, area_km2 (one row per grid cell)
    centroids : table with ellipse_id, centroid_lon, centroid_lat, cell_row
                (row of the centroid's nearest cell in `cells`)
    Returns (per_centroid table, zone table).
    """
    area = cells["area_km2"].to_numpy(float)
    cell_lon = cells["lon"].to_numpy(float)
    cell_lat = cells["lat"].to_numpy(float)
    area_weighted_rank, zone_cells, zone_area_fraction = zones_and_area_weighted_rank(
        intensity, area, zone_fractions)
    cell_rank = pd.Series(intensity).rank(method="average", pct=True).to_numpy()
    expected_share = intensity * area
    log_density = np.log((expected_share / expected_share.sum()) / area)

    rows = centroids["cell_row"].to_numpy(int)
    lon = centroids["centroid_lon"].to_numpy(float)
    lat = centroids["centroid_lat"].to_numpy(float)
    per_centroid = pd.DataFrame({
        "ellipse_id": centroids["ellipse_id"].to_numpy(),
        "area_weighted_rank": area_weighted_rank[rows],
        "cell_rank": cell_rank[rows],
    })
    for fraction in zone_fractions:
        per_centroid[f"in_{zone_label(fraction)}_zone"] = zone_cells[fraction][rows]
    for fraction in zone_fractions:
        per_centroid[f"distance_to_{zone_label(fraction)}_zone_km"] = distance_to_zone_km(
            lon, lat, rows, zone_cells[fraction], cell_lon, cell_lat)
    per_centroid["log_density"] = log_density[rows]

    zones = pd.DataFrame([{
        "zone_area_fraction": fraction,
        "n_cells_in_zone": int(zone_cells[fraction].sum()),
        "zone_area_fraction_covered": zone_area_fraction[fraction],
        "area_overshoot_fraction": zone_area_fraction[fraction] - fraction,
        "lowest_area_weighted_rank_in_zone":
            float(area_weighted_rank[zone_cells[fraction]].min()),
    } for fraction in zone_fractions])
    return per_centroid, zones


def summarize_groups(per_centroid, zone_fractions):
    """Group summaries for all centroids, daily-largest and event-largest.

    per_centroid needs the score columns plus is_daily_largest and
    is_event_largest.
    """
    groups = [("all_centroids", np.ones(len(per_centroid), dtype=bool)),
              ("daily_largest", per_centroid["is_daily_largest"].to_numpy(bool)),
              ("event_largest", per_centroid["is_event_largest"].to_numpy(bool))]
    rows = []
    for name, selected in groups:
        group = per_centroid[selected]
        row = {"group": name, "n_centroids": len(group),
               "mean_area_weighted_rank": group["area_weighted_rank"].mean(),
               "median_area_weighted_rank": group["area_weighted_rank"].median(),
               "min_area_weighted_rank": group["area_weighted_rank"].min(),
               "max_area_weighted_rank": group["area_weighted_rank"].max(),
               "mean_cell_rank": group["cell_rank"].mean()}
        for fraction in zone_fractions:
            label = zone_label(fraction)
            row[f"n_in_{label}_zone"] = int(group[f"in_{label}_zone"].sum())
            row[f"share_in_{label}_zone"] = group[f"in_{label}_zone"].mean()
        for fraction in zone_fractions:
            label = zone_label(fraction)
            row[f"mean_distance_to_{label}_zone_km"] = group[f"distance_to_{label}_zone_km"].mean()
        row["median_distance_to_top20_zone_km"] = group["distance_to_top20_zone_km"].median()
        row["mean_log_density"] = group["log_density"].mean()
        rows.append(row)
    return pd.DataFrame(rows)
