"""Prepare the centroid locations and temperature covariates for the LGCP.

Reads: Processed temperatures, catalogs, and prepared maximum-area tables.
Writes: supporting_files/additional_results/modeling/lgcp/
Run: python scripts/modeling/lgcp/01_prepare_input.py
"""
from __future__ import annotations

import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
from sklearn.neighbors import NearestNeighbors

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402

# Fixed settings for preparing the covariates and spatial window.
CSV_BLOCK_ROWS = 2_000_000   # rows per block; fixes the floating-point summation order
WINDOW_MARGIN_KM = 80.0      # a pixel is inside if within 25 + 80 km of a cell centre

INPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp"
POWER_LAW_INPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "power_law"


def tmax_as_read_by_pandas(values):
    """Parse daily Tmax from decimal text for the manuscript covariates.

    Write each value as its shortest exact decimal, then read it with
    pandas' default number reader. This can change the last binary digit.
    Keep this parsing step and float64 precision for the covariates; the
    float32 rounding used for centroid weights does not apply here.
    """
    text = "v\n" + "\n".join(repr(float(v)) for v in values)
    return pd.read_csv(io.StringIO(text))["v"].to_numpy()


def temperature_covariates(config):
    """Mean and population standard deviation of daily Tmax per grid cell."""
    tmax_input = f.read_tmax_input(config)
    f.check_study_calendar(tmax_input["dates"], config)
    tmax = tmax_input["tmax"]                       # (date, lat, lon)
    lat, lon = tmax_input["lat"], tmax_input["lon"]
    n_dates, n_lat, n_lon = tmax.shape

    # Sum values in latitude, longitude, then date order.
    tmax_in_csv_order = np.transpose(tmax, (1, 2, 0)).reshape(-1)
    records = pd.DataFrame({
        "lat1": np.repeat(lat, n_lon * n_dates),
        "lon1": np.tile(np.repeat(lon, n_dates), n_lat),
        "val": tmax_as_read_by_pandas(tmax_in_csv_order),
    })

    # Accumulate the count, sum and squared sum within each block.
    n_days, sum_tmax, sum_tmax_squared = {}, {}, {}
    for start in range(0, len(records), CSV_BLOCK_ROWS):
        block = records.iloc[start:start + CSV_BLOCK_ROWS].dropna(subset=["val"])
        by_cell = block.groupby(["lon1", "lat1"])["val"]
        block_count = by_cell.count()
        block_sum = by_cell.sum()
        block_sum_squared = by_cell.apply(lambda v: float(np.sum(np.square(v))))
        for key, value in block_count.items():
            n_days[key] = n_days.get(key, 0) + int(value)
        for key, value in block_sum.items():
            sum_tmax[key] = sum_tmax.get(key, 0.0) + float(value)
        for key, value in block_sum_squared.items():
            sum_tmax_squared[key] = sum_tmax_squared.get(key, 0.0) + float(value)

    rows = []
    for (cell_lon, cell_lat), count in n_days.items():
        mean = sum_tmax[(cell_lon, cell_lat)] / count
        variance = max(0.0, sum_tmax_squared[(cell_lon, cell_lat)] / count - mean * mean)
        rows.append((float(cell_lon), float(cell_lat), int(count), mean,
                     float(np.sqrt(variance))))
    cells = pd.DataFrame(rows, columns=["lon", "lat", "n_days", "mean_tmax", "std_tmax"])
    cells = cells.sort_values(["lat", "lon"]).reset_index(drop=True)
    return cells, tmax_input


def main():
    started = time.time()
    config = f.load_config()
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    grid_cells = f.read_grid_cells(config)
    ellipse_file = f.catalog_dir(config, "04_extreme_heatwave_ellipses") / "ellipses.csv"
    event_file = f.catalog_dir(config, "05_compound_heatwave_events") / "compound_events.csv"
    event_days_file = event_file.with_name("event_days.csv")
    daily_max_file = POWER_LAW_INPUT_DIR / "daily_max_area.csv"
    event_max_file = POWER_LAW_INPUT_DIR / "event_max_area.csv"
    for needed in (daily_max_file, event_max_file):
        if not needed.exists():
            raise FileNotFoundError(f"{needed} is missing: run "
                                    "scripts/modeling/power_law/01_prepare_input.py first")

    # 1. temperature covariates
    cells, tmax_input = temperature_covariates(config)
    if len(cells) != config["study_area"]["expected_cells"]:
        raise ValueError(f"{len(cells)} grid cells, expected 1800")
    # The cells come out in the order of grid_cells.csv (latitude, then longitude).
    if not (np.array_equal(cells["lat"], grid_cells["lat"]) and
            np.array_equal(cells["lon"], grid_cells["lon"])):
        raise ValueError("covariate cells are not in grid_cells.csv order")
    cells.insert(0, "cell_id", grid_cells["cell_id"].to_numpy())
    print(f"temperature covariates: {len(cells)} cells, "
          f"{cells['n_days'].min()}-{cells['n_days'].max()} days each "
          f"({time.time() - started:.0f} s)")

    # 2. projection
    # Centre the projection on the grid-cell domain (45E, 28N).
    lon0 = float((cells["lon"].min() + cells["lon"].max()) / 2.0)
    lat0 = float((cells["lat"].min() + cells["lat"].max()) / 2.0)
    projection = (f"+proj=laea +lat_0={lat0} +lon_0={lon0} +x_0=0 +y_0=0 "
                  f"+ellps=WGS84 +units=m +no_defs")
    to_map = Transformer.from_crs("EPSG:4326", CRS.from_proj4(projection), always_xy=True)
    to_lonlat = Transformer.from_crs(to_map.target_crs, "EPSG:4326", always_xy=True)

    def project_km(lon, lat):
        x, y = to_map.transform(np.asarray(lon), np.asarray(lat))
        return np.asarray(x) / 1000.0, np.asarray(y) / 1000.0

    cells["x_km"], cells["y_km"] = project_km(cells["lon"], cells["lat"])
    cells = cells[["cell_id", "lon", "lat", "x_km", "y_km", "n_days",
                   "mean_tmax", "std_tmax"]]

    # 3. centroids
    # Use pandas' default number reader for centroid coordinates. It can
    # change the last binary digit; the manuscript's projected coordinates
    # use those parsed values. Read areas exactly with round_trip. (The
    # power-law samples instead keep the default reader's areas; see
    # power_law/01_prepare_input.py.)
    ellipses = f.attach_compound_events(pd.read_csv(ellipse_file),
                                       pd.read_csv(event_days_file))
    ellipses["area_km2"] = pd.read_csv(ellipse_file,
                                       float_precision="round_trip")["area_km2"]
    events = pd.read_csv(event_file)
    if len(ellipses) != 753:
        raise ValueError(f"{len(ellipses)} ellipses, expected 753")
    centroid_x, centroid_y = project_km(ellipses["centroid_lon_weighted"],
                                        ellipses["centroid_lat_weighted"])
    nearest = NearestNeighbors(n_neighbors=1).fit(cells[["x_km", "y_km"]].to_numpy())
    distance, index = nearest.kneighbors(np.column_stack([centroid_x, centroid_y]))
    nearest_cell = cells.iloc[index.ravel()].reset_index(drop=True)

    # These shared ellipse lists label concentration groups; no fitted
    # power-law parameters are used by the LGCP model.
    daily_max = pd.read_csv(daily_max_file)
    event_max = pd.read_csv(event_max_file)
    centroids = pd.DataFrame({
        "ellipse_id": ellipses["ellipse_id"],
        "cluster_id": ellipses["cluster_id"],
        "compound_event_id": ellipses["compound_event_id"],
        "date": ellipses["date"],
        "event_type": ellipses["compound_event_id"].map(
            events.set_index("compound_event_id")["event_type"]),
        "centroid_lon": ellipses["centroid_lon_weighted"],
        "centroid_lat": ellipses["centroid_lat_weighted"],
        "x_km": centroid_x,
        "y_km": centroid_y,
        "nearest_cell_id": nearest_cell["cell_id"],
        "nearest_cell_lon": nearest_cell["lon"],
        "nearest_cell_lat": nearest_cell["lat"],
        "nearest_cell_distance_km": distance.ravel(),
        "mean_tmax": nearest_cell["mean_tmax"],
        "std_tmax": nearest_cell["std_tmax"],
        "area_km2": ellipses["area_km2"],
        "is_daily_largest": ellipses["ellipse_id"].isin(daily_max["ellipse_id"]),
        "is_event_largest": ellipses["ellipse_id"].isin(event_max["ellipse_id"]),
    })
    if centroids["is_daily_largest"].sum() != 395 or centroids["is_event_largest"].sum() != 51:
        raise ValueError("daily/event-largest flags do not give 395 and 51 centroids")
    if (centroids["is_event_largest"] & ~centroids["is_daily_largest"]).any():
        raise ValueError("an event-largest centroid is not a daily-largest centroid")

    # 4. 25-km raster and study window
    pixel_km = float(config["lgcp"]["raster_resolution_km"])
    x_first, x_last = cells["x_km"].min() - pixel_km, cells["x_km"].max() + pixel_km
    y_first, y_last = cells["y_km"].min() - pixel_km, cells["y_km"].max() + pixel_km
    pixel_x = np.arange(x_first + pixel_km / 2, x_last, pixel_km)
    pixel_y = np.arange(y_first + pixel_km / 2, y_last, pixel_km)
    grid_x, grid_y = np.meshgrid(pixel_x, pixel_y)          # shape (ny, nx)
    pixels = np.column_stack([grid_x.ravel(), grid_y.ravel()])
    pixel_distance, pixel_index = nearest.kneighbors(pixels)
    inside = pixel_distance.ravel() <= (pixel_km + WINDOW_MARGIN_KM)
    pixel_cell = cells.iloc[pixel_index.ravel()].reset_index(drop=True)
    pixel_lon, pixel_lat = to_lonlat.transform(pixels[:, 0] * 1000.0, pixels[:, 1] * 1000.0)
    raster = pd.DataFrame({
        "ix": np.tile(np.arange(len(pixel_x)), len(pixel_y)),
        "iy": np.repeat(np.arange(len(pixel_y)), len(pixel_x)),
        "x_km": pixels[:, 0],
        "y_km": pixels[:, 1],
        "inside": inside.astype(int),
        "nearest_cell_id": pixel_cell["cell_id"],
        "lon": pixel_cell["lon"],
        "lat": pixel_cell["lat"],
        "mean_tmax": pixel_cell["mean_tmax"],
        "std_tmax": pixel_cell["std_tmax"],
        "pixel_lon": pixel_lon,
        "pixel_lat": pixel_lat,
    })
    if not (distance.ravel() <= pixel_km + WINDOW_MARGIN_KM).all():
        raise ValueError("a centroid lies outside the study window")

    window = {
        "projection": projection,
        "projection_centre_lon": lon0,
        "projection_centre_lat": lat0,
        "units": "km",
        "pixel_km": pixel_km,
        "nx": int(len(pixel_x)),
        "ny": int(len(pixel_y)),
        "first_pixel_centre_x_km": float(pixel_x[0]),
        "first_pixel_centre_y_km": float(pixel_y[0]),
        "n_pixels_inside": int(inside.sum()),
        "inside_rule_km": pixel_km + WINDOW_MARGIN_KM,
        "inside_rule": "pixel centre within 105 km of the nearest grid-cell centre",
        "n_centroids": int(len(centroids)),
        "n_cells": int(len(cells)),
        "tmax_period": [str(tmax_input["dates"][0]), str(tmax_input["dates"][-1])],
    }

    # write
    covariate_file = INPUT_DIR / "temperature_covariates.csv"
    centroid_file = INPUT_DIR / "centroids.csv"
    raster_file = INPUT_DIR / "covariate_raster_25km.csv"
    window_file = INPUT_DIR / "study_window.json"
    cells.to_csv(covariate_file, index=False, lineterminator="\n")
    centroids.to_csv(centroid_file, index=False, lineterminator="\n")
    raster.to_csv(raster_file, index=False, lineterminator="\n")
    f.write_json(window_file, window)

    print(f"finished in {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
