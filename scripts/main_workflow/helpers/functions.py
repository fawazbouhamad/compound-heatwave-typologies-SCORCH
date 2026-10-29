"""Calculations shared by the analysis, models and supplementary scripts.

Study settings and file helpers come first. Calculations are grouped by stage,
with DBSCAN summaries alongside clustering and geographic helpers at the end.
This file is imported, not run.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parents[3]
NOISE = -1   # DBSCAN label for cells that belong to no cluster


# Study settings
def load_config():
    """Return a fresh copy of the fixed manuscript settings.

    Values marked informational describe the study and are not read by scripts.
    """
    return {
        "paths": {
            "input_dir": "input",
            "supporting_input_dir": "supporting_files/additional_inputs",
            "supporting_results_dir": "supporting_files/additional_results",
            "output_dir": "output",
        },
        # Bounds are informational. The input supplies the grid; counts are checked.
        "study_area": {
            "lat_min": 10.0,
            "lat_max": 46.0,
            "lon_min": 20.0,
            "lon_max": 70.0,
            "grid_step_deg": 1.0,
            "warm_season_months": [4, 5, 6, 7, 8, 9],
            "first_year": 1940,
            "last_year": 2025,
            "expected_dates": 15738,
            "expected_cells": 1800,
        },
        # Linear percentile; Tmax >= threshold. Isolated gaps are bridged.
        "local_heatwaves": {
            "threshold_quantile": 0.95,
            "min_duration_days": 3,
            "min_exceeding_days": 3,
        },
        # Includes zero-count days; the higher method selects an observed value.
        "regional_days": {
            "quantile": 0.975,
            "quantile_method": "higher",
        },
        # Neighbourhood distances are in grid-cell units.
        "clustering": {
            "eps_values": [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            "min_samples_values": [4, 5, 6, 7, 8, 9, 10, 11, 12],
        },
        # Only centroid weights use float32; other steps retain input precision.
        "ellipses": {
            "scale_factor": 1.25,
            "eigenvalue_floor_km2": 1.0,
            "weighted_centroid_tmax_precision": "float32",
        },
        "representative_events": {
            "figure_5_type_1": ["2002-07-31"],
            "figure_5_type_2": ["1983-07-13"],
            "figure_6a_type_3": ["2017-07-17", "2017-07-18"],
            "figure_6b_type_3": ["1991-06-04", "1991-06-05"],
            "figure_7_type_4": [
                "2001-08-07", "2001-08-08", "2001-08-09",
                "2001-08-10", "2001-08-11", "2001-08-12",
            ],
            "figure_s1_type_3": ["2016-06-07", "2016-06-08"],
        },
        # Clauset fit: all unique cutoffs except the largest; smaller cutoff on ties.
        # Uses the MLE exponent and full two-sided KS distance.
        # The significance level is informational.
        "power_law": {
            "ks_mode": "full",
            "n_bootstrap": 5000,
            "seed_daily_maxima": 20261020,
            "seed_event_maxima": 20261120,
            "significance_level": 0.1,
        },
        # Projection is informational; the LGCP preparation derives it from the data.
        # Lambert azimuthal equal-area, WGS84, centred at 45 E / 28 N.
        "lgcp": {
            "projection": (
                "+proj=laea +lat_0=28 +lon_0=45 +x_0=0 +y_0=0 "
                "+ellps=WGS84 +units=m +no_defs"
            ),
            "raster_resolution_km": 25,
            "trend_terms": ["lon", "lat", "mean_tmax", "std_tmax"],
            "seed": 20260704,
        },
        # Scale factors are tested on the two axes independently.
        "appendix_b": {
            "scale_factors": [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5],
        },
        "appendix_c": {
            "n_folds": 5,
            "fold_seed": 20260704,
            "zone_area_fractions": [0.1, 0.2, 0.3, 0.5],
        },
        "supplementary_power_law": {
            "n_bootstrap": 5000,
            "base_seed": 20260918,
        },
    }


# File paths
def project_path(relative_or_absolute):
    """Turn a configured repository path into an absolute path."""
    path = Path(relative_or_absolute)
    return path if path.is_absolute() else (REPO_DIR / path).resolve()


def input_dir(config):
    return project_path(config["paths"]["input_dir"])


def supporting_input_dir(config):
    return project_path(config["paths"]["supporting_input_dir"])


def supporting_results_dir(config):
    return project_path(config["paths"]["supporting_results_dir"])


def catalog_dir(config, stage_folder):
    """Folder of one catalog stage, created if needed."""
    folder = project_path(config["paths"]["output_dir"]) / "catalogs" / stage_folder
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def output_subdir(config, name):
    folder = project_path(config["paths"]["output_dir"]) / name
    folder.mkdir(parents=True, exist_ok=True)
    return folder


# Run records
def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def software_versions():
    versions = {"python": platform.python_version()}
    for name in ("numpy", "pandas", "scipy", "sklearn", "xarray", "netCDF4",
                 "matplotlib", "cartopy", "yaml"):
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", "unknown")
        except ImportError:
            versions[name] = "not installed"
    return versions


def write_json(path, payload):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(payload, fh, indent=2, default=str)
        fh.write("\n")


def now_utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def relative_to_repo(path):
    path = Path(path).resolve()
    try:
        return path.relative_to(REPO_DIR).as_posix()
    except ValueError:
        return str(path)


def write_run_record(config, script_name, inputs, outputs, notes=None):
    """Record what a script read and wrote, with checksums and versions."""
    record = {
        "script": script_name,
        "finished_utc": now_utc(),
        "software": software_versions(),
        "inputs": {relative_to_repo(p): sha256_of(p) for p in inputs},
        "outputs": {relative_to_repo(p): sha256_of(p) for p in outputs},
    }
    if notes:
        record["notes"] = notes
    folder = REPO_DIR / ".scorch/logs/run-records"
    write_json(folder / f"{Path(script_name).stem}.json", record)


# Processed temperature input
def read_tmax_input(config):
    """Read input/tmax_processed.nc.

    Returns a dict with
      dates : array of "YYYY-MM-DD" strings (the 15,738 warm-season dates)
      lat, lon : cell-centre coordinates in degrees (ascending)
      tmax : float64 array (date, lat, lon), degrees Celsius
      path : the file that was read
    """
    from netCDF4 import Dataset, num2date

    path = input_dir(config) / "tmax_processed.nc"
    with Dataset(path) as ds:
        time = ds.variables["time"]
        dates = np.array([d.strftime("%Y-%m-%d")
                          for d in num2date(time[:], time.units, time.calendar)])
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        variable = ds.variables["tmax"]
        variable.set_auto_mask(False)
        tmax = np.asarray(variable[:], dtype=np.float64)
        units = variable.units
    if units != "degC":
        raise ValueError(f"tmax units are {units!r}, expected degC")
    if not np.all(np.isfinite(tmax)):
        raise ValueError("tmax contains missing or non-finite values")
    return {"dates": dates, "lat": lat, "lon": lon, "tmax": tmax, "path": path}


def centroid_temperature_weights(values):
    """Prepare the single-precision temperature weights used for centroids.

    Convert to decimal text, read with pandas' default parser, then round to
    float32. Keeping this exact sequence preserves the manuscript's weights.
    """
    import io
    text = "v\n" + "\n".join(repr(float(v)) for v in np.ravel(values))
    parsed = pd.read_csv(io.StringIO(text))["v"].to_numpy()
    return parsed.astype(np.float32).reshape(np.shape(values))


def read_grid_cells(config):
    return pd.read_csv(supporting_input_dir(config) / "maps" / "grid_cells.csv")


def expected_study_dates(config):
    """The April-September days of 1940-2025 as "YYYY-MM-DD" strings."""
    area = config["study_area"]
    days = pd.date_range(f"{area['first_year']}-01-01",
                         f"{area['last_year']}-12-31", freq="D")
    days = days[days.month.isin(area["warm_season_months"])]
    return days.strftime("%Y-%m-%d").to_numpy()


def check_study_calendar(dates, config):
    """Stop unless the dates are exactly the expected warm-season calendar."""
    expected = expected_study_dates(config)
    if len(dates) != config["study_area"]["expected_dates"] or \
            not np.array_equal(np.asarray(dates), expected):
        raise ValueError("the input dates are not the expected warm-season calendar")


# Stage 1: cell thresholds, exceedance and local heatwaves
def cell_thresholds(tmax, quantile):
    """Empirical percentile of each cell's Tmax over all study dates.

    Use linear interpolation between ordered values, with float64 Tmax.
    """
    n_time, n_lat, n_lon = tmax.shape
    values = pd.DataFrame({
        "cell": np.repeat(np.arange(n_lat * n_lon), n_time),
        "val": tmax.reshape(n_time, -1).T.ravel(),
    })
    thresholds = values.groupby("cell", sort=True)["val"].quantile(quantile)
    return thresholds.to_numpy().reshape(n_lat, n_lon)


def label_heatwaves(exceed_series, min_len=3, min_ones=3,
                    allow_unlimited_single_zeros=True):
    """Label local heatwaves in one 0/1 daily exceedance series.

    * split the series wherever two consecutive days do not exceed,
    * trim each piece so it starts and ends on an exceeding day,
    * keep it when it spans at least ``min_len`` days and contains at least
      ``min_ones`` exceeding days; single non-exceeding days inside it are
      part of the heatwave.

    Returns (labels, episodes): the heatwave number of every day (0 = not in
    a heatwave) and a list of (number, start, end, length, n_exceeding,
    n_bridged) tuples.
    """
    s = pd.Series(exceed_series).astype(int)
    idx = pd.DatetimeIndex(s.index)
    vals = s.values
    zero = vals == 0

    # Position of the first zero of every "00" pair.
    double_zero_breaks = np.where(zero[:-1] & zero[1:])[0]
    boundaries = [-1] + double_zero_breaks.tolist() + [len(s) - 1]
    raw_segments = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        i0, i1 = a + 1, b
        if i0 <= i1:
            raw_segments.append((i0, i1))

    hw_id = 0
    labels = np.zeros(len(s), dtype=int)
    episodes = []
    for i0, i1 in raw_segments:
        seg_vals = vals[i0:i1 + 1]
        if seg_vals.sum() == 0:
            continue
        ones_pos = np.where(seg_vals == 1)[0]
        j0 = i0 + int(ones_pos[0])
        j1 = i0 + int(ones_pos[-1])
        event_vals = vals[j0:j1 + 1]
        ones_count = int(event_vals.sum())
        zeros_count = int((event_vals == 0).sum())
        length_days = int(j1 - j0 + 1)
        if allow_unlimited_single_zeros:
            ok = length_days >= min_len and ones_count >= min_ones
        else:
            ok = (length_days >= min_len and ones_count >= min_ones
                  and zeros_count <= 1)
        if ok:
            hw_id += 1
            labels[j0:j1 + 1] = hw_id
            episodes.append((hw_id, idx[j0].date(), idx[j1].date(),
                             length_days, ones_count, zeros_count))
    return labels, episodes


def label_cell_by_season(exceed, dates, min_len, min_ones):
    """Apply label_heatwaves to one cell, one warm season at a time.

    30 September and the following 1 April are not consecutive days, so a
    heatwave may not continue across that gap. Each season is labelled
    separately, and heatwave numbers keep counting upward through the
    record so every heatwave in a cell has its own number (1, 2, 3, ...).
    """
    years = np.array([int(d[:4]) for d in dates])
    labels = np.zeros(len(exceed), dtype=np.int32)
    episodes = []
    offset = 0
    for year in np.unique(years):
        rows = np.where(years == year)[0]
        season = pd.Series(exceed[rows], index=pd.to_datetime(dates[rows]))
        season_labels, season_episodes = label_heatwaves(
            season, min_len=min_len, min_ones=min_ones)
        labels[rows] = np.where(season_labels > 0, season_labels + offset, 0)
        for number, start, end, length, n_exceed, n_bridged in season_episodes:
            episodes.append((number + offset, start, end, length, n_exceed, n_bridged))
        offset += len(season_episodes)
    return labels, episodes


# Stage 2: regionally extensive days and consecutive-day events
def regional_threshold(daily_heatwave_cells, quantile, method):
    """Regional threshold from heatwave-cell counts on all study dates.

    The "higher" quantile method returns the observed count at or just above
    the quantile position.
    """
    counts = pd.Series(daily_heatwave_cells)
    return float(counts.quantile(quantile, interpolation=method))


def group_consecutive_days(selected_dates):
    """Identify consecutive-date sequences by their start date.

    Step 3 pools DBSCAN settings within each sequence. Compound event IDs
    are assigned only in step 5, after clustering and ellipse fitting.
    """
    dates = pd.DatetimeIndex(sorted(pd.to_datetime(list(selected_dates))))
    gaps = np.diff(dates.values).astype("timedelta64[D]").astype(int)
    new_sequence = np.concatenate([[True], gaps > 1])
    starts = np.maximum.accumulate(np.where(new_sequence, np.arange(len(dates)), 0))
    return pd.DataFrame({"date": dates.strftime("%Y-%m-%d"),
                         "sequence_start_date": dates[starts].strftime("%Y-%m-%d")})


def attach_compound_events(records, event_days):
    """Link dated records to the completed step-5 event-day catalog in memory."""
    if "compound_event_id" in records:
        raise ValueError("records already contain compound event IDs")
    lookup = event_days.set_index("date")["compound_event_id"]
    if not lookup.index.is_unique:
        raise ValueError("the event-day catalog contains duplicate dates")
    result = records.copy()
    result["compound_event_id"] = result["date"].map(lookup)
    if result["compound_event_id"].isna().any():
        raise ValueError("a record date is missing from the step-5 event-day catalog")
    result["compound_event_id"] = result["compound_event_id"].astype(int)
    return result


# Stage 3: DBSCAN clustering of each selected day
def longitude_major_order(lon, lat):
    """Sort cells by longitude, then latitude (west to east, south to north).

    The order matters for DBSCAN: a border cell within reach of two clusters
    joins the cluster that reaches it first. Labels follow the order in which
    clusters are first found.
    """
    return np.lexsort((np.asarray(lat), np.asarray(lon)))


def grid_distance_matrix(lon, lat, grid_step=1.0):
    """Euclidean distances between cells in grid-step units.

    One unit is one grid step (1 degree).
    """
    from scipy.spatial.distance import cdist
    pts = np.column_stack((np.asarray(lon, dtype=float) / float(grid_step),
                           np.asarray(lat, dtype=float) / float(grid_step)))
    return cdist(pts, pts, metric="euclidean").astype(float)


def run_dbscan(distances, eps, min_samples):
    """DBSCAN on a precomputed distance matrix (scikit-learn).

    A cell is a core cell when at least min_samples cells, itself included,
    lie within eps. Returns (labels, is_core): labels are 0, 1, 2, ... in
    order of first appearance and -1 for noise.
    """
    from sklearn.cluster import DBSCAN
    model = DBSCAN(eps=float(eps), min_samples=int(min_samples),
                   metric="precomputed").fit(distances)
    is_core = np.zeros(len(model.labels_), dtype=bool)
    is_core[model.core_sample_indices_] = True
    return model.labels_, is_core


def count_clusters(labels):
    return len(set(np.asarray(labels).tolist()) - {NOISE})


def parameter_search(distances, eps_values, min_samples_values):
    """Number of clusters for every (eps, min_samples) pair on one day."""
    rows = []
    for eps in eps_values:
        for min_samples in min_samples_values:
            labels, _ = run_dbscan(distances, eps, min_samples)
            rows.append({"eps": float(eps), "min_samples": int(min_samples),
                         "n_clusters": count_clusters(labels),
                         "n_noise_cells": int((labels == NOISE).sum())})
    return rows


def round_half_up(x):
    """Nearest integer with halves rounded up: floor(x + 0.5).

    This is the rounding used for the daily minimum sample size. It is not
    the same as Python's round(), which rounds halves to the nearest even
    number.
    """
    return int(np.floor(float(x) + 0.5))


def select_daily_parameters(search_rows):
    """Choose one day's DBSCAN settings from its parameter search.

    Find the most frequent nonzero cluster count (all counts if every pair
    gives zero clusters). If several counts tie, keep the pairs of all of them;
    daily eps is the mean retained eps and daily min_samples the mean
    retained min_samples rounded half up (never below 1).
    """
    counts = [r["n_clusters"] for r in search_rows]
    nonzero = [c for c in counts if c > 0]
    pool = nonzero if nonzero else counts
    frequency = pd.Series(pool).value_counts()
    top = int(frequency.max())
    modal_counts = sorted(int(v) for v, n in frequency.items() if n == top)
    kept = [r for r in search_rows if r["n_clusters"] in modal_counts]
    mean_min_samples = float(np.mean([r["min_samples"] for r in kept]))
    return {
        "modal_cluster_counts": modal_counts,
        "modal_frequency": top,
        "n_pairs_retained": len(kept),
        "daily_eps": float(np.mean([r["eps"] for r in kept])),
        "daily_min_samples_mean": mean_min_samples,
        "daily_min_samples": max(1, round_half_up(mean_min_samples)),
        "retained_pairs": {(r["eps"], r["min_samples"]) for r in kept},
    }


# DBSCAN summaries for Appendix A
DBSCAN_GROUPS = ["All Types", "Type 1", "Type 2", "Type 3", "Type 4"]


def values_by_group(frame, column):
    groups = {"All Types": frame[column].to_numpy(float)}
    for group in DBSCAN_GROUPS[1:]:
        groups[group] = frame.loc[frame["type"] == group, column].to_numpy(float)
    return groups


def box_statistics(values):
    """Box-plot statistics; quartiles by numpy's default linear interpolation."""
    return {
        "n": len(values),
        "mean": values.mean(),
        "median": float(np.median(values)),
        "min": values.min(),
        "max": values.max(),
        "q1": np.percentile(values, 25),
        "q3": np.percentile(values, 75),
    }


def figure_a2_panels(per_day, per_event):
    """(values by group, y-axis label, population) for panels (a)-(d)."""
    min_samples_label = "Selected DBSCAN min_samples"
    eps_label = "Selected DBSCAN eps (grid-cell units)"
    days_note = f"{len(per_day)} event-days"
    events_note = f"{len(per_event)} events"
    return {
        "a": (values_by_group(per_day, "min_samples"), min_samples_label, days_note),
        "b": (values_by_group(per_event, "min_samples"), min_samples_label, events_note),
        "c": (values_by_group(per_day, "eps"), eps_label, days_note),
        "d": (values_by_group(per_event, "eps"), eps_label, events_note),
    }


def statistics_table(panels):
    rows = []
    for letter, (groups, ylabel, population) in panels.items():
        for group in DBSCAN_GROUPS:
            rows.append({
                "panel": letter,
                "metric": ylabel,
                "population": population,
                "group": group,
                **box_statistics(groups[group]),
            })
    return pd.DataFrame(rows)


# Stage 4: PCA ellipse geometry and temperature-weighted centroid
KM_PER_DEG_LON_AT_EQUATOR = 111.320
KM_PER_DEG_LAT = 110.574


def project_lonlat_to_km(lon, lat, lon0, lat0):
    """Local equirectangular projection to km around (lon0, lat0).

    One degree of latitude is 110.574 km and one degree of longitude is
    111.320 km times cos(lat0).
    """
    km_per_deg_lon = KM_PER_DEG_LON_AT_EQUATOR * math.cos(math.radians(lat0))
    x = (np.asarray(lon, dtype=float) - lon0) * km_per_deg_lon
    y = (np.asarray(lat, dtype=float) - lat0) * KM_PER_DEG_LAT
    return np.column_stack((x, y))


def fit_pca(points_xy):
    """Sample-covariance PCA of km coordinates.

    Returns (centre, eigenvectors, eigenvalues) with eigenvalues in
    descending order (clipped at zero); eigenvector column 0 is the major
    axis. A single-cell cluster gets unit eigenvalues.
    """
    points_xy = np.asarray(points_xy, dtype=float)
    if points_xy.shape[0] == 1:
        return points_xy[0].copy(), np.eye(2), np.array([1.0, 1.0])
    center = points_xy.mean(axis=0)
    cov = np.cov(points_xy - center, rowvar=False, bias=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    return center, eigvecs[:, order], np.clip(eigvals[order], 0.0, None)


def ellipse_semi_axes_km(eigvals, scale_factor, eigenvalue_floor_km2):
    """Semi-axes a = s*sqrt(max(lambda1, floor)), b = s*sqrt(max(lambda2, floor))."""
    a = scale_factor * math.sqrt(max(float(eigvals[0]), eigenvalue_floor_km2))
    b = scale_factor * math.sqrt(max(float(eigvals[1]), eigenvalue_floor_km2))
    return a, b


def major_axis_orientation_east_ccw(eigvecs):
    """Major-axis angle counterclockwise from east, in [0, 180) degrees."""
    return math.degrees(math.atan2(float(eigvecs[1, 0]), float(eigvecs[0, 0]))) % 180.0


def azimuth_from_north(orientation_east_ccw):
    """Convert to the manuscript convention: clockwise from north, (-90, 90].

    theta_north = 90 - theta_east, wrapped so that an exact east-west axis
    is +90.
    """
    out = np.mod(90.0 - np.asarray(orientation_east_ccw, dtype=float) + 90.0, 180.0) - 90.0
    return np.where(out <= -90.0, out + 180.0, out)


def ellipse_geometry(lon, lat, scale_factor, eigenvalue_floor_km2):
    """PCA ellipse of one cluster (unweighted, centred on the cell mean).

    L1 = 2 s sqrt(max(lambda1, floor)), L2 = 2 s sqrt(max(lambda2, floor)),
    area = pi L1 L2 / 4, ratio = L2 / L1.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon0, lat0 = float(lon.mean()), float(lat.mean())
    center_xy, eigvecs, eigvals = fit_pca(project_lonlat_to_km(lon, lat, lon0, lat0))
    a, b = ellipse_semi_axes_km(eigvals, scale_factor, eigenvalue_floor_km2)
    orientation = major_axis_orientation_east_ccw(eigvecs)
    return {
        "n_member_cells": int(lon.size),
        "centroid_lon_unweighted": lon0,
        "centroid_lat_unweighted": lat0,
        "eigenvalue_major_km2": float(eigvals[0]),
        "eigenvalue_minor_km2": float(eigvals[1]),
        "major_axis_km": 2.0 * a,
        "minor_axis_km": 2.0 * b,
        "area_km2": math.pi * a * b,
        "axis_ratio": b / a,
        "orientation_east_ccw_deg": orientation,
        "azimuth_north_cw_deg": float(azimuth_from_north(orientation)),
        "major_axis_east": math.cos(math.radians(orientation)),
        "major_axis_north": math.sin(math.radians(orientation)),
        "_center_xy": center_xy, "_eigvecs": eigvecs, "_eigvals": eigvals,
    }


def points_inside_ellipse(lon, lat, center_lon, center_lat, center_xy, eigvecs,
                          eigvals, scale_factor, eigenvalue_floor_km2,
                          scale_factor_minor=None):
    """True for points inside the ellipse (boundary included, 1e-12 slack).

    A different scale factor can be given for the minor axis, as used in
    the Appendix B sensitivity analysis.
    """
    if scale_factor_minor is None:
        scale_factor_minor = scale_factor
    rel = project_lonlat_to_km(lon, lat, center_lon, center_lat) - center_xy
    pc = rel @ eigvecs
    a = scale_factor * math.sqrt(max(float(eigvals[0]), eigenvalue_floor_km2))
    b = scale_factor_minor * math.sqrt(max(float(eigvals[1]), eigenvalue_floor_km2))
    return ((pc[:, 0] / a) ** 2 + (pc[:, 1] / b) ** 2) <= 1.0 + 1e-12


def weighted_centroid(lon, lat, tmax_c):
    """Temperature-weighted centroid: sum(lon*T)/sum(T), sum(lat*T)/sum(T).

    T is the raw daily Tmax in degrees Celsius of the cluster's own member
    cells. Stops if any weight is missing, zero or negative (there is no
    silent fallback to the unweighted centroid).
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    w = np.asarray(tmax_c, dtype=float)
    if w.size == 0 or not np.all(np.isfinite(w)) or np.any(w <= 0.0):
        raise ValueError("weighted centroid needs finite, positive Tmax weights")
    wsum = float(np.sum(w))
    return float(np.sum(lon * w) / wsum), float(np.sum(lat * w) / wsum), wsum


# Stage 5: event classification, axial statistics and trend tests
TYPE_NAMES = {1: "Widespread (Isolated)", 2: "Spatially Clustered",
              3: "Temporally Clustered", 4: "Compound Clustering (Multi-Type)"}


def classify_event(daily_cluster_counts):
    """Type 1-4 from the number of clusters on each day of an event.

    Type 1: one day, one cluster.  Type 2: one day, several clusters.
    Type 3: several days, either one cluster every day or several clusters
    every day.  Type 4: several days mixing one-cluster and several-cluster
    days.
    """
    counts = list(daily_cluster_counts)
    if len(counts) <= 1:
        return 1 if counts and counts[0] == 1 else 2
    all_one = all(c == 1 for c in counts)
    all_multi = all(c > 1 for c in counts)
    return 3 if (all_one or all_multi) else 4


# Axial statistics and trend tests
def axial_mean_and_resultant(azimuths_deg):
    """Doubled-angle (axial) mean of orientations and its resultant length.

    mu = 0.5 * atan2(mean sin 2theta, mean cos 2theta), wrapped to (-90, 90];
    R = hypot(mean cos 2theta, mean sin 2theta). Axes 180 degrees apart are
    the same axis, so an ordinary arithmetic mean is not used.
    """
    arr = np.asarray(azimuths_deg, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    two = np.deg2rad(2.0 * arr)
    c = float(np.sum(np.cos(two)) / arr.size)
    s = float(np.sum(np.sin(two)) / arr.size)
    r = float(np.hypot(c, s))
    if r <= 1e-12:
        raise ValueError("axial mean is undefined for these orientations")
    mu = 0.5 * float(np.degrees(np.arctan2(s, c)))
    mu = float(np.mod(mu + 90.0, 180.0) - 90.0)
    return (mu + 180.0 if mu <= -90.0 else mu), r


def mann_kendall_test(y):
    """Two-sided Mann-Kendall test with tie-corrected variance and continuity
    correction."""
    from scipy.stats import norm
    y = np.asarray(y, dtype=float)
    n = len(y)
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(y[j] - y[k])
    _, counts = np.unique(y, return_counts=True)
    var_s = n * (n - 1) * (2 * n + 5)
    var_s = (var_s - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18.0
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    if p < 0.05 and z > 0:
        trend = "increasing"
    elif p < 0.05 and z < 0:
        trend = "decreasing"
    else:
        trend = "no significant trend"
    return s, z, p, trend


def sens_slope(x, y):
    """Median of all pairwise slopes."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slopes = []
    for i in range(len(x) - 1):
        for j in range(i + 1, len(x)):
            if x[j] != x[i]:
                slopes.append((y[j] - y[i]) / (x[j] - x[i]))
    return np.median(slopes)


# Geographic classification
def country_of_points(lon, lat, boundaries_path):
    """Name of the country (Natural Earth ADMIN field) containing each point.

    Uses the even-odd ray-casting test on every polygon ring, the rule used
    for the country counts reported in Section 4 (points outside every
    country, for example over the sea, get "Unassigned").
    """
    import shapefile
    points = np.column_stack((np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)))
    names = np.array(["Unassigned"] * len(points), dtype=object)
    reader = shapefile.Reader(str(boundaries_path))
    admin_field = [field[0] for field in reader.fields[1:]].index("ADMIN")
    for record in reader.iterShapeRecords():
        vertices = np.array(record.shape.points)
        parts = list(record.shape.parts) + [len(vertices)]
        inside = np.zeros(len(points), dtype=bool)
        for a, b in zip(parts[:-1], parts[1:]):
            inside ^= _points_in_ring(vertices[a:b], points)
        names[inside] = record.record[admin_field]
    return names


def _points_in_ring(ring, points):
    out = np.zeros(len(points), dtype=bool)
    if len(ring) < 3:
        return out
    xmin, ymin = ring.min(axis=0)
    xmax, ymax = ring.max(axis=0)
    subset = np.flatnonzero((points[:, 0] >= xmin) & (points[:, 0] <= xmax)
                            & (points[:, 1] >= ymin) & (points[:, 1] <= ymax))
    if len(subset) == 0:
        return out
    x, y = points[subset].T
    flags = np.zeros(len(subset), dtype=bool)
    xj, yj = ring[-1]
    for xi, yi in ring:
        if yi != yj:
            flags ^= ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        xj, yj = xi, yi
    out[subset] = flags
    return out
