# SPDX-License-Identifier: GPL-3.0-only
"""Reading inputs and writing results.

The processed temperature input is one NetCDF file (``tmax_processed.nc``)
with dimensions ``(time, lat, lon)``, a ``tmax`` variable in degrees Celsius
and a CF ``time`` coordinate holding every warm-season date of the study.
``docs/io_reference.md`` describes the file and every output table.
"""
from __future__ import annotations

import hashlib
import io
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as cfg
from . import __version__

GRID_CELL_COLUMNS = ["cell_id", "lat", "lon", "lat_index", "lon_index"]


class InputError(ValueError):
    """An input file is missing or does not have the expected content."""


# Checksums, versions and run records
def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def software_versions():
    versions = {"python": platform.python_version(), "scorch": __version__}
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


def write_run_record(config, script_name, inputs, outputs, notes=None):
    """Record what a stage read and wrote, with checksums and versions.

    ``script_name`` is a script path or a stage name; the record is saved as
    ``<run_records_dir>/<stem>.json``.
    """
    record = {
        "script": str(script_name),
        "finished_utc": now_utc(),
        "software": software_versions(),
        "inputs": {cfg.relative_to_repo(p, config): sha256_of(p) for p in inputs},
        "outputs": {cfg.relative_to_repo(p, config): sha256_of(p) for p in outputs},
    }
    if notes:
        record["notes"] = notes
    folder = cfg.run_records_dir(config)
    path = folder / f"{Path(script_name).stem}.json"
    write_json(path, record)
    return path


def write_csv(table, path, index=False):
    """Write a table as CSV with LF line endings (the convention of every catalog)."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=index, lineterminator="\n")
    return Path(path)


# Processed temperature input
def read_netcdf_tmax(path):
    """Read a processed temperature NetCDF file.

    Returns a dict with
      dates : array of "YYYY-MM-DD" strings, one per time step
      lat, lon : cell-centre coordinates in degrees (ascending)
      tmax : float64 array (date, lat, lon), degrees Celsius
      path : the file that was read
    """
    try:
        from netCDF4 import Dataset, num2date
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise InputError("the netCDF4 package is not installed; install SCORCH with "
                         "'pip install -e .' (see docs/installation.md)") from error

    path = Path(path)
    if not path.is_file():
        raise InputError(
            f"input file not found: {path}\n"
            "Expected the processed temperature file tmax_processed.nc. For the study, "
            "download SCORCH_data.zip from https://doi.org/10.5281/zenodo.21717874, "
            "extract it and copy scorch_data/tmax_processed.nc into the input folder "
            "(or pass --input-dir / set SCORCH_INPUT_DIR). To try the software "
            "without the dataset, run 'scorch demo'.")
    with Dataset(path) as ds:
        for name in ("time", "lat", "lon", "tmax"):
            if name not in ds.variables:
                raise InputError(f"{path.name} has no variable '{name}'. It needs the "
                                 "variables time, lat, lon and tmax(time, lat, lon); "
                                 "see docs/io_reference.md")
        variable = ds.variables["tmax"]
        if tuple(variable.dimensions) != ("time", "lat", "lon"):
            raise InputError(f"tmax must have dimensions (time, lat, lon); found "
                             f"{tuple(variable.dimensions)}")
        time = ds.variables["time"]
        units = getattr(time, "units", None)
        calendar = getattr(time, "calendar", "standard")
        if not units:
            raise InputError("the time variable needs CF units such as "
                             "'days since 1940-01-01 00:00:00'")
        dates = np.array([d.strftime("%Y-%m-%d")
                          for d in num2date(time[:], units, calendar)])
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        variable.set_auto_mask(False)
        tmax = np.asarray(variable[:], dtype=np.float64)
        tmax_units = getattr(variable, "units", None)
    if tmax_units != "degC":
        raise InputError(f"tmax units are {tmax_units!r}, expected 'degC' (degrees Celsius)")
    if not np.all(np.isfinite(tmax)):
        raise InputError("tmax contains missing or non-finite values; every cell needs a "
                         "value on every date (SCORCH does not fill gaps)")
    for name, values in (("lat", lat), ("lon", lon)):
        if values.ndim != 1 or len(values) < 2 or not np.all(np.diff(values) > 0):
            raise InputError(f"{name} must be a one-dimensional, strictly increasing coordinate")
    if len(np.unique(dates)) != len(dates):
        raise InputError("the time coordinate contains repeated dates")
    return {"dates": dates, "lat": lat, "lon": lon, "tmax": tmax, "path": path}


def read_tmax_input(config):
    """Read the processed temperature file named by the settings."""
    return read_netcdf_tmax(cfg.input_file(config))


def describe_netcdf(path):
    """A short description of a NetCDF input (for ``scorch check-input``)."""
    data = read_netcdf_tmax(path)
    tmax = data["tmax"]
    return {
        "file": str(data["path"]),
        "sha256": sha256_of(data["path"]),
        "dates": int(len(data["dates"])),
        "first_date": str(data["dates"][0]),
        "last_date": str(data["dates"][-1]),
        "lat_cells": int(len(data["lat"])),
        "lon_cells": int(len(data["lon"])),
        "cells": int(len(data["lat"]) * len(data["lon"])),
        "lat_range": [float(data["lat"][0]), float(data["lat"][-1])],
        "lon_range": [float(data["lon"][0]), float(data["lon"][-1])],
        "tmax_min_c": float(tmax.min()),
        "tmax_max_c": float(tmax.max()),
    }


def centroid_temperature_weights(values):
    """Prepare the single-precision temperature weights used for centroids.

    Convert to decimal text, read with pandas' default parser, then round to
    float32. Keeping this exact sequence preserves the manuscript's weights.
    """
    text = "v\n" + "\n".join(repr(float(v)) for v in np.ravel(values))
    parsed = pd.read_csv(io.StringIO(text))["v"].to_numpy()
    return parsed.astype(np.float32).reshape(np.shape(values))


# Grid cells
def read_grid_cells(config):
    """The grid-cell table (cell_id, lat, lon, lat_index, lon_index, ...)."""
    path = cfg.supporting_input_dir(config) / "maps" / "grid_cells.csv"
    if not path.is_file():
        raise InputError(f"grid-cell table not found: {path}\n"
                         "It lists every grid cell with cell_id, lat, lon, lat_index and "
                         "lon_index. For the study it is part of the repository; for your own "
                         "grid, write one with scorch.io.grid_cells_table() or see "
                         "docs/io_reference.md")
    grid = pd.read_csv(path)
    missing = [c for c in GRID_CELL_COLUMNS if c not in grid.columns]
    if missing:
        raise InputError(f"{path} lacks the columns {missing}")
    return grid


def grid_cells_table(lat, lon, grid_step_deg=1.0):
    """Build the grid-cell table for ascending cell-centre coordinates.

    Cells are numbered 1, 2, ... by latitude row and then longitude, the
    numbering used by the study grid.
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat_index, lon_index = np.meshgrid(np.arange(len(lat)), np.arange(len(lon)), indexing="ij")
    half = float(grid_step_deg) / 2.0
    table = pd.DataFrame({
        "cell_id": np.arange(1, lat.size * lon.size + 1),
        "lat": lat[lat_index.ravel()],
        "lon": lon[lon_index.ravel()],
        "lat_index": lat_index.ravel(),
        "lon_index": lon_index.ravel(),
    })
    table["lat_south"] = table["lat"] - half
    table["lat_north"] = table["lat"] + half
    table["lon_west"] = table["lon"] - half
    table["lon_east"] = table["lon"] + half
    return table


def check_grid_matches(grid, lat, lon):
    """Stop unless the grid table describes exactly the input coordinates."""
    expected = grid_cells_table(lat, lon)
    same = (len(grid) == len(expected)
            and np.array_equal(grid["lat_index"].to_numpy(), expected["lat_index"].to_numpy())
            and np.array_equal(grid["lon_index"].to_numpy(), expected["lon_index"].to_numpy())
            and np.allclose(grid["lat"].to_numpy(float), expected["lat"].to_numpy(), atol=1e-9)
            and np.allclose(grid["lon"].to_numpy(float), expected["lon"].to_numpy(), atol=1e-9))
    if not same:
        raise InputError("grid_cells.csv does not match the latitude/longitude coordinates of "
                         "the input file (different cells, order or coordinates)")


# Study calendar
def expected_study_dates(config):
    """The warm-season days of the study years as "YYYY-MM-DD" strings."""
    area = config["study_area"]
    days = pd.date_range(f"{int(area['first_year'])}-01-01",
                         f"{int(area['last_year'])}-12-31", freq="D")
    days = days[days.month.isin([int(m) for m in area["warm_season_months"]])]
    return days.strftime("%Y-%m-%d").to_numpy()


def check_study_calendar(dates, config):
    """Stop unless the dates are exactly the expected warm-season calendar."""
    expected = expected_study_dates(config)
    expected_count = config["study_area"].get("expected_dates")
    if expected_count is not None and len(expected) != int(expected_count):
        raise InputError(f"the settings expect {expected_count} dates but the years "
                         f"{config['study_area']['first_year']}-{config['study_area']['last_year']} "
                         f"and months {config['study_area']['warm_season_months']} give "
                         f"{len(expected)}; correct study_area in the settings")
    dates = np.asarray(dates)
    if len(dates) != len(expected) or not np.array_equal(dates, expected):
        first, last = (dates[0], dates[-1]) if len(dates) else ("none", "none")
        raise InputError(
            "the input dates are not the expected warm-season calendar: expected every day "
            f"of months {config['study_area']['warm_season_months']} from "
            f"{config['study_area']['first_year']} to {config['study_area']['last_year']} "
            f"({len(expected)} dates); the file has {len(dates)} dates from {first} to {last}. "
            "Adjust study_area in the settings file to match your data.")
