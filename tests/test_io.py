# SPDX-License-Identifier: GPL-3.0-only
"""Input reading: helpful errors, calendar checks and the grid table."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import config as cfg, io
from conftest import write_small_netcdf

DATES = pd.date_range("2000-04-01", "2000-09-30", freq="D")
LAT = np.array([10.5, 11.5])
LON = np.array([20.5, 21.5, 22.5])


def small_tmax():
    return np.random.default_rng(1).normal(30, 2, size=(len(DATES), 2, 3))


def test_read_valid_file(tmp_path):
    path = write_small_netcdf(tmp_path / "t.nc", DATES, LAT, LON, small_tmax())
    data = io.read_netcdf_tmax(path)
    assert data["tmax"].shape == (183, 2, 3) and data["tmax"].dtype == np.float64
    assert data["dates"][0] == "2000-04-01" and data["dates"][-1] == "2000-09-30"
    description = io.describe_netcdf(path)
    assert description["cells"] == 6 and description["dates"] == 183


def test_missing_file_points_to_the_dataset_record(tmp_path):
    with pytest.raises(io.InputError, match="10.5281/zenodo.21717874"):
        io.read_netcdf_tmax(tmp_path / "nothing.nc")


def test_wrong_units_are_rejected(tmp_path):
    path = write_small_netcdf(tmp_path / "k.nc", DATES, LAT, LON, small_tmax() + 273.15, units="K")
    with pytest.raises(io.InputError, match="expected 'degC'"):
        io.read_netcdf_tmax(path)


def test_missing_values_are_rejected(tmp_path):
    tmax = small_tmax()
    tmax[3, 1, 1] = np.nan
    path = write_small_netcdf(tmp_path / "nan.nc", DATES, LAT, LON, tmax)
    with pytest.raises(io.InputError, match="missing or non-finite"):
        io.read_netcdf_tmax(path)


def test_missing_variable_and_wrong_dimension_order(tmp_path):
    path = write_small_netcdf(tmp_path / "novar.nc", DATES, LAT, LON, small_tmax(), drop_variable="tmax")
    with pytest.raises(io.InputError, match="no variable 'tmax'"):
        io.read_netcdf_tmax(path)
    path = write_small_netcdf(tmp_path / "dims.nc", DATES, LAT, LON, small_tmax(),
                              dims=("lat", "lon", "time"))
    with pytest.raises(io.InputError, match=r"dimensions \(time, lat, lon\)"):
        io.read_netcdf_tmax(path)


def test_study_calendar_check():
    config = cfg.load_config(overrides={"study_area": {"first_year": 2000, "last_year": 2000,
                                                       "expected_dates": None, "expected_cells": None}},
                             environment={})
    dates = np.array([d.strftime("%Y-%m-%d") for d in DATES])
    io.check_study_calendar(dates, config)                     # the right calendar passes
    with pytest.raises(io.InputError, match="Adjust study_area"):
        io.check_study_calendar(dates[:-1], config)
    manuscript = cfg.load_config(environment={})
    assert len(io.expected_study_dates(manuscript)) == 15738
    bad = cfg.load_config(overrides={"study_area": {"expected_dates": 100}}, environment={})
    with pytest.raises(io.InputError, match="expect 100 dates"):
        io.check_study_calendar(dates, bad)


def test_grid_table_and_match_check():
    grid = io.grid_cells_table(LAT, LON)
    assert grid["cell_id"].tolist() == [1, 2, 3, 4, 5, 6]
    assert grid.loc[1, ["lat", "lon", "lat_index", "lon_index"]].tolist() == [10.5, 21.5, 0, 1]
    assert grid.loc[5, ["lat_south", "lat_north", "lon_west", "lon_east"]].tolist() == [11.0, 12.0, 22.0, 23.0]
    io.check_grid_matches(grid, LAT, LON)
    with pytest.raises(io.InputError, match="does not match"):
        io.check_grid_matches(grid, LAT, LON + 1.0)


def test_centroid_temperature_weights_round_to_float32():
    values = np.array([30.123456789, 41.9876543219])
    weights = io.centroid_temperature_weights(values)
    assert weights.dtype == np.float32
    assert np.allclose(weights, values.astype(np.float32))
