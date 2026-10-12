# SPDX-License-Identifier: GPL-3.0-only
"""Shared helpers for the SCORCH tests."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import config as cfg

SEASON_2000 = pd.date_range("2000-04-01", "2000-09-30", freq="D")


def exceedance_series(index, ones_at):
    """A 0/1 series over ``index`` with ones at the given positions."""
    values = np.zeros(len(index), dtype=int)
    values[list(ones_at)] = 1
    return pd.Series(values, index=index)


@pytest.fixture
def manuscript_config(tmp_path):
    """Manuscript settings with every folder pointed at a temporary directory."""
    return cfg.load_config(overrides={"paths": {
        "base_dir": str(tmp_path), "input_dir": "input", "supporting_input_dir": "input",
        "supporting_results_dir": "output/supporting_results", "output_dir": "output",
        "run_records_dir": "output/run_records"}}, environment={})


def write_small_netcdf(path, dates, lat, lon, tmax, units="degC", time_units=None,
                       drop_variable=None, dims=("time", "lat", "lon")):
    """Write a small processed-temperature NetCDF file for the I/O tests."""
    from netCDF4 import Dataset, date2num
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", len(dates))
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        if drop_variable != "time":
            tv = ds.createVariable("time", "i4", ("time",))
            tv.units = time_units or "days since 2000-01-01 00:00:00"
            tv.calendar = "standard"
            tv[:] = date2num(list(pd.to_datetime(dates).to_pydatetime()), units=tv.units,
                             calendar=tv.calendar)
        for name, values in (("lat", lat), ("lon", lon)):
            if drop_variable == name:
                continue
            v = ds.createVariable(name, "f8", (name,))
            v[:] = values
        if drop_variable != "tmax":
            v = ds.createVariable("tmax", "f8", dims)
            if units is not None:
                v.units = units
            v[:] = np.transpose(tmax, [("time", "lat", "lon").index(d) for d in dims])
    return path
