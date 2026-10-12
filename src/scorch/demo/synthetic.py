# SPDX-License-Identifier: GPL-3.0-only
"""A small synthetic temperature dataset with known compound heatwave events.

The demo grid is a 20 x 16 block of one-degree cells (320 cells) and the
demo period is April-September of four years (732 dates). The temperatures
are invented: a smooth seasonal cycle, cooler cells to the north, and mild
random weather noise. On top of that, warm "blobs" of cells are made hot
(+10 degC) on chosen days so that SCORCH finds one example of each event
type and two longer examples:

* Type 1: a single hot day with one connected blob
* Type 2: a single hot day with two separate blobs
* Type 3: four consecutive days with one blob each day, and three days with
  two blobs each day
* Type 4: multi-day events that start with one blob and gain a second one

None of this is real climate data. The coordinates are arbitrary and the
figure produced by ``scorch demo`` is labelled as synthetic.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .. import io

DEMO_SEED = 20261011
LAT = np.arange(30.5, 46.0, 1.0)          # 16 rows of cells, 30.5 to 45.5
LON = np.arange(20.5, 40.0, 1.0)          # 20 columns of cells, 20.5 to 39.5
FIRST_YEAR, LAST_YEAR = 1991, 1994
WARM_SEASON_MONTHS = [4, 5, 6, 7, 8, 9]
ANOMALY_C = 10.0                          # added to every blob cell on its hot days
# Background climate: a gentle seasonal cycle and mild, mostly local weather noise,
# so that the random weather alone rarely puts many cells in a heatwave at once.
SEASONAL_AMPLITUDE_C = 2.0
REGIONAL_NOISE = {"phi": 0.6, "marginal_sd": 0.4}
LOCAL_NOISE = {"phi": 0.5, "marginal_sd": 1.5}


@dataclass(frozen=True)
class Blob:
    """A disc of grid cells made hot from ``start`` for ``days`` days."""
    lon: float
    lat: float
    radius: float
    start: str
    days: int

    def dates(self):
        return [d.strftime("%Y-%m-%d")
                for d in pd.date_range(self.start, periods=self.days, freq="D")]

    def mask(self, lat=LAT, lon=LON):
        grid_lat, grid_lon = np.meshgrid(lat, lon, indexing="ij")
        return (grid_lon - self.lon) ** 2 + (grid_lat - self.lat) ** 2 <= self.radius ** 2


@dataclass(frozen=True)
class DesignedEvent:
    name: str
    expected_type: int
    start: str
    end: str
    blobs: tuple
    explanation: str

    @property
    def expected_days(self):
        return [d.strftime("%Y-%m-%d")
                for d in pd.date_range(self.start, self.end, freq="D")]


# Blobs that start before (or end after) an event's selected days keep the
# day count high only on the selected days: the day before and the day after
# have one blob of 29 cells, below the regional threshold.
DESIGNED_EVENTS = (
    DesignedEvent(
        "Type 1 example", 1, "1991-07-10", "1991-07-10",
        (Blob(25.5, 36.5, 3, "1991-07-08", 3), Blob(30.5, 36.5, 3, "1991-07-10", 3)),
        "Two overlapping blobs are both hot on 10 July only, forming one large "
        "connected cluster on that single day."),
    DesignedEvent(
        "Type 3 example (one cluster per day)", 3, "1991-08-12", "1991-08-15",
        (Blob(29.5, 38.5, 4, "1991-08-12", 4),),
        "One large blob stays hot for four consecutive days: one cluster every day."),
    DesignedEvent(
        "Type 2 example", 2, "1992-06-20", "1992-06-20",
        (Blob(24.5, 34.5, 3, "1992-06-18", 3), Blob(35.5, 42.5, 3, "1992-06-20", 3)),
        "Two blobs far apart are both hot on 20 June only: two separate clusters "
        "on a single day."),
    DesignedEvent(
        "Type 4 example", 4, "1992-08-03", "1992-08-06",
        (Blob(33.5, 36.5, 4, "1992-08-03", 4), Blob(23.5, 43.5, 3, "1992-08-04", 3)),
        "A large blob is hot for four days; a distant second blob joins on days 2-4, "
        "so the event mixes one-cluster and two-cluster days."),
    DesignedEvent(
        "Type 3 example (two clusters per day)", 3, "1993-07-15", "1993-07-17",
        (Blob(24.5, 33.5, 3, "1993-07-15", 3), Blob(35.5, 43.5, 3, "1993-07-15", 3)),
        "Two distant blobs are hot on the same three days: two clusters every day."),
    DesignedEvent(
        "Type 4 example (longer)", 4, "1994-07-20", "1994-07-25",
        (Blob(26.5, 35.5, 4, "1994-07-20", 6), Blob(36.5, 43.5, 3, "1994-07-22", 3)),
        "A large blob is hot for six days; a distant second blob appears on days 3-5."),
)


def demo_settings():
    """Settings that differ from the manuscript for the synthetic grid."""
    return {
        "study_area": {
            "lat_min": float(LAT[0] - 0.5), "lat_max": float(LAT[-1] + 0.5),
            "lon_min": float(LON[0] - 0.5), "lon_max": float(LON[-1] + 0.5),
            "grid_step_deg": 1.0,
            "warm_season_months": list(WARM_SEASON_MONTHS),
            "first_year": FIRST_YEAR, "last_year": LAST_YEAR,
            "expected_dates": None, "expected_cells": int(LAT.size * LON.size),
        },
        # The Section 4 statistics describe the Middle East; they do not apply here.
        "summaries": {"in_text_statistics": False},
    }


def study_dates():
    days = pd.date_range(f"{FIRST_YEAR}-01-01", f"{LAST_YEAR}-12-31", freq="D")
    return days[days.month.isin(WARM_SEASON_MONTHS)]


def _ar1(rng, n, phi, marginal_sd, shape=()):
    """An AR(1) series of length n with the given marginal standard deviation."""
    innovation_sd = marginal_sd * np.sqrt(1.0 - phi ** 2)
    series = np.empty((n, *shape))
    series[0] = rng.normal(0.0, marginal_sd, shape)
    for t in range(1, n):
        series[t] = phi * series[t - 1] + rng.normal(0.0, innovation_sd, shape)
    return series


def make_synthetic_tmax(seed=DEMO_SEED, include_events=True):
    """Daily maximum temperatures (date, lat, lon) in degC for the demo grid.

    ``include_events=False`` gives the background weather alone (used to check
    that the designed events, not the noise, form the regional days).
    """
    rng = np.random.default_rng(seed)
    dates = study_dates()
    n_lat, n_lon = LAT.size, LON.size
    tmax = np.empty((len(dates), n_lat, n_lon))
    # Position of each date within its season (0 on 1 April, 182 on 30 September).
    season_start = pd.to_datetime(dates.year.astype(str) + "-04-01")
    day_of_season = (dates - season_start).days.to_numpy()
    seasonal = 30.0 + SEASONAL_AMPLITUDE_C * np.sin(np.pi * day_of_season / 182.0)
    north_cooling = -0.5 * (LAT - LAT[0])
    base = seasonal[:, None, None] + north_cooling[None, :, None] + np.zeros((1, 1, n_lon))
    # Weather noise is generated season by season (the seasons are not consecutive).
    for year in range(FIRST_YEAR, LAST_YEAR + 1):
        rows = np.where(dates.year == year)[0]
        regional = _ar1(rng, len(rows), shape=(), **REGIONAL_NOISE)
        local = _ar1(rng, len(rows), shape=(n_lat, n_lon), **LOCAL_NOISE)
        tmax[rows] = base[rows] + regional[:, None, None] + local
    # The designed events.
    date_index = {d.strftime("%Y-%m-%d"): k for k, d in enumerate(dates)}
    for event in (DESIGNED_EVENTS if include_events else ()):
        for blob in event.blobs:
            mask = blob.mask()
            for day in blob.dates():
                tmax[date_index[day]][mask] += ANOMALY_C
    return dates, LAT.copy(), LON.copy(), np.round(tmax, 3)


def write_netcdf(path, dates, lat, lon, tmax, seed):
    """Write the synthetic temperatures in the layout of the study input."""
    from netCDF4 import Dataset, date2num
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", len(dates))
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        tv = ds.createVariable("time", "i4", ("time",))
        tv.units = f"days since {FIRST_YEAR}-01-01 00:00:00"
        tv.calendar = "standard"
        tv.standard_name = "time"
        tv.long_name = f"date (April-September warm-season days, {FIRST_YEAR}-{LAST_YEAR})"
        tv.axis = "T"
        tv[:] = date2num(list(dates.to_pydatetime()), units=tv.units, calendar=tv.calendar)
        for name, values, units, standard in (("lat", lat, "degrees_north", "latitude"),
                                              ("lon", lon, "degrees_east", "longitude")):
            v = ds.createVariable(name, "f8", (name,))
            v.units = units
            v.standard_name = standard
            v.long_name = f"{standard} of grid-cell centre (synthetic grid)"
            v[:] = values
        v = ds.createVariable("tmax", "f8", ("time", "lat", "lon"), zlib=True, complevel=4,
                              chunksizes=(min(183, len(dates)), len(lat), len(lon)))
        v.units = "degC"
        v.standard_name = "air_temperature"
        v.long_name = "daily maximum 2-m air temperature, synthetic demo data"
        v[:] = tmax
        ds.title = "SCORCH demo: synthetic daily maximum temperature (not observations)"
        ds.comment = ("Invented temperatures for trying the SCORCH software: a seasonal cycle, "
                      "a north-south gradient, random weather noise and designed hot "
                      "events. Generated by scorch.demo.synthetic.make_synthetic_tmax "
                      f"with seed {seed}. The coordinates do not describe any real place.")
        ds.Conventions = "CF-1.10"
    return path


def write_inputs(input_dir, seed=DEMO_SEED):
    """Write tmax_processed.nc and maps/grid_cells.csv for the demo; return their paths."""
    input_dir = Path(input_dir)
    dates, lat, lon, tmax = make_synthetic_tmax(seed)
    netcdf = write_netcdf(input_dir / "tmax_processed.nc", dates, lat, lon, tmax, seed)
    grid = io.grid_cells_table(lat, lon)
    grid_path = input_dir / "maps" / "grid_cells.csv"
    io.write_csv(grid, grid_path)
    return {"netcdf": netcdf, "grid_cells": grid_path}


def compare_with_design(events):
    """One row per designed event with the type SCORCH found (or 'not found').

    Also returns the events that were not designed and the set of designed days.
    """
    rows = []
    for design in DESIGNED_EVENTS:
        found = events[(events["start_date"] == design.start) & (events["end_date"] == design.end)]
        if len(found) == 1:
            row = found.iloc[0]
            detected, counts, event_id = int(row["event_type"]), str(row["daily_cluster_counts"]), int(row["compound_event_id"])
        else:
            detected, counts, event_id = None, "", None
        rows.append({
            "designed_event": design.name,
            "dates": design.start if design.start == design.end else f"{design.start} to {design.end}",
            "expected_type": design.expected_type,
            "detected_type": detected,
            "compound_event_id": event_id,
            "daily_cluster_counts": counts,
            "matches_design": detected == design.expected_type,
        })
    comparison = pd.DataFrame(rows)
    designed_days = {d for e in DESIGNED_EVENTS for d in e.expected_days}
    extra = events[~events["start_date"].isin({e.start for e in DESIGNED_EVENTS})]
    return comparison, extra, designed_days
