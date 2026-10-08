"""Find local heatwaves using each cell's warm-season 95th percentile.

Reads input/tmax_processed.nc and supporting_files/additional_inputs/maps/grid_cells.csv.
Writes the local-event catalog; saves flags and thresholds as supporting results.
Run: python scripts/main_workflow/01_detect_heatwaves.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions as f  # noqa: E402


def detect_all_cells(tmax, thresholds, dates, settings):
    """Exceedance, heatwave membership and heatwave summaries for all cells.

    Returns
      exceed     int8 array (date, lat, lon): 1 when Tmax >= threshold
      number     int32 array (date, lat, lon): heatwave number within the
                 cell (0 = not in a heatwave), counted 1, 2, 3, ... through
                 the record in each cell
      episodes   list of per-heatwave rows
    """
    n_time, n_lat, n_lon = tmax.shape
    # A value equal to the threshold counts as exceedance.
    exceed = (tmax >= thresholds[None, :, :]).astype(np.int8)
    number = np.zeros((n_time, n_lat, n_lon), dtype=np.int32)
    episodes = []
    for i in range(n_lat):
        for j in range(n_lon):
            labels, cell_episodes = f.label_cell_by_season(
                exceed[:, i, j], dates,
                min_len=settings["min_duration_days"],
                min_ones=settings["min_exceeding_days"])
            number[:, i, j] = labels
            for hw_number, start, end, length, n_exceed, n_bridged in cell_episodes:
                episodes.append((i, j, hw_number, str(start), str(end),
                                 length, n_exceed, n_bridged))
    return exceed, number, episodes


def summarize_events(episodes, tmax, dates, lat, lon, grid):
    """One row per local heatwave, with its peak date and peak Tmax.

    The peak is the date of the highest Tmax during the heatwave (bridged
    days included). If two days share the highest value, the earlier date is
    taken.
    """
    date_index = {d: k for k, d in enumerate(dates)}
    cell_id = grid.set_index(["lat_index", "lon_index"])["cell_id"]
    rows = []
    for i, j, hw_number, start, end, length, n_exceed, n_bridged in episodes:
        t0, t1 = date_index[start], date_index[end]
        series = tmax[t0:t1 + 1, i, j]
        k = int(np.argmax(series))          # first occurrence of the maximum
        rows.append({
            "cell_id": int(cell_id.loc[(i, j)]),
            "lat": lat[i],
            "lon": lon[j],
            "event_number_in_cell": hw_number,
            "start_date": start,
            "end_date": end,
            "duration_days": length,
            "exceeding_days": n_exceed,
            "bridged_days": n_bridged,
            "peak_date": dates[t0 + k],
            "peak_tmax_c": float(series[k]),
        })
    events = pd.DataFrame(rows).sort_values(["cell_id", "start_date"])
    events.insert(0, "local_event_id", np.arange(1, len(events) + 1))
    return events.reset_index(drop=True)


def event_id_array(number, events, grid):
    """Replace within-cell heatwave numbers by the catalog-wide local_event_id."""
    ids = np.zeros(number.shape, dtype=np.int32)
    lookup = events.merge(grid[["cell_id", "lat_index", "lon_index"]], on="cell_id")
    for (i, j), group in lookup.groupby(["lat_index", "lon_index"]):
        mapping = np.zeros(group["event_number_in_cell"].max() + 1, dtype=np.int32)
        mapping[group["event_number_in_cell"].to_numpy()] = group["local_event_id"].to_numpy()
        ids[:, i, j] = mapping[number[:, i, j]]
    return ids


def write_flags(path, dates, lat, lon, exceed, heatwave, local_event_id):
    tmp = path.with_suffix(".nc.tmp")
    n_time, n_lat, n_lon = exceed.shape
    chunks = (183, n_lat, n_lon)
    with Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("time", n_time)
        ds.createDimension("lat", n_lat)
        ds.createDimension("lon", n_lon)
        tv = ds.createVariable("time", "i4", ("time",))
        tv.units = "days since 1940-01-01 00:00:00"
        tv.calendar = "standard"
        tv.long_name = "date (April-September warm-season days, 1940-2025)"
        tv[:] = date2num(list(pd.to_datetime(dates).to_pydatetime()),
                         units=tv.units, calendar=tv.calendar)
        for name, values, units in (("lat", lat, "degrees_north"),
                                    ("lon", lon, "degrees_east")):
            v = ds.createVariable(name, "f8", (name,))
            v.units = units
            v.long_name = ("latitude" if name == "lat" else "longitude") + " of grid-cell centre"
            v[:] = values

        v = ds.createVariable("exceedance", "i1", ("time", "lat", "lon"),
                              zlib=True, complevel=4, chunksizes=chunks)
        v.long_name = "1 when daily Tmax >= the cell's 95th-percentile threshold, else 0"
        v.flag_values = np.array([0, 1], dtype=np.int8)
        v.flag_meanings = "below_threshold at_or_above_threshold"
        v[:] = exceed

        v = ds.createVariable("heatwave", "i1", ("time", "lat", "lon"),
                              zlib=True, complevel=4, chunksizes=chunks)
        v.long_name = ("1 when the cell is part of a local heatwave on that date "
                       "(including bridged single non-exceeding days), else 0")
        v.flag_values = np.array([0, 1], dtype=np.int8)
        v.flag_meanings = "not_in_heatwave in_heatwave"
        v[:] = heatwave

        v = ds.createVariable("local_event_id", "i4", ("time", "lat", "lon"),
                              zlib=True, complevel=4, chunksizes=chunks)
        v.long_name = "local_event_id of heatwave_events.csv (0 = not in a heatwave)"
        v[:] = local_event_id

        ds.title = "SCORCH Catalog 1: daily exceedance and local heatwave flags"
        ds.comment = ("Every date and cell of the study is present; there are no "
                      "missing values. Exceedance and heatwave membership differ: "
                      "a bridged day has exceedance 0 and heatwave 1, and an "
                      "exceedance that never forms a heatwave has exceedance 1 "
                      "and heatwave 0.")
        ds.history = "written by scripts/main_workflow/01_detect_heatwaves.py"
    tmp.replace(path)


def main():
    started = time.time()
    config = f.load_config()
    settings = config["local_heatwaves"]
    out = f.catalog_dir(config, "01_heatwave_events")
    support = f.supporting_results_dir(config) / "main_workflow/heatwave_detection"
    support.mkdir(parents=True, exist_ok=True)

    data = f.read_tmax_input(config)
    grid = f.read_grid_cells(config)
    dates, lat, lon, tmax = data["dates"], data["lat"], data["lon"], data["tmax"]
    f.check_study_calendar(dates, config)
    if tmax.shape[1] * tmax.shape[2] != config["study_area"]["expected_cells"]:
        raise SystemExit("unexpected number of grid cells")

    print("calculating cell thresholds")
    thresholds = f.cell_thresholds(tmax, settings["threshold_quantile"])

    print("detecting local heatwaves in 1,800 cells")
    exceed, number, episodes = detect_all_cells(tmax, thresholds, dates, settings)
    events = summarize_events(episodes, tmax, dates, lat, lon, grid)
    local_event_id = event_id_array(number, events, grid)
    heatwave = (local_event_id > 0).astype(np.int8)

    threshold_table = grid[["cell_id", "lat", "lon"]].copy()
    threshold_table["threshold_tmax_c"] = thresholds.ravel()
    threshold_table["n_exceeding_days"] = exceed.reshape(len(dates), -1).sum(axis=0)
    threshold_table["n_heatwave_days"] = heatwave.reshape(len(dates), -1).sum(axis=0)
    threshold_table["n_local_events"] = (events.groupby("cell_id").size()
                                         .reindex(threshold_table["cell_id"], fill_value=0)
                                         .to_numpy())

    flags_path = support / "heatwave_flags.nc"
    events_path = out / "heatwave_events.csv"
    thresholds_path = support / "cell_thresholds.csv"
    write_flags(flags_path, dates, lat, lon, exceed, heatwave, local_event_id)
    events.to_csv(events_path, index=False, lineterminator="\n")
    threshold_table.to_csv(thresholds_path, index=False, lineterminator="\n")

    # Consistency checks within this stage.
    if not np.array_equal(heatwave, (number > 0).astype(np.int8)):
        raise SystemExit("heatwave membership and event numbers disagree")
    if int(events["duration_days"].sum()) != int(heatwave.sum()):
        raise SystemExit("event durations do not add up to the heatwave cell-days")
    if (events["start_date"].str[:4] != events["end_date"].str[:4]).any():
        raise SystemExit("a heatwave crosses the gap between warm seasons")

    print(f"local heatwaves: {len(events):,}")
    print(f"heatwave cell-days: {int(heatwave.sum()):,}; "
          f"exceedance cell-days: {int(exceed.sum()):,}")
    f.write_run_record(config, __file__,
                       [data["path"], f.supporting_input_dir(config) / "maps" / "grid_cells.csv"],
                       [flags_path, events_path, thresholds_path],
                       notes={"local_events": int(len(events)),
                              "runtime_seconds": round(time.time() - started, 1)})


if __name__ == "__main__":
    main()
