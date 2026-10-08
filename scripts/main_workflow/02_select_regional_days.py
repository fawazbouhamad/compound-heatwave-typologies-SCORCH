"""Select days with regionally extensive heatwaves and their heatwave cells.

Reads step-1 heatwave flags; writes selected days and cells in Catalog 2.
The full daily summary is saved with the supporting results.
Run: python scripts/main_workflow/02_select_regional_days.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions as f  # noqa: E402


def read_flags(path):
    with Dataset(path) as ds:
        time = ds.variables["time"]
        dates = np.array([d.strftime("%Y-%m-%d")
                          for d in num2date(time[:], time.units, time.calendar)])
        out = {"dates": dates,
               "lat": np.asarray(ds.variables["lat"][:]),
               "lon": np.asarray(ds.variables["lon"][:])}
        for name in ("exceedance", "heatwave", "local_event_id"):
            ds.variables[name].set_auto_mask(False)
            out[name] = ds.variables[name][:]
    return out


def main():
    config = f.load_config()
    settings = config["regional_days"]
    support = f.supporting_results_dir(config) / "main_workflow"
    daily_dir = support / "regional_selection"
    daily_dir.mkdir(parents=True, exist_ok=True)
    out = f.catalog_dir(config, "02_regionally_extensive_heatwaves")
    flags_path = support / "heatwave_detection/heatwave_flags.nc"
    flags = read_flags(flags_path)
    grid = f.read_grid_cells(config)
    dates = flags["dates"]
    f.check_study_calendar(dates, config)
    n_time = len(dates)
    n_cells = flags["heatwave"].shape[1] * flags["heatwave"].shape[2]

    # Count the heatwave cells on each day, including days with zero cells.
    heatwave = flags["heatwave"].reshape(n_time, -1).astype(bool)
    exceed = flags["exceedance"].reshape(n_time, -1).astype(bool)
    n_heatwave = heatwave.sum(axis=1)
    n_exceed = exceed.sum(axis=1)
    n_both = (heatwave & exceed).sum(axis=1)

    threshold = f.regional_threshold(n_heatwave, settings["quantile"],
                                     settings["quantile_method"])
    selected = n_heatwave >= threshold

    daily = pd.DataFrame({
        "date": dates,
        "n_domain_cells": n_cells,
        "n_exceeding_cells": n_exceed,
        "n_heatwave_cells": n_heatwave,
        "n_heatwave_and_exceeding_cells": n_both,
        "heatwave_fraction_of_domain": n_heatwave / n_cells,
        # Share of exceeding cells that are also heatwave cells (Figure 2b);
        # left empty on days without any exceeding cell.
        "heatwave_share_of_exceeding_cells": np.where(
            n_exceed > 0, n_both / np.where(n_exceed > 0, n_exceed, 1), np.nan),
        "regional_threshold_cells": int(threshold),
        "selected": selected.astype(int),
    })

    selected_days = daily[daily["selected"] == 1][[
        "date", "n_heatwave_cells", "heatwave_fraction_of_domain",
        "n_exceeding_cells", "n_heatwave_and_exceeding_cells",
        "heatwave_share_of_exceeding_cells"]].copy()

    # Heatwave cells of every selected date, west to east and south to north
    # within each date (the cell order used by the clustering in step 3).
    lat_index, lon_index = np.meshgrid(np.arange(len(flags["lat"])),
                                       np.arange(len(flags["lon"])), indexing="ij")
    cell_id = grid.set_index(["lat_index", "lon_index"])["cell_id"]
    ids = cell_id.loc[list(zip(lat_index.ravel(), lon_index.ravel()))].to_numpy()
    cell_lat = flags["lat"][lat_index.ravel()]
    cell_lon = flags["lon"][lon_index.ravel()]
    local_ids = flags["local_event_id"].reshape(n_time, -1)
    rows = []
    for t in np.where(selected)[0]:
        on = np.where(heatwave[t])[0]
        on = on[f.longitude_major_order(cell_lon[on], cell_lat[on])]
        rows.append(pd.DataFrame({
            "date": dates[t],
            "cell_id": ids[on],
            "lat": cell_lat[on],
            "lon": cell_lon[on],
            "local_event_id": local_ids[t, on],
        }))
    cells = pd.concat(rows, ignore_index=True)

    daily_path = daily_dir / "daily_summary.csv"
    selected_path = out / "selected_days.csv"
    cells_path = out / "selected_day_cells.csv"
    daily.to_csv(daily_path, index=False, lineterminator="\n")
    selected_days.to_csv(selected_path, index=False, lineterminator="\n")
    cells.to_csv(cells_path, index=False, lineterminator="\n")

    # Consistency checks within this stage.
    per_day = cells.groupby("date").size().reindex(selected_days["date"]).to_numpy()
    if (per_day != selected_days["n_heatwave_cells"].to_numpy()).any():
        raise SystemExit("selected-day cell lists do not match the daily counts")
    if (cells["local_event_id"] <= 0).any():
        raise SystemExit("a selected heatwave cell has no local event")

    print(f"regional threshold: {int(threshold)} heatwave cells "
          f"({threshold / n_cells:.2%} of the domain)")
    print(f"selected days: {int(selected.sum())}")
    f.write_run_record(config, __file__, [flags_path, f.supporting_input_dir(config) / "maps" / "grid_cells.csv"],
                       [daily_path, selected_path, cells_path],
                       notes={"regional_threshold_cells": int(threshold),
                              "selected_days": int(selected.sum())})


if __name__ == "__main__":
    main()
