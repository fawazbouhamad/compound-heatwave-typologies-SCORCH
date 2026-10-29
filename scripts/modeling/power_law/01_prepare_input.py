"""Prepare the daily and event maximum ellipse areas.

Reads: Ellipse and compound-event catalogs.
Writes: supporting_files/additional_results/modeling/power_law/
Run: python scripts/modeling/power_law/01_prepare_input.py
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402

ELLIPSES_CSV = REPO_DIR / "output/catalogs/04_extreme_heatwave_ellipses/ellipses.csv"
EVENTS_CSV = REPO_DIR / "output/catalogs/05_compound_heatwave_events/compound_events.csv"
EVENT_DAYS_CSV = EVENTS_CSV.with_name("event_days.csv")
INPUT_DIR = REPO_DIR / "supporting_files/additional_results/modeling/power_law"


def reread_with_pandas_default(values):
    """Write values with pandas.to_csv and read them back with the default
    pandas reader used for the manuscript area samples."""
    buffer = io.StringIO()
    pd.DataFrame({"area_km2": values}).to_csv(buffer, index=False)
    buffer.seek(0)
    return pd.read_csv(buffer)["area_km2"].to_numpy()


def exact_float_column(path, column):
    """Read one column with an exact (correctly rounded) number parser."""
    table = pd.read_csv(path, usecols=[column], dtype={column: str})
    return np.array([float(text) for text in table[column]])


def write_exact_csv(table, path):
    """Floats with 17 significant digits; everything else as text."""
    columns = [table[c].map(lambda v: "%.17g" % v) if table[c].dtype == float
               else table[c].astype(str) for c in table.columns]
    lines = [",".join(table.columns)]
    lines += [",".join(row) for row in zip(*columns)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def units_in_last_place(a, b):
    """How many representable doubles lie between a and b (0 = identical)."""
    ia = np.asarray(a, dtype="<f8").view("<i8").astype(np.int64)
    ib = np.asarray(b, dtype="<f8").view("<i8").astype(np.int64)
    return np.abs(ia - ib)


def main():
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Select maxima using values read by the default pandas number reader.
    ellipses = f.attach_compound_events(pd.read_csv(ELLIPSES_CSV),
                                       pd.read_csv(EVENT_DAYS_CSV))
    exact_area = exact_float_column(ELLIPSES_CSV, "area_km2")

    in_catalog_order = ellipses.sort_values(
        ["date", "cluster_number_in_day"], kind="stable").index
    assert (in_catalog_order == ellipses.index).all(), "catalog is not in date order"

    # daily maxima
    # Each date belongs to one event, so grouping by (event, date) gives
    # one maximum per selected day, in date order.
    daily_rows = ellipses.groupby(["compound_event_id", "date"])["area_km2"].idxmax()
    daily = ellipses.loc[daily_rows.to_numpy()].sort_values("date")
    assert (daily.index == daily_rows.to_numpy()).all()
    daily_ties = int((ellipses.groupby("date")["area_km2"]
                      .apply(lambda a: (a == a.max()).sum()) > 1).sum())

    # event maxima
    event_rows = ellipses.groupby("compound_event_id")["area_km2"].idxmax()
    event = ellipses.loc[event_rows.to_numpy()]
    event_ties = int((ellipses.groupby("compound_event_id")["area_km2"]
                      .apply(lambda a: (a == a.max()).sum()) > 1).sum())

    # The compound-event catalog records the same largest ellipse per event.
    events = pd.read_csv(EVENTS_CSV).set_index("compound_event_id")
    assert (events.loc[event["compound_event_id"], "largest_ellipse_id"].to_numpy()
            == event["ellipse_id"].to_numpy()).all()

    print("daily maxima: %d selected days, %d with a tie for the largest area"
          % (len(daily), daily_ties))
    print("event maxima: %d events, %d with a tie for the largest area"
          % (len(event), event_ties))

    # parse the selected areas at the precision used for the fits
    for name, table, columns in (
            ("daily_max_area.csv", daily,
             ["date", "compound_event_id", "ellipse_id", "cluster_id", "area_km2"]),
            ("event_max_area.csv", event,
             ["compound_event_id", "ellipse_id", "cluster_id", "date", "area_km2"])):
        out = table[columns].reset_index(drop=True).copy()
        catalog_value = exact_area[table.index.to_numpy()]
        out["area_km2"] = reread_with_pandas_default(out["area_km2"].to_numpy())
        ulp = units_in_last_place(out["area_km2"], catalog_value)
        print("%s: %d of %d areas differ from the exact catalog value "
              "(largest difference %d units in the last place, %.2g km2)"
              % (name, (ulp > 0).sum(), len(out), ulp.max(),
                 np.abs(out["area_km2"].to_numpy() - catalog_value).max()))
        path = INPUT_DIR / name
        write_exact_csv(out, path)
        # The file must give back exactly these numbers.
        assert np.array_equal(exact_float_column(path, "area_km2"),
                              out["area_km2"].to_numpy())


if __name__ == "__main__":
    main()
