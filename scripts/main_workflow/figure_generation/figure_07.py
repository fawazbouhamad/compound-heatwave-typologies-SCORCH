"""Figure 7: representative Type 4 event.

Writes the full PNG and its panels to output/figures/figure_7/.
Run: python scripts/main_workflow/figure_generation/figure_07.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("SCORCH_REQUIRE_SLIDE_FONTS", "0")
WORKFLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW))
sys.path.insert(0, str(WORKFLOW / "helpers"))

import matplotlib as mpl
import pandas as pd
import functions as f

import figure_style as Q
from snapshot_layout import check_event, _render, _grid, _time_label, _save

LEGEND_AXIS_LABELS = (r"$L_1$", r"$L_2$")

def make_figure_07(cluster_cells, ellipses, compound_events, dates,
                   out_stem):
    """dates: the six days of Event 10 in order."""
    check_event(compound_events, dates, 10, 4)
    fig, axes = _grid(2, 3, (15.6, 8.9), top=0.855, bottom=0.095,
                      hspace=0.34)   # room for the t+3..t+5 labels
    for k, (ax, date) in enumerate(zip(axes.reshape(-1), dates)):
        _render(ax, cluster_cells, ellipses, date)
        _time_label(fig, ax, "t" if k == 0 else f"t+{k}")
    Q.type_heading(fig, 0.5, 0.945, 4, fontsize=26)
    Q.unified_snapshot_legend(fig, axis_labels=LEGEND_AXIS_LABELS)
    return _save(fig, out_stem)


def main():
    config = f.load_config()
    inputs = {
        "cluster_cells": f.catalog_dir(config, "03_extreme_heatwave_clusters") / "cluster_cells.csv",
        "ellipses": f.catalog_dir(config, "04_extreme_heatwave_ellipses") / "ellipses.csv",
        "events": f.catalog_dir(config, "05_compound_heatwave_events") / "compound_events.csv",
    }
    cluster_cells = pd.read_csv(inputs["cluster_cells"], usecols=["date", "lon", "lat", "dbscan_label"])
    ellipses = pd.read_csv(inputs["ellipses"])
    events = pd.read_csv(inputs["events"])
    dates = config["representative_events"]
    out_stem = f.output_subdir(config, "figures") / "figure_7"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_07(cluster_cells, ellipses, events, dates["figure_7_type_4"], out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
