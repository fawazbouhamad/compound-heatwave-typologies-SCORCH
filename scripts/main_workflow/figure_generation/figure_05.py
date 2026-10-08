"""Figure 5: representative Type 1 and Type 2 events.

Writes the full PNG and its panels to output/figures/figure_5/.
Run: python scripts/main_workflow/figure_generation/figure_05.py
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
from snapshot_layout import check_event, _render, _grid, _save

LEGEND_AXIS_LABELS = (r"$L_1$", r"$L_2$")

def make_figure_05(cluster_cells, ellipses, compound_events, dates,
                   out_stem):
    """dates: {"type_1": [date], "type_2": [date]}."""
    check_event(compound_events, dates["type_1"], 12, 1)
    check_event(compound_events, dates["type_2"], 1, 2)
    fig, axes = _grid(1, 2, (12.4, 5.4), top=0.86, bottom=0.14)
    _render(axes[0, 0], cluster_cells, ellipses, dates["type_1"][0])
    _render(axes[0, 1], cluster_cells, ellipses, dates["type_2"][0])
    for j, type_number in ((0, 1), (1, 2)):
        b = axes[0, j].get_position()
        Q.type_heading(fig, (b.x0 + b.x1) / 2.0, 0.93, type_number,
                       fontsize=24)
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
    out_stem = f.output_subdir(config, "figures") / "figure_5"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_05(cluster_cells, ellipses, events, {"type_1": dates["figure_5_type_1"], "type_2": dates["figure_5_type_2"]}, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
