"""Figure 6: representative Type 3 events.

Writes the full PNG and its panels to output/figures/figure_6/.
Run: python scripts/main_workflow/figure_generation/figure_06.py
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

ROW_LABEL_GAP_PT = 14.0

def _row_labels_measured(fig, axes):
    """(a)/(b) to the left of the two rows, 14 pt clear of the leftmost
    tick labels (at least 12 pt is required), at one shared x position and
    vertically centred on each row."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72.0
    fig_w_px = fig.get_size_inches()[0] * fig.dpi

    row_left_px = [axes[i, 0].get_tightbbox(renderer).x0 for i in range(2)]
    label_right_px = min(row_left_px) - ROW_LABEL_GAP_PT * px_per_pt
    fx = label_right_px / fig_w_px

    texts = []
    for i, letter in ((0, "(a)"), (1, "(b)")):
        boxes = [axes[i, j].get_position() for j in range(2)]
        yc = (min(b.y0 for b in boxes) + max(b.y1 for b in boxes)) / 2.0
        texts.append((i, fig.text(fx, yc, letter, ha="right", va="center",
                                  **Q.PANEL_LABEL)))

    fig.canvas.draw()
    for i, text in texts:
        box = text.get_window_extent(renderer)
        assert (row_left_px[i] - box.x1) / px_per_pt >= 12.0
        assert box.x0 >= 0.0

def make_figure_06(cluster_cells, ellipses, compound_events, dates,
                   out_stem):
    """dates: {"event_a": [day t, day t+1], "event_b": [day t, day t+1]}."""
    check_event(compound_events, dates["event_a"], 29, 3)
    check_event(compound_events, dates["event_b"], 3, 3)
    fig, axes = _grid(2, 2, (11.6, 9.4), top=0.865, bottom=0.095)
    days = [dates["event_a"], dates["event_b"]]
    for i in range(2):
        for j in range(2):
            _render(axes[i, j], cluster_cells, ellipses, days[i][j])
    Q.type_heading(fig, 0.5, 0.945, 3, fontsize=26)
    _time_label(fig, axes[0, 0], "t")           # once, above the top row
    _time_label(fig, axes[0, 1], "t+1")
    _row_labels_measured(fig, axes)
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
    out_stem = f.output_subdir(config, "figures") / "figure_6"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_06(cluster_cells, ellipses, events, {"event_a": dates["figure_6a_type_3"], "event_b": dates["figure_6b_type_3"]}, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
