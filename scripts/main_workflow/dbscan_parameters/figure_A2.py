"""Figure A2: selected DBSCAN settings by event type.

Writes the full PNG and its panels to output/figures/figure_A2/.
Run: python scripts/main_workflow/dbscan_parameters/figure_A2.py
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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figure_style
from figure_exports import save_figure
from functions import DBSCAN_GROUPS as GROUPS, figure_a2_panels

GROUP_COLORS = {"All Types": "#666666", "Type 1": "#d7191c", "Type 2": "#f57c00",
                "Type 3": "#d9a300", "Type 4": "#6a3d9a"}

def draw_boxes(ax, groups, ylabel):
    stats, colors = [], []
    for group in GROUPS:
        v = np.sort(groups[group])
        stats.append(dict(label=group, med=np.median(v), mean=v.mean(),
                          q1=np.percentile(v, 25), q3=np.percentile(v, 75),
                          whislo=v.min(), whishi=v.max(), fliers=[]))
        colors.append(GROUP_COLORS[group])
    artists = ax.bxp(stats, showmeans=True, showfliers=False, widths=0.52,
                     meanprops=dict(marker="D", markerfacecolor="white",
                                    markeredgecolor="black", markersize=6),
                     patch_artist=True)
    for box, median, color in zip(artists["boxes"], artists["medians"], colors):
        box.set_facecolor(color)
        box.set_alpha(0.55)
        box.set_edgecolor(color)
        median.set_color(color)
        median.set_linewidth(2.6)
    for i, (w1, w2, c1, c2) in enumerate(zip(artists["whiskers"][::2], artists["whiskers"][1::2],
                                             artists["caps"][::2], artists["caps"][1::2])):
        for artist in (w1, w2, c1, c2):
            artist.set_color(colors[i])
            artist.set_linewidth(1.2)
    ax.set_ylabel(ylabel, fontsize=11.5)
    ax.grid(True, axis="y", alpha=0.22, linewidth=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9.5, length=3)

def make_figure_a2(panels, out_stem):
    plt.rcParams.update({"font.family": "Arial"})
    fig = plt.figure(figsize=(12.6, 8.6))
    x0 = {0: 0.075, 1: 0.575}
    y0 = {0: 0.585, 1: 0.145}
    axes = {}
    for letter, (row, col) in (("a", (0, 0)), ("b", (0, 1)), ("c", (1, 0)), ("d", (1, 1))):
        ax = fig.add_axes([x0[col], y0[row], 0.40, 0.345])
        groups, ylabel, _ = panels[letter]
        draw_boxes(ax, groups, ylabel)
        axes[letter] = ax
    fig.canvas.draw()
    for letter, ax in axes.items():
        figure_style.label_above(fig, ax, f"({letter})")
    handles = [Line2D([0], [0], color="0.2", lw=2.6, label="median"),
               Line2D([0], [0], marker="D", color="w", markerfacecolor="white",
                      markeredgecolor="black", markersize=6, label="mean"),
               Line2D([0], [0], color="0.2", lw=1.2, label="whiskers = min / max")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.015), ncol=3,
               frameon=False, fontsize=10.5, handlelength=2.2, columnspacing=1.8)
    paths = save_figure(fig, out_stem, dpi=600)
    plt.close(fig)
    return paths


def main():
    config = f.load_config()
    inputs = {
        "days": f.supporting_results_dir(config) / "main_workflow/dbscan_parameters/figure_A2_values_per_day.csv",
        "events": f.supporting_results_dir(config) / "main_workflow/dbscan_parameters/figure_A2_values_per_event.csv",
    }
    panels = figure_a2_panels(pd.read_csv(inputs["days"]), pd.read_csv(inputs["events"]))
    out_stem = f.output_subdir(config, "figures") / "figure_A2"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_a2(panels, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
