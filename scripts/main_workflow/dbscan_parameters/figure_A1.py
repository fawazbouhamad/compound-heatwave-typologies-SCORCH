"""Figure A1: DBSCAN parameter-selection matrices.

Writes the full PNG and its panels to output/figures/figure_A1/.
Run: python scripts/main_workflow/dbscan_parameters/figure_A1.py
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
import figure_style
from figure_exports import save_figure

TEXT_WIDTH_IN = 9026.0 / 1440.0
RETAINED_FACE = "#bcd9ec"
RETAINED_EDGE = "#2a6ea6"

def draw_selection_matrix(ax, matrix, retained, show_ylabel, show_xlabel):
    """One day's matrix: cluster count per (min_samples, eps); retained pairs shaded."""
    min_samples_values = [int(m) for m in matrix.index]
    eps_values = [float(e) for e in matrix.columns]
    counts = matrix.to_numpy(float)
    n_rows, n_cols = counts.shape
    retained_pairs = {(int(r.min_samples), float(r.eps)) for r in retained.itertuples()}
    for i in range(n_rows):
        for j in range(n_cols):
            is_retained = (min_samples_values[i], eps_values[j]) in retained_pairs
            if is_retained:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           facecolor=RETAINED_FACE, edgecolor=RETAINED_EDGE,
                                           linewidth=1.4, zorder=2))
            ax.text(j, i, f"{counts[i, j]:.0f}", ha="center", va="center", fontsize=8.5,
                    zorder=3, fontweight="bold" if is_retained else "normal", color="black")
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.set_xticks(range(n_cols), [f"{e:g}" for e in eps_values], fontsize=9)
    ax.set_yticks(range(n_rows), [str(m) for m in min_samples_values], fontsize=9)
    if show_xlabel:
        ax.set_xlabel("eps (grid-cell units)", fontsize=11)
    if show_ylabel:
        ax.set_ylabel("min_samples", fontsize=11)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="0.85", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

def make_figure_a1(days, output_dir, out_stem):
    """Six panels, (a)-(f), two per row; drawn from the saved matrix CSV files."""
    plt.rcParams.update({"font.family": "Arial"})
    # Layout in inches (panel size, margins, gaps), rounded to whole pixels at 600 dpi.
    panel_w, panel_h = 3.55, 3.24
    left, col_gap, right = 0.85, 0.75, 0.20
    label_space, row_gap, top = 0.42, 0.30, 0.12
    xlabel_space, bottom = 0.68, 0.10
    fig_w = left + 2 * panel_w + col_gap + right
    fig_h = top + 3 * (label_space + panel_h) + 2 * row_gap + xlabel_space + bottom
    fig_w = round(fig_w * 600) / 600.0
    fig_h = round(fig_h * 600) / 600.0

    retained_all = pd.read_csv(output_dir / "figure_A1_retained_pairs.csv")
    fig = plt.figure(figsize=(fig_w, fig_h))
    axes = {}
    for k, day in enumerate(days):
        row, col = divmod(k, 2)
        x = left + col * (panel_w + col_gap)
        y = bottom + xlabel_space + (2 - row) * (panel_h + label_space + row_gap)
        matrix = pd.read_csv(output_dir / f"figure_A1_matrix_{day}.csv", index_col=0)
        retained = retained_all[retained_all["date"] == day]
        ax = fig.add_axes([x / fig_w, y / fig_h, panel_w / fig_w, panel_h / fig_h])
        draw_selection_matrix(ax, matrix, retained, show_ylabel=(col == 0),
                              show_xlabel=(row == 2))
        axes[chr(ord("a") + k)] = ax
    fig.canvas.draw()
    for letter, ax in axes.items():
        figure_style.label_above(fig, ax, f"({letter})")
    paths = save_figure(fig, out_stem, dpi=600,
                        physical_width=TEXT_WIDTH_IN)
    plt.close(fig)
    return paths


def main():
    config = f.load_config()
    results = f.supporting_results_dir(config) / "main_workflow/dbscan_parameters"
    days = config["representative_events"]["figure_7_type_4"]
    with mpl.rc_context():
        mpl.rcdefaults()
        paths = make_figure_a1(days, results, f.output_subdir(config, "figures") / "figure_A1")
    for path in paths:
        print(f.relative_to_repo(path))
    return paths


if __name__ == "__main__":
    main()
