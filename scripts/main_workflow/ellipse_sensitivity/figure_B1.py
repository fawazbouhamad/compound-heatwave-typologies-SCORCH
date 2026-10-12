"""Figure B1: ellipse purity across axis scales.

Writes the full PNG and its panels to output/figures/figure_B1/.
Run: python scripts/main_workflow/ellipse_sensitivity/figure_B1.py
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

def draw_matrix(ax, fig, matrix, colorbar_label):
    """One heatmap panel: values printed in each cell, diagonal outlined."""
    rows_scale = [float(s) for s in matrix.index]
    cols_scale = [float(c) for c in matrix.columns]
    values = matrix.to_numpy(float)
    image = ax.imshow(values, cmap="viridis", aspect="equal")
    vmin, vmax = np.nanmin(values), np.nanmax(values)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            # White text on the dark part of the colour scale.
            frac = (values[i, j] - vmin) / (vmax - vmin) if vmax > vmin else 1
            ax.text(j, i, "%.0f" % values[i, j], ha="center", va="center",
                    fontsize=7.2, color="white" if frac < 0.55 else "black")
        if i < values.shape[1]:   # equal scale factors on both axes
            ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor="magenta", linewidth=1.6,
                                       zorder=5))
    ax.set_xticks(range(len(cols_scale)), [f"{s:g}" for s in cols_scale], fontsize=9)
    ax.set_yticks(range(len(rows_scale)), [f"{s:g}" for s in rows_scale], fontsize=9)
    # Axis wording as printed in the manuscript figure (sigma = scale factor).
    ax.set_xlabel(r"$\sigma_{PC2}$ (Minor Axis)", fontsize=11)
    ax.set_ylabel(r"$\sigma_{PC1}$ (Major Axis)", fontsize=11)
    box = ax.get_position()
    colorbar_ax = fig.add_axes([box.x1 + 0.012, box.y0, 0.012, box.height])
    colorbar = fig.colorbar(image, cax=colorbar_ax)
    colorbar.set_label(colorbar_label, fontsize=10)
    colorbar.ax.tick_params(labelsize=8.5)

def make_figure(output_dir, out_stem):
    """Figure B.1 drawn from the saved matrix CSV files (layout of the manuscript figure)."""
    plt.rcParams.update({"font.family": "Arial"})
    panels = [("a", "minimum", "Minimum (%)"), ("b", "maximum", "Maximum (%)"),
              ("c", "mean", "Average (%)"), ("d", "sum", "Combined score")]
    positions = {"a": (0, 0), "b": (0, 1), "c": (1, 0), "d": (1, 1)}
    x0 = {0: 0.075, 1: 0.575}
    y0 = {0: 0.565, 1: 0.075}
    fig = plt.figure(figsize=(11.6, 10.6))
    axes = {}
    for letter, name, label in panels:
        matrix = pd.read_csv(output_dir / f"purity_matrix_{name}.csv", index_col=0)
        row, col = positions[letter]
        ax = fig.add_axes([x0[col], y0[row], 0.335, 0.365])
        draw_matrix(ax, fig, matrix, label)
        axes[letter] = ax
    fig.canvas.draw()
    for letter, ax in axes.items():
        figure_style.label_above(fig, ax, f"({letter})")
    paths = save_figure(fig, out_stem, dpi=600)
    plt.close(fig)

    return paths


def main():
    config = f.load_config()
    results = f.supporting_results_dir(config) / "main_workflow/ellipse_sensitivity"
    with mpl.rc_context():
        mpl.rcdefaults()
        paths = make_figure(results, f.output_subdir(config, "figures") / "figure_B1")
    for path in paths:
        print(f.relative_to_repo(path))
    return paths


if __name__ == "__main__":
    main()
