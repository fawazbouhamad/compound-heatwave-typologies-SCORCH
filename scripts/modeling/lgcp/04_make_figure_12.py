"""Draw Figure 12 from the fitted LGCP.

Reads: supporting_files/additional_results/modeling/lgcp/
Writes: output/figures/figure_12/
Run: python scripts/modeling/lgcp/04_make_figure_12.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("SCORCH_REQUIRE_SLIDE_FONTS", "0")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.lines import Line2D

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "figure_generation"))
from figure_exports import save_figure  # noqa: E402
import figure_03 as ref  # noqa: E402  (map style)
import cartopy.crs as ccrs  # noqa: E402

OUTPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp"
FIGURES_DIR = REPO_DIR / "output" / "figures"

INTENSITY_COLUMN = "mean_intensity_per_million_km2"
CBAR_LABEL = "Relative centroid-concentration rank, R(s)"
# Colour classes: the full 0 to 1 rank range in steps of 0.1.
RANK_BOUNDS = np.round(np.arange(0.0, 1.0001, 0.1), 1)


def main():
    started = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    pred = pd.read_csv(OUTPUT_DIR / "predicted_intensity_1deg.csv")
    assert pred[INTENSITY_COLUMN].notna().all(), "missing predicted intensities"
    n_cells = len(pred)
    # Rank R(s) of the fitted intensity among all cells (0-1).
    pred["concentration_rank_R"] = pred[INTENSITY_COLUMN].rank(method="average", pct=True)
    # The same ranks as written by 02_fit_lgcp.R (that file keeps 15 digits,
    # so the average ranks, rank x 1800, are compared).
    if not np.array_equal(np.rint(pred["concentration_rank_R"] * 2 * n_cells),
                          np.rint(pred["concentration_rank"] * 2 * n_cells)):
        raise ValueError("R(s) differs from the rank written by 02_fit_lgcp.R")

    lons = np.sort(pred["lon"].unique())
    lats = np.sort(pred["lat"].unique())
    grid2d = (pred.pivot(index="lat", columns="lon", values="concentration_rank_R")
              .reindex(index=lats, columns=lons).to_numpy())

    # float_precision="round_trip" reads the stored coordinates exactly.
    m = pd.read_csv(OUTPUT_DIR / "centroids.csv", float_precision="round_trip")
    is_daily = m["is_daily_largest"].to_numpy(bool)
    is_event = m["is_event_largest"].to_numpy(bool)
    print(f"cells = {n_cells}; centroids = {len(m)}; daily-largest = {is_daily.sum()}; "
          f"event-largest = {is_event.sum()}")

    # single map, drawing code of the manuscript figure
    fig, ax = plt.subplots(1, 1, figsize=(7.2, 5.6),
                           subplot_kw={"projection": ccrs.PlateCarree()})
    ref.setup_panel_ax(ax)

    norm = BoundaryNorm(RANK_BOUNDS, ncolors=ref.HEAT_CMAP.N, clip=True)
    im = ax.pcolormesh(lons, lats, grid2d, transform=ccrs.PlateCarree(),
                       cmap=ref.HEAT_CMAP, norm=norm, shading="auto",
                       alpha=0.78, zorder=3)
    tf = dict(transform=ccrs.PlateCarree())

    # 1) every centroid: small black dot
    ax.scatter(m["centroid_lon"], m["centroid_lat"],
               s=ref.POINT_SIZE, c="black", alpha=ref.POINT_ALPHA,
               edgecolors="none", zorder=8, **tf)
    # 2) largest ellipse of each day: open ring around the dot
    ax.scatter(m.loc[is_daily, "centroid_lon"], m.loc[is_daily, "centroid_lat"],
               s=70, facecolors="none", edgecolors="black", linewidths=0.9,
               alpha=0.85, zorder=9, **tf)
    # 3) largest ellipse of each event: red cross on a white halo, on top
    ax.scatter(m.loc[is_event, "centroid_lon"], m.loc[is_event, "centroid_lat"],
               s=62, marker="x", c="white", linewidths=3.2,
               alpha=1.0, zorder=10, **tf)
    ax.scatter(m.loc[is_event, "centroid_lon"], m.loc[is_event, "centroid_lat"],
               s=55, marker="x", c="#b2182b", linewidths=1.6,
               alpha=1.0, zorder=11, **tf)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markeredgecolor="black", markersize=5, label="Centroid"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="black", markeredgewidth=1.0, markersize=9,
               label="Daily-largest structure centroid"),
        Line2D([0], [0], marker="x", color="#b2182b", markeredgewidth=1.6,
               markersize=7, linestyle="none",
               label="Event-largest structure centroid"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.045),
              ncol=3, fontsize=8.2, frameon=True, fancybox=False, framealpha=1.0,
              edgecolor="black", facecolor="white", handlelength=2.0,
              columnspacing=1.4, borderpad=0.65)

    # external vertical colour bar with decimal ticks
    plt.tight_layout(rect=[0, 0.04, 0.88, 0.98])
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.025, pos.y0, 0.020, pos.height])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical", ticks=RANK_BOUNDS)
    cbar.set_label(CBAR_LABEL, fontsize=10)
    cbar.ax.set_yticklabels([f"{b:.1f}" for b in RANK_BOUNDS])
    cbar.ax.tick_params(labelsize=9, length=2.5, width=0.6)
    cbar.outline.set_linewidth(0.6)

    save_figure(fig, FIGURES_DIR / "figure_12", dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
