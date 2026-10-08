"""Figure 2: local thresholds and daily heatwave extent.

Panel (a) maps each cell's warm-season 95th-percentile Tmax for 1940–2025.
Panel (b) compares the fraction of domain cells in heatwaves with the fraction
of threshold-exceeding cells in heatwaves, using hexagonal bins on a log scale.
The red line marks 371 heatwave cells out of 1,800 (the 97.5th percentile).
The two panels are drawn at 600 dpi and combined with panel letters above.

Writes the full PNG and its panels to output/figures/figure_2/.
Run: python scripts/main_workflow/figure_generation/figure_02.py
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

import io

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import figure_style as Q
from figure_exports import save_figure  # noqa: E402

LON_MIN, LON_MAX = 20.0, 70.0
LAT_MIN, LAT_MAX = 10.0, 46.0
N_DOMAIN_CELLS = 1800
HEXBIN_GRIDSIZE = 50
PANEL_DPI = 600


def _png_image(fig):
    """Save a panel as a 600 dpi PNG (cropped to its content) and read the
    pixels back, exactly as the manuscript figure was assembled."""
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=PANEL_DPI, bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return plt.imread(buffer)


def draw_panel_a(thresholds):
    """Map of the local 95th-percentile Tmax thresholds (degC)."""
    lats = np.sort(thresholds["lat"].unique())
    lons = np.sort(thresholds["lon"].unique())
    grid = (thresholds.pivot(index="lat", columns="lon",
                             values="threshold_tmax_c")
            .reindex(index=lats, columns=lons).to_numpy())

    # Cell edges for pcolormesh (cell centres are at .5 degrees).
    lon_e = np.arange(lons.min() - 0.5, lons.max() + 0.51, 1.0)
    lat_e = np.arange(lats.min() - 0.5, lats.max() + 0.51, 1.0)

    fig, ax = plt.subplots(figsize=(9.0, 6.0),
                           subplot_kw={"projection": ccrs.PlateCarree()})
    im = ax.pcolormesh(lon_e, lat_e, grid, transform=ccrs.PlateCarree(),
                       cmap="YlOrRd", edgecolors="black", linewidth=0.25)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.7)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.5)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.0,
                      xlocs=np.arange(30, 70, 10),
                      ylocs=np.arange(20, 41, 10))
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {"size": 11}
    gl.ylabel_style = {"size": 11}

    cbar = fig.colorbar(im, ax=ax, orientation="vertical",
                        fraction=0.035, pad=0.02)
    cbar.set_label("Tmax p95 (°C)", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    fig.tight_layout()
    return _png_image(fig)


def draw_panel_b(daily, regional_threshold_cells):
    """Hexbin of the daily heatwave fraction against the heatwave share."""
    threshold_fraction = regional_threshold_cells / N_DOMAIN_CELLS

    # Plotted days: at least one heatwave cell and one exceeding cell.
    sub = daily[(daily["n_heatwave_cells"] > 0)
                & (daily["n_exceeding_cells"] > 0)].copy()
    x = (sub["n_heatwave_cells"] / N_DOMAIN_CELLS).to_numpy(float)
    y = (sub["n_heatwave_and_exceeding_cells"]
         / sub["n_exceeding_cells"]).to_numpy(float)
    assert 0.0 <= float(y.min()) and float(y.max()) <= 1.0

    diagnostics = {
        "n_plotted_days": len(sub),
        "min_share": float(y.min()),
        "max_share": float(y.max()),
        "median_share": float(np.median(y)),
        "p95_share": float(np.percentile(y, 95)),
        "p97_5_share": float(np.percentile(y, 97.5)),
        "p99_share": float(np.percentile(y, 99)),
        "n_below_0_3": int((y < 0.3).sum()),
        "n_equal_0": int((y == 0.0).sum()),
        "n_equal_1": int((y == 1.0).sum()),
        "n_above_1": int((y > 1.0).sum()),
        "n_days_without_exceeding_cells": int(
            (daily["n_exceeding_cells"] == 0).sum()),
        "regional_threshold_cells": int(regional_threshold_cells),
        "regional_threshold_fraction": threshold_fraction,
    }

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    hb = ax.hexbin(x, y, gridsize=HEXBIN_GRIDSIZE, cmap="viridis",
                   norm=LogNorm(), mincnt=1, linewidths=0.1,
                   extent=(0.0, 0.6, 0.0, 1.0))
    diagnostics["max_days_in_one_hexagon"] = int(hb.get_array().max())

    ax.axvline(threshold_fraction, color="red", linestyle="--",
               linewidth=1.8, zorder=5)
    ax.text(threshold_fraction - 0.006, 0.03, "P97.5", color="red",
            rotation=90, ha="right", va="bottom", fontsize=13,
            fontweight="bold", transform=ax.get_xaxis_transform())

    ax.set_xlim(0.0, 0.6)
    ax.set_ylim(0.0, 1.01)
    ax.set_xlabel("Fraction of Domain Grid Cells with Heatwave Label",
                  fontsize=14, fontweight="bold")
    ax.set_ylabel("Heatwave-Labeled Fraction of\nTmax-Exceedance Grid Cells",
                  fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=12)

    cax = make_axes_locatable(ax).append_axes("right", size="4.5%", pad=0.1)
    cbar = fig.colorbar(hb, cax=cax)
    cbar.set_label("Number of days", fontsize=13)
    cbar.ax.tick_params(labelsize=11)
    fig.tight_layout()
    return _png_image(fig), diagnostics


def make_figure_02(thresholds, daily, out_stem):
    """Save Figure 2 and its panel PNGs at 600 dpi; return their paths."""
    regional = daily["regional_threshold_cells"].unique()
    assert len(regional) == 1
    regional_threshold_cells = float(regional[0])

    panel_a = draw_panel_a(thresholds)
    panel_b, _ = draw_panel_b(daily, regional_threshold_cells)

    images = [panel_a, panel_b]
    h = max(im.shape[0] / im.shape[1] for im in images) / 2.0
    fig_w = 14.0
    fig = plt.figure(figsize=(fig_w, fig_w * h * 1.06))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.02,
                           left=0.005, right=0.995, top=0.93, bottom=0.005)
    for i, (letter, im) in enumerate(zip(("a", "b"), images)):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(im)
        ax.set_axis_off()
        Q.label_above(fig, ax, f"({letter})")

    written = save_figure(fig, out_stem, dpi=600)
    plt.close(fig)

    return written

def main():
    config = f.load_config()
    inputs = {
        "thresholds": f.supporting_results_dir(config) / "main_workflow/heatwave_detection/cell_thresholds.csv",
        "daily": f.supporting_results_dir(config) / "main_workflow/regional_selection/daily_summary.csv",
    }
    thresholds = pd.read_csv(inputs["thresholds"])
    daily = pd.read_csv(inputs["daily"])
    out_stem = f.output_subdir(config, "figures") / "figure_2"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_02(thresholds, daily, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
