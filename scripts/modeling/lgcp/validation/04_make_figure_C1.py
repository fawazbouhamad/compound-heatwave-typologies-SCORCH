"""Draw Figure C1 from the five-fold validation results.

Reads: Prepared LGCP inputs and Appendix C results.
Writes: output/figures/figure_C1/
Run: python scripts/modeling/lgcp/validation/04_make_figure_C1.py
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

REPO_DIR = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "figure_generation"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "modeling" / "lgcp"))
import functions as f  # noqa: E402
import figure_style as Q  # noqa: E402
from figure_exports import save_figure  # noqa: E402
import figure_03 as ref  # noqa: E402  (map style)
import zone_scoring as zs  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402

OUTPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp" / "validation"
FIGURES_DIR = REPO_DIR / "output" / "figures"
LGCP_INPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp"

# Layout of manuscript Figure C1.
RANK_BOUNDS = np.round(np.arange(0.0, 1.0001, 0.1), 1)
CBAR_LABEL = "Relative centroid-concentration rank, R(s)"
AX_W = 4.30
AX_H = AX_W * 36.0 / 50.0
L_M, COL_GAP, ROW_GAP, TOP_M, BOT_M = 1.05, 1.25, 0.75, 0.45, 0.25
CB_GUTTER = 1.05
FIG_W = L_M + 2 * AX_W + 1 * COL_GAP + CB_GUTTER + 0.15
FIG_H = TOP_M + 3 * AX_H + 2 * ROW_GAP + BOT_M + 0.55
FIG_W = round(FIG_W * 500) / 500.0
FIG_H = round(FIG_H * 500) / 500.0
TEXT_W_IN = 9026.0 / 1440.0          # manuscript text width, 6.2681 in = 159.21 mm
ZONE_COLUMNS = [("Top 10%", "distance_to_top10_zone_km", "mean_distance_to_top10_zone_km"),
                ("Top 20%", "distance_to_top20_zone_km", "mean_distance_to_top20_zone_km"),
                ("Top 30%", "distance_to_top30_zone_km", "mean_distance_to_top30_zone_km"),
                ("Top 50%", "distance_to_top50_zone_km", "mean_distance_to_top50_zone_km")]


def area_weighted_rank_grid(cells, intensity, area):
    """Area-weighted rank R(s) (fraction of domain area with intensity <= the
    cell's), as a latitude x longitude array for plotting."""
    total = float(area.sum())
    ascending = np.argsort(intensity)
    rank = np.empty_like(intensity)
    rank[ascending] = np.cumsum(area[ascending]) / total
    table = cells.copy()
    table["risk_q"] = rank
    lons = np.sort(table["lon"].unique())
    lats = np.sort(table["lat"].unique())
    grid2d = (table.pivot(index="lat", columns="lon", values="risk_q")
              .reindex(index=lats, columns=lons).to_numpy())
    return lons, lats, grid2d, rank


def main():
    started = time.time()
    config = f.load_config()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial"})

    grid = pd.read_csv(LGCP_INPUT_DIR / "temperature_covariates.csv").sort_values(
        "cell_id").reset_index(drop=True)
    grid["area_km2"] = zs.cell_area_km2(grid["lat"])
    cen = pd.read_csv(LGCP_INPUT_DIR / "centroids.csv")
    n_cen = len(cen)
    assert n_cen == 753, n_cen

    pc = pd.read_csv(OUTPUT_DIR / "validation_per_centroid.csv").sort_values(
        "ellipse_id").reset_index(drop=True)
    assert len(pc) == n_cen and pc["fold"].nunique() == config["appendix_c"]["n_folds"]
    assert np.array_equal(pc["ellipse_id"], cen["ellipse_id"])
    # Check that the folds match the seeded split in the shared study settings.
    rng = np.random.default_rng(config["appendix_c"]["fold_seed"])
    folds = np.array_split(rng.permutation(n_cen), config["appendix_c"]["n_folds"])
    fold_of = np.empty(n_cen, int)
    for k, idx in enumerate(folds, start=1):
        fold_of[idx] = k
    assert (pc["fold"].to_numpy(int) == fold_of).all()

    cen["is_daily_max"] = pc["is_daily_largest"].astype(bool).to_numpy()
    cen["is_event_max"] = pc["is_event_largest"].astype(bool).to_numpy()
    fm = pd.read_csv(OUTPUT_DIR / "fold_mean_intensity_1deg.csv")
    summ = pd.read_csv(OUTPUT_DIR / "validation_group_summary.csv")
    s_all = summ[summ["group"] == "all_centroids"].iloc[0]

    area = grid["area_km2"].to_numpy(float)
    norm = BoundaryNorm(RANK_BOUNDS, ncolors=ref.HEAT_CMAP.N, clip=True)

    # slot 0 = distance box plot, slots 1-5 = folds 1-5
    positions = {}
    for k in range(6):
        r, c = divmod(k, 2)
        x = L_M + c * (AX_W + COL_GAP)
        y = BOT_M + 0.55 + (2 - r) * (AX_H + ROW_GAP)
        positions[k] = (x, y)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    axes, im = {}, None
    for k, test_idx in enumerate(folds):
        fold_no = k + 1
        test_mask = np.zeros(n_cen, bool)
        test_mask[test_idx] = True
        test = cen.loc[test_mask]

        sub = fm[fm["fold"] == fold_no].set_index("cell_id").reindex(grid["cell_id"])
        assert sub["mean_intensity_per_km2"].notna().all()
        lam = sub["mean_intensity_per_km2"].to_numpy(float)
        lons, lats, grid2d, _ = area_weighted_rank_grid(grid, lam, area)

        x, y = positions[fold_no]
        ax = fig.add_axes([x / FIG_W, y / FIG_H, AX_W / FIG_W, AX_H / FIG_H],
                          projection=ccrs.PlateCarree())
        ref.setup_panel_ax(ax)
        im = ax.pcolormesh(lons, lats, grid2d, transform=ccrs.PlateCarree(),
                           cmap=ref.HEAT_CMAP, norm=norm, shading="auto",
                           alpha=0.78, zorder=3)
        tf = dict(transform=ccrs.PlateCarree())
        is_daily = test["is_daily_max"].to_numpy(bool)
        is_event = test["is_event_max"].to_numpy(bool)
        ax.scatter(test["centroid_lon"], test["centroid_lat"], s=ref.POINT_SIZE,
                   c="black", alpha=ref.POINT_ALPHA, edgecolors="none", zorder=8, **tf)
        ax.scatter(test.loc[is_daily, "centroid_lon"], test.loc[is_daily, "centroid_lat"],
                   s=70, facecolors="none", edgecolors="black", linewidths=0.9,
                   alpha=0.85, zorder=9, **tf)
        ax.scatter(test.loc[is_event, "centroid_lon"], test.loc[is_event, "centroid_lat"],
                   s=62, marker="x", c="white", linewidths=3.2, zorder=10, **tf)
        ax.scatter(test.loc[is_event, "centroid_lon"], test.loc[is_event, "centroid_lat"],
                   s=55, marker="x", c="#b2182b", linewidths=1.6, zorder=11, **tf)
        axes[fold_no] = ax
        print("fold %d -> panel (%s): n=%d  daily-largest=%d  event-largest=%d"
              % (fold_no, chr(ord("a") + fold_no), len(test),
                 int(is_daily.sum()), int(is_event.sum())))

    # panel (a): distance box plots
    x, y = positions[0]
    axf = fig.add_axes([x / FIG_W, y / FIG_H, AX_W / FIG_W, AX_H / FIG_H])
    stats = []
    for label, col, sumcol in ZONE_COLUMNS:
        v = pc[col].to_numpy(float)
        m = v.mean()
        assert abs(m - float(s_all[sumcol])) < 1e-9, (label, m, s_all[sumcol])
        stats.append(dict(label=label, med=np.median(v), mean=m,
                          q1=np.percentile(v, 25), q3=np.percentile(v, 75),
                          whislo=v.min(), whishi=v.max(), fliers=[]))
        print("%s  mean=%.3f km  median=%.3f km  max=%.1f km" % (label, m, np.median(v), v.max()))
    arts = axf.bxp(stats, showmeans=True, showfliers=False, widths=0.5,
                   meanprops=dict(marker="D", markerfacecolor="white",
                                  markeredgecolor="black", markersize=5),
                   patch_artist=True)
    for box, med in zip(arts["boxes"], arts["medians"]):
        box.set_facecolor("0.88")
        box.set_edgecolor("0.25")
        med.set_color("0.15")
        med.set_linewidth(2.0)
    for a in arts["whiskers"] + arts["caps"]:
        a.set_color("0.25")
        a.set_linewidth(1.0)
    axf.set_ylabel("Distance to nearest top-concentration zone (km)", fontsize=9.5)
    axf.tick_params(axis="both", labelsize=8.2)
    axf.grid(True, axis="y", alpha=0.22, linewidth=0.45)
    axf.spines["top"].set_visible(False)
    axf.spines["right"].set_visible(False)
    axes[0] = axf

    fig.canvas.draw()
    for k, ax in axes.items():
        Q.label_above(fig, ax, "(%s)" % chr(ord("a") + k))
    sizes = {(round(a.get_position().width, 6), round(a.get_position().height, 6))
             for a in axes.values()}
    assert len(sizes) == 1, sizes

    cb_y0 = BOT_M + 0.55
    cb_h = 3 * AX_H + 2 * ROW_GAP
    cax = fig.add_axes([(FIG_W - CB_GUTTER + 0.18) / FIG_W, cb_y0 / FIG_H,
                        0.16 / FIG_W, cb_h / FIG_H])
    cbar = fig.colorbar(im, cax=cax, orientation="vertical", ticks=RANK_BOUNDS)
    cbar.set_label(CBAR_LABEL, fontsize=10)
    cbar.ax.tick_params(labelsize=9, length=2.5, width=0.6)
    cbar.outline.set_linewidth(0.6)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markeredgecolor="black", markersize=5, label="Centroid"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="black", markeredgewidth=1.0, markersize=9,
               label="Daily-largest structure centroid"),
        Line2D([0], [0], marker="x", color="#b2182b", markeredgewidth=1.6,
               markersize=7, linestyle="none", label="Event-largest structure centroid"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.012),
               ncol=3, fontsize=8.2, frameon=True, fancybox=False,
               framealpha=1.0, edgecolor="black", facecolor="white",
               handlelength=2.0, columnspacing=1.4, borderpad=0.65)

    save_figure(fig, FIGURES_DIR / "figure_C1", dpi=500,
                physical_width=TEXT_W_IN)
    plt.close(fig)

    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
