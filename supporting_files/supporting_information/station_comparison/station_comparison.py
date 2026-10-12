"""Compare Aswan station temperatures with the containing ERA5 cell and draw Figure S1.

Reads: supporting_files/supporting_information/station_comparison/input/ghcnd_aswan_EG000062414.csv,
processed ERA5 data, and catalogs.
Writes: supporting_files/supporting_information/station_comparison/results/ and supporting_files/supporting_information/figures/
Run: python supporting_files/supporting_information/station_comparison/station_comparison.py
"""
from __future__ import annotations

import io
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("SCORCH_REQUIRE_SLIDE_FONTS", "0")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.dates as mdates  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from PIL import Image  # noqa: E402

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402
import figure_style as Q  # noqa: E402
from figure_exports import save_figure, _panel  # noqa: E402
import snapshot_maps as SN  # noqa: E402

Image.MAX_IMAGE_PIXELS = None

INPUT_DIR = REPO_DIR / "supporting_files/supporting_information/station_comparison/input"
OUTPUT_DIR = REPO_DIR / "supporting_files/supporting_information/station_comparison/results"
FIGURES_DIR = REPO_DIR / "supporting_files" / "supporting_information" / "figures"

PANEL_LAYOUT = {
        "a": _panel((0, 0.06, 0.52, 0.486), header=(0, 0, 1, 0.06),
                    footer=(0, 0.486, 1, 0.53)),
        "b": _panel((0.52, 0.06, 1, 0.486), header=(0, 0, 1, 0.06),
                    footer=(0, 0.486, 1, 0.53)),
        "c": _panel((0, 0.55, 0.63, 1)),
        "d": _panel((0.63, 0.55, 1, 1)),
    }

# Station comparison settings.
STATION_ID = "EG000062414"
STATION_LAT, STATION_LON = 23.97, 32.78
STATION_FILE = "ghcnd_aswan_EG000062414.csv"
EVENT_ID = 25
EVENT_START, EVENT_END = "2016-06-07", "2016-06-08"
WINDOW_START, WINDOW_END = "2016-05-31", "2016-06-12"
THRESHOLD_PERCENTILE = 95
PERCENTILE_METHOD = "linear"       # numpy's default, as in the main analysis
DPI = 600

# Colours of the station panels.
COL_GHCND = "#1f1f1f"
COL_ERA5 = "#c44e52"


def save_csv(table, path, written):
    table.to_csv(path, index=False, lineterminator="\n")
    written.append(path)


def png_pixels(fig, **savefig_kwargs):
    """Save a figure as PNG in memory and read the pixels back (the manuscript
    figure was assembled from PNG files)."""
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", **savefig_kwargs)
    plt.close(fig)
    buffer.seek(0)
    return buffer


# Panels (a, b): the two days of Event 25
def draw_event_maps(days, cluster_cells, ellipses, compound_events):
    for date in days:
        row = compound_events[(compound_events["start_date"] <= date)
                              & (compound_events["end_date"] >= date)]
        assert len(row) == 1
        assert int(row["compound_event_id"].iloc[0]) == EVENT_ID
        assert int(row["event_type"].iloc[0]) == 3

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), squeeze=False,
                             subplot_kw={"projection": SN.PROJ})
    fig.subplots_adjust(top=0.80, bottom=0.14, left=0.055, right=0.985,
                        wspace=0.05)
    drawn = []
    for j, date in enumerate(days):
        ax = axes[0, j]
        for row in SN.render_day(ax, cluster_cells[cluster_cells["date"] == date],
                                 ellipses[ellipses["date"] == date]):
            drawn.append({"date": date, **row})
        Q.date_label(ax, date)
    assert len(drawn) == 4    # two clusters on each day
    Q.type_heading(fig, 0.5, 0.93, 3, fontsize=24)
    for ax, text in ((axes[0, 0], "t"), (axes[0, 1], "t+1")):
        b = ax.get_position()
        fig.text((b.x0 + b.x1) / 2.0, b.y1 + 0.012, rf"$\mathbf{{{text}}}$",
                 ha="center", va="bottom", fontsize=15, color="black")
    Q.unified_snapshot_legend(fig, subscript_axis_labels=True)
    return png_pixels(fig, dpi=DPI)


# Station comparison
def read_station_series(config):
    """Validate station Tmax and return degC values (NaN where TMAX is blank)."""
    obs = pd.read_csv(INPUT_DIR / STATION_FILE,
                      dtype={"STATION": str, "NAME": str, "DATE": str})
    if set(obs["STATION"]) != {STATION_ID}:
        raise ValueError(f"the station file should hold only {STATION_ID}")
    # Dates: every row a valid calendar date, in order, none repeated.
    obs["DATE"] = pd.to_datetime(obs["DATE"], format="%Y-%m-%d")
    n_duplicates = int(obs["DATE"].duplicated().sum())
    if n_duplicates or not obs["DATE"].is_monotonic_increasing:
        raise ValueError("station dates are repeated or out of order")
    # Units: the Climate Data Online export gives TMAX in whole degrees
    # Fahrenheit. Converted once; a file already in Celsius would show
    # decimals or an implausibly cool Aswan warm season.
    tmax_f = obs["TMAX"].dropna()
    if not (tmax_f == np.round(tmax_f)).all():
        raise ValueError("station TMAX is not in whole degrees Fahrenheit")
    obs["station_tmax_c"] = (obs["TMAX"] - 32) * 5 / 9
    warm = obs[obs["DATE"].dt.month.isin(config["study_area"]["warm_season_months"])]
    warm_mean_c = float(warm["station_tmax_c"].mean())
    if not (0 < obs["station_tmax_c"].min() and obs["station_tmax_c"].max() < 55
            and 30 < warm_mean_c < 45):
        raise ValueError("station temperatures are implausible for Aswan in degC")

    return obs[["DATE", "station_tmax_c"]]


def read_era5_cell(config):
    """Tmax of the 1-degree study cell that contains the station, at the
    stored precision, and the row of that cell in grid_cells.csv."""
    cells = f.read_grid_cells(config)
    # Strict inequalities: a station on a cell edge would stop the script.
    inside = cells[(cells["lat_south"] < STATION_LAT) & (STATION_LAT < cells["lat_north"])
                   & (cells["lon_west"] < STATION_LON) & (STATION_LON < cells["lon_east"])]
    if len(inside) != 1:
        raise ValueError(f"{len(inside)} grid cells contain the station, expected one")
    cell = inside.iloc[0]

    data = f.read_tmax_input(config)
    f.check_study_calendar(data["dates"], config)
    i, j = int(cell["lat_index"]), int(cell["lon_index"])
    if (data["lat"][i], data["lon"][j]) != (cell["lat"], cell["lon"]):
        raise ValueError("grid_cells.csv and tmax_processed.nc disagree on the cell centre")
    era5 = pd.DataFrame({"DATE": pd.to_datetime(data["dates"]),
                         "era5_tmax_c": data["tmax"][:, i, j].copy()})
    return era5, cell


def threshold_baseline(obs, era5, config):
    """The dates of both threshold calculations: April-September 1940-2025
    with a value in both series. A date missing in either series is left out
    of both; nothing is filled in."""
    area = config["study_area"]
    pairs = obs.merge(era5, on="DATE", how="inner")
    in_study_period = (pairs["DATE"].dt.month.isin(area["warm_season_months"])
                       & pairs["DATE"].dt.year.between(area["first_year"], area["last_year"]))
    pairs = pairs[in_study_period].dropna(subset=["station_tmax_c", "era5_tmax_c"])
    return pairs.sort_values("DATE").reset_index(drop=True)


def compare_station(obs, era5, pairs):
    # Both thresholds from the same paired dates, with the same method.
    station_p95 = float(np.percentile(pairs["station_tmax_c"], THRESHOLD_PERCENTILE,
                                      method=PERCENTILE_METHOD))
    era5_p95 = float(np.percentile(pairs["era5_tmax_c"], THRESHOLD_PERCENTILE,
                                   method=PERCENTILE_METHOD))

    # The window: every ERA5 day, with the station value where there is one.
    window = era5[(era5["DATE"] >= WINDOW_START) & (era5["DATE"] <= WINDOW_END)]
    window = window.merge(obs, on="DATE", how="left").reset_index(drop=True)
    has_station = window["station_tmax_c"].notna()
    window["used_in_metrics"] = (has_station & window["era5_tmax_c"].notna()).astype(int)
    window["event_day"] = ((window["DATE"] >= EVENT_START)
                           & (window["DATE"] <= EVENT_END)).astype(int)
    # Exceedance as in the main analysis: at or above the threshold.
    window["era5_at_or_above_p95"] = (window["era5_tmax_c"] >= era5_p95).astype(int)
    window["station_at_or_above_p95"] = ((window["station_tmax_c"] >= station_p95)
                                         .astype("Int64").mask(~has_station))
    both = window[window["used_in_metrics"] == 1]
    event = window[window["event_day"] == 1]

    station = both["station_tmax_c"].to_numpy(float)
    grid = both["era5_tmax_c"].to_numpy(float)
    stats = {
        "rmse": float(np.sqrt(np.mean((station - grid) ** 2))),
        "mae": float(np.mean(np.abs(station - grid))),
        "bias": float(np.mean(grid - station)),            # ERA5 minus station
        "r": float(both["station_tmax_c"].corr(both["era5_tmax_c"])),
        "station_p95": station_p95,
        "era5_p95": era5_p95,
        "n_baseline": len(pairs),
        "baseline_first": pairs["DATE"].min().strftime("%Y-%m-%d"),
        "baseline_last": pairs["DATE"].max().strftime("%Y-%m-%d"),
        "station_baseline_at_or_above": int((pairs["station_tmax_c"] >= station_p95).sum()),
        "era5_baseline_at_or_above": int((pairs["era5_tmax_c"] >= era5_p95).sum()),
        "station_event_days_at_or_above": int(event["station_at_or_above_p95"].sum()),
        "era5_event_days_at_or_above": int(event["era5_at_or_above_p95"].sum()),
        "n_event_days_with_station": int(event["station_tmax_c"].notna().sum()),
        "n_both": len(both),
        "n_era5_window": int(window["era5_tmax_c"].notna().sum()),
        "n_station_window": int(has_station.sum()),
    }
    return window, both, stats


STATION_RC = {
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.fontsize": 8.8,
}


def draw_timeseries(window, stats):
    """Panel (c): daily Tmax of the station and the ERA5 cell."""
    fig, ax = plt.subplots(figsize=(7.0, 3.9), dpi=300)

    # Station: straight segments join the observed days, also across days
    # without an observation (a plotting choice only; no value is created
    # for those days), with markers on the observed days.
    obs_valid = window.dropna(subset=["station_tmax_c"])
    ax.plot(obs_valid["DATE"], obs_valid["station_tmax_c"],
            color=COL_GHCND, linewidth=2.2, linestyle="-", alpha=0.98,
            zorder=5)
    ax.plot(obs_valid["DATE"], obs_valid["station_tmax_c"], linestyle="None",
            marker="s", markersize=4.6, markerfacecolor=COL_GHCND,
            markeredgecolor="white", markeredgewidth=0.5,
            label="GHCNd Tmax", zorder=6)

    ax.plot(window["DATE"], window["era5_tmax_c"], color=COL_ERA5,
            linewidth=2.0, linestyle="-", marker="s", markersize=4.2,
            markerfacecolor=COL_ERA5, markeredgecolor="white",
            markeredgewidth=0.45, alpha=0.96, label="ERA5 Tmax", zorder=4)

    ax.axhline(stats["station_p95"], color=COL_GHCND, linestyle="--",
               linewidth=1.5, alpha=0.75, label="GHCNd P95")
    ax.axhline(stats["era5_p95"], color=COL_ERA5, linestyle=":",
               linewidth=1.7, alpha=0.82, label="ERA5 P95")

    # The days of the regional compound event (Event 25), not a station
    # heatwave.
    ax.axvspan(pd.to_datetime(EVENT_START), pd.to_datetime(EVENT_END),
               color="red", alpha=0.12, zorder=0)

    ax.set_ylabel("Daily Tmax (°C)")
    ax.set_xlabel("Date")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.32)

    # The station legend entry shows both the line and the marker.
    handles, labels = ax.get_legend_handles_labels()
    handles[labels.index("GHCNd Tmax")] = Line2D(
        [0], [0], color=COL_GHCND, linewidth=2.2, linestyle="-", marker="s",
        markersize=4.6, markerfacecolor=COL_GHCND, markeredgecolor="white",
        markeredgewidth=0.5)
    ax.legend(handles, labels, frameon=False, loc="lower left",
              bbox_to_anchor=(0.08, 0.02))

    ax.set_xlim(pd.to_datetime(WINDOW_START), pd.to_datetime(WINDOW_END))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    fig.tight_layout()
    return png_pixels(fig, dpi=DPI, bbox_inches="tight")


def draw_scatter(both, stats):
    """Panel (d): ERA5 against station Tmax on the days with both values."""
    fig, ax = plt.subplots(figsize=(4.3, 4.1), dpi=300)
    ax.scatter(both["station_tmax_c"], both["era5_tmax_c"], s=42,
               linewidth=0.55, edgecolor="black", alpha=0.88)
    min_val = min(both["station_tmax_c"].min(), both["era5_tmax_c"].min()) - 1
    max_val = max(both["station_tmax_c"].max(), both["era5_tmax_c"].max()) + 1
    ax.plot([min_val, max_val], [min_val, max_val], linestyle="--",
            linewidth=1.4)
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_xlabel("GHCNd Tmax (°C)")
    ax.set_ylabel("ERA5 Tmax (°C)")
    text = (f"n = {stats['n_both']} days\n"
            f"r = {stats['r']:.2f}\n"
            f"RMSE = {stats['rmse']:.2f} °C\n"
            f"Bias = {stats['bias']:.2f} °C")
    ax.text(0.05, 0.95, text, transform=ax.transAxes, ha="left", va="top",
            fontsize=9, bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.88))
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.32)
    fig.tight_layout()
    return png_pixels(fig, dpi=DPI, bbox_inches="tight")


def top_spine_row(im, minfrac=0.30):
    """First image row with a long horizontal dark run (a panel's top spine)."""
    a = np.asarray(im, dtype=float) / 255.0
    lum = a[..., :3] @ np.array([0.299, 0.587, 0.114])
    dark = (lum < 0.35) & (a[..., 3] > 0.5)
    width = dark.shape[1]
    for r in range(dark.shape[0]):
        idx = np.flatnonzero(dark[r])
        if idx.size == 0:
            continue
        brk = np.where(np.diff(idx) > 1)[0]
        st = np.r_[idx[0], idx[brk + 1]]
        en = np.r_[idx[brk], idx[-1]]
        if any((b - a_ + 1) >= minfrac * width for a_, b in zip(st, en)):
            return r
    raise RuntimeError("no top spine found in a station panel")


def station_strip(timeseries_png, scatter_png, gap_px=60):
    """Panels (c) and (d) side by side at their native pixel size, aligned
    on their top spines, on a white background."""
    left = Image.open(timeseries_png).convert("RGBA")
    right = Image.open(scatter_png).convert("RGBA")
    ly, ry = top_spine_row(left), top_spine_row(right)
    off_l, off_r = max(0, ry - ly), max(0, ly - ry)
    w = left.size[0] + gap_px + right.size[0]
    h = max(left.size[1] + off_l, right.size[1] + off_r)
    canvas = Image.new("RGBA", (w, h), (255, 255, 255, 255))
    canvas.alpha_composite(left, (0, off_l))
    canvas.alpha_composite(right, (left.size[0] + gap_px, off_r))
    buffer = io.BytesIO()
    canvas.save(buffer, "PNG", optimize=False, compress_level=6)
    buffer.seek(0)
    return buffer


# Composite
def ink_row_blocks(im, thr=0.99):
    """Blocks of consecutive image rows that contain any non-white ink."""
    mask = np.any(im[..., :3] < thr, axis=2)
    if im.shape[2] == 4:
        mask &= im[..., 3] > 0.01
    rows = np.flatnonzero(mask.any(axis=1))
    blocks, start = [], rows[0]
    for i in range(1, len(rows)):
        if rows[i] != rows[i - 1] + 1:
            blocks.append((start, rows[i - 1]))
            start = rows[i]
    blocks.append((start, rows[-1]))
    return blocks


def frame_spans(im, minfrac=0.15):
    """x-spans of the panel frames: the long dark runs of the first image
    row that has them (the top spines)."""
    L = im[..., :3] @ np.array([0.299, 0.587, 0.114])
    dark = L < 0.35
    if im.shape[2] == 4:
        dark &= im[..., 3] > 0.5
    W = dark.shape[1]
    for r in range(dark.shape[0]):
        idx = np.flatnonzero(dark[r])
        if idx.size == 0:
            continue
        sp = np.where(np.diff(idx) > 1)[0]
        st = np.r_[idx[0], idx[sp + 1]]
        en = np.r_[idx[sp], idx[-1]]
        runs = [(int(a), int(b)) for a, b in zip(st, en)
                if b - a + 1 >= minfrac * W]
        if runs:
            return r, runs
    raise RuntimeError("no panel frame found")


FIG_W = 12.4            # design width in inches (the width of the maps)
MAPS_DPI = 600.0
LABEL_LINE_IN = 0.42    # line reserved for the panel letters
BLOCK_GAP_IN = 0.20     # gap between the maps and the station panels
# The Supplement text width (A4 with 1 inch margins: 9026 twips); the PNG is
# declared at the resolution that makes it exactly this wide.
TEXT_W_IN = 9026.0 / 1440.0


def compose_figure(maps_png, strip_png, out_stem):
    A = plt.imread(maps_png)
    B = plt.imread(strip_png)

    # Split the maps image between the "Type 3" heading and the t / t+1 line.
    ab = ink_row_blocks(A)
    split_a = (ab[0][1] + ab[1][0]) // 2
    A_top, A_rest = A[:split_a], A[split_a:]

    row_a, runs_a = frame_spans(A_rest)
    assert len(runs_a) == 2, runs_a
    cxa = [(a + b) / 2.0 for a, b in runs_a]
    row_b, runs_b = frame_spans(B)
    assert len(runs_b) == 2, runs_b
    cxb = [(a + b) / 2.0 for a, b in runs_b]

    h_top = A_top.shape[0] / MAPS_DPI
    h_rest = A_rest.shape[0] / MAPS_DPI
    h_b = FIG_W * B.shape[0] / B.shape[1]
    fig_h = (h_top + LABEL_LINE_IN + h_rest + BLOCK_GAP_IN + LABEL_LINE_IN
             + h_b + 0.06)
    fig = plt.figure(figsize=(FIG_W, fig_h))

    y = fig_h - h_top
    ax1 = fig.add_axes([0, y / fig_h, 1.0, h_top / fig_h])
    y -= LABEL_LINE_IN
    y_lab_ab = y / fig_h
    y -= h_rest
    ax2 = fig.add_axes([0, y / fig_h, 1.0, h_rest / fig_h])
    y -= BLOCK_GAP_IN + LABEL_LINE_IN
    y_lab_cd = y / fig_h
    y -= h_b
    ax3 = fig.add_axes([0, y / fig_h, 1.0, h_b / fig_h])
    for ax, im in ((ax1, A_top), (ax2, A_rest), (ax3, B)):
        ax.imshow(im, interpolation="bilinear")
        ax.set_axis_off()

    # Panel letters centred on each panel frame, on their own line.
    pad = Q._pt_to_fig(fig, Q.LABEL_PAD_PT, "y")
    for cx, letter in zip(cxa, ("(a)", "(b)")):
        fig.text(cx / A_rest.shape[1], y_lab_ab + pad, letter,
                 ha="center", va="bottom", **Q.PANEL_LABEL)
    for cx, letter in zip(cxb, ("(c)", "(d)")):
        fig.text(cx / B.shape[1], y_lab_cd + pad, letter,
                 ha="center", va="bottom", **Q.PANEL_LABEL)

    paths = save_figure(fig, out_stem, dpi=DPI, physical_width=TEXT_W_IN, panel_layout=PANEL_LAYOUT)
    plt.close(fig)
    return paths


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    config = f.load_config()
    catalogs = f.project_path(config["paths"]["output_dir"]) / "catalogs"
    cluster_cells = pd.read_csv(
        catalogs / "03_extreme_heatwave_clusters" / "cluster_cells.csv",
        usecols=["date", "lon", "lat", "dbscan_label"])
    ellipses = pd.read_csv(
        catalogs / "04_extreme_heatwave_ellipses" / "ellipses.csv")
    compound_events = pd.read_csv(
        catalogs / "05_compound_heatwave_events" / "compound_events.csv")
    days = config["representative_events"]["figure_s1_type_3"]
    written = []

    with mpl.rc_context():
        maps_png = draw_event_maps(days, cluster_cells, ellipses, compound_events)

    obs = read_station_series(config)
    era5, cell = read_era5_cell(config)
    pairs = threshold_baseline(obs, era5, config)
    window, both, stats = compare_station(obs, era5, pairs)

    with mpl.rc_context(STATION_RC):
        timeseries_png = draw_timeseries(window, stats)
        scatter_png = draw_scatter(both, stats)
    strip_png = station_strip(timeseries_png, scatter_png)

    with mpl.rc_context({"font.family": "Arial"}):
        written += compose_figure(maps_png, strip_png, FIGURES_DIR / "figure_S1")

    # The shared threshold baseline.
    baseline = pairs.assign(date=pairs["DATE"].dt.strftime("%Y-%m-%d"))
    save_csv(baseline[["date", "station_tmax_c", "era5_tmax_c"]],
             OUTPUT_DIR / "threshold_baseline_pairs.csv", written)

    daily = window.assign(date=window["DATE"].dt.strftime("%Y-%m-%d"))
    save_csv(daily[["date", "era5_tmax_c", "station_tmax_c", "used_in_metrics",
                    "event_day", "era5_at_or_above_p95", "station_at_or_above_p95"]],
             OUTPUT_DIR / "station_era5_daily_tmax.csv", written)

    common = {"baseline_first_date": stats["baseline_first"],
              "baseline_last_date": stats["baseline_last"],
              "n_baseline_days": stats["n_baseline"],
              "percentile": f"{THRESHOLD_PERCENTILE} (numpy, {PERCENTILE_METHOD})"}
    thresholds = pd.DataFrame([
        {"series": "station GHCND:" + STATION_ID, "cell_id": None,
         "lat": STATION_LAT, "lon": STATION_LON, **common,
         "warm_season_p95_tmax_c": stats["station_p95"],
         "n_baseline_days_at_or_above_p95": stats["station_baseline_at_or_above"]},
        {"series": "ERA5 1-degree study cell containing the station",
         "cell_id": int(cell["cell_id"]),
         "lat": float(cell["lat"]), "lon": float(cell["lon"]), **common,
         "warm_season_p95_tmax_c": stats["era5_p95"],
         "n_baseline_days_at_or_above_p95": stats["era5_baseline_at_or_above"]},
    ])
    thresholds["cell_id"] = thresholds["cell_id"].astype("Int64")   # blank for the station
    save_csv(thresholds, OUTPUT_DIR / "warm_season_thresholds.csv", written)

    metrics = pd.DataFrame([
        ("n_days_with_both_values", stats["n_both"], "days"),
        ("pearson_r", stats["r"], ""),
        ("rmse_c", stats["rmse"], "degC"),
        ("mae_c", stats["mae"], "degC"),
        ("mean_bias_era5_minus_station_c", stats["bias"], "degC"),
        ("station_warm_season_p95_c", stats["station_p95"], "degC"),
        ("era5_warm_season_p95_c", stats["era5_p95"], "degC"),
        ("n_threshold_baseline_days", stats["n_baseline"], "days"),
        ("station_event_days_at_or_above_p95",
         stats["station_event_days_at_or_above"], "days"),
        ("era5_event_days_at_or_above_p95",
         stats["era5_event_days_at_or_above"], "days"),
        ("n_event_days_with_station_value", stats["n_event_days_with_station"], "days"),
        ("n_era5_days_in_window", stats["n_era5_window"], "days"),
        ("n_station_days_in_window", stats["n_station_window"], "days"),
    ], columns=["metric", "value", "unit"],
        dtype=object)   # keeps the counts as whole numbers
    save_csv(metrics, OUTPUT_DIR / "station_comparison_metrics.csv", written)

    print(f"ERA5 cell {int(cell['cell_id'])} ({cell['lat']} N, {cell['lon']} E); "
          f"baseline {stats['n_baseline']} dates {stats['baseline_first']} to "
          f"{stats['baseline_last']}; P95 station {stats['station_p95']:.4f}, "
          f"ERA5 {stats['era5_p95']:.4f} degC")
    print(f"window: n = {stats['n_both']}, r = {stats['r']:.4f}, "
          f"RMSE = {stats['rmse']:.4f} degC, bias = {stats['bias']:.4f} degC; "
          f"event days at or above P95: station {stats['station_event_days_at_or_above']}, "
          f"ERA5 {stats['era5_event_days_at_or_above']}")
    for path in written:
        print("wrote", f.relative_to_repo(path))


if __name__ == "__main__":
    main()
