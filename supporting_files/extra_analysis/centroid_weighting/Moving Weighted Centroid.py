#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare centroid weighting methods using animated maps.

Compare the unweighted centroid with Tmax, exceedance, percent-above-P95 and
hybrid weights. Each animation shows the tracks, weighting equation, daily shift,
and running mean and maximum shift for all heatwave cells on each day
(heatwave_id > 0).

Required CSV columns: date, lat1, lon1, val, thr_p95, heatwave_id.
Exceedance is calculated as val - thr_p95; any supplied exceed column is ignored.
The settings below retain paths from the original workspace.
"""

import os
import math
import shutil
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)

try:
    import imageio.v2 as imageio
except Exception:
    import imageio

# Optional map projection and geographic features.
HAVE_CARTOPY = False
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False


# Input, dates and drawing settings
MASTER_CSV = r"C:\Users\fawaw\OneDrive - University of Florida\Coding\Data\p95_exceed_heatwaves_1deg_flat\master_exceed_heatwaves_long.csv"
OUT_DIR    = r"C:\Users\fawaw\OneDrive - University of Florida\Coding\Outputs\najibi_meeting_centroid_gifs"

START_DATE = "2010-08-01"
END_DATE   = "2010-08-21"

MONTHS_IN_SCOPE = [4, 5, 6, 7, 8, 9]

# GIF settings
GIF_FPS = 0.002              # Animation frame rate.
GIF_LOOP = 0                 # infinite loop
FRAME_DPI = 180
OVERWRITE_FRAMES_DIR = False # Keep existing frame directories.

# Plot styling
GRIDLINE_STEP = 5
PAST_PT_SIZE = 18
CURR_PT_SIZE = 85
TRACK_LW = 2.2
CELL_PT_SIZE = 14
CELL_ALPHA = 0.28

# Arrow styling
ARROW_WIDTH = 0.0022
ARROW_HEADWIDTH = 5.5
ARROW_HEADLENGTH = 7.0
ARROW_HEADAxisLENGTH = 6.0

# Map padding
LON_PAD = 2.0
LAT_PAD = 2.0

SHOW_HW_CELLS = True

EPS = 1e-12


# Helper functions
def distance_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    if not (np.isfinite(lon1) and np.isfinite(lat1) and np.isfinite(lon2) and np.isfinite(lat2)):
        return np.nan
    mean_lat = 0.5 * (lat1 + lat2)
    dx = (lon1 - lon2) * math.cos(math.radians(mean_lat))
    dy = (lat1 - lat2)
    return math.sqrt(dx * dx + dy * dy) * 111.32


def safe_weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)

    good = np.isfinite(values) & np.isfinite(weights)
    if not np.any(good):
        return np.nan

    values = values[good]
    weights = weights[good]
    weights = np.where(weights > 0, weights, 0.0)

    sw = float(np.sum(weights))
    if sw <= EPS:
        return float(np.mean(values))
    return float(np.sum(values * weights) / sw)


def weighted_centroid(df_sub: pd.DataFrame, weight_col: str) -> Tuple[float, float]:
    lon = df_sub["lon1"].to_numpy(dtype=float)
    lat = df_sub["lat1"].to_numpy(dtype=float)
    w   = df_sub[weight_col].to_numpy(dtype=float)
    return safe_weighted_mean(lon, w), safe_weighted_mean(lat, w)


def ensure_frames_dir(path: str) -> None:
    if os.path.isdir(path):
        if OVERWRITE_FRAMES_DIR:
            shutil.rmtree(path)
            os.makedirs(path, exist_ok=True)
    else:
        os.makedirs(path, exist_ok=True)


def setup_axes(lon_min: float, lon_max: float, lat_min: float, lat_max: float):
    if HAVE_CARTOPY:
        fig, axes = plt.subplots(
            1, 2, figsize=(14.2, 6.8),
            subplot_kw={"projection": ccrs.PlateCarree()}
        )
    else:
        fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.8))

    for ax in axes:
        if HAVE_CARTOPY:
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND.with_scale("50m"), alpha=0.12)
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6)
            ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.6)
            ax.gridlines(
                crs=ccrs.PlateCarree(),
                draw_labels=False,
                linewidth=0.35,
                alpha=0.35,
                xlocs=np.arange(math.floor(lon_min), math.ceil(lon_max) + 1, GRIDLINE_STEP),
                ylocs=np.arange(math.floor(lat_min), math.ceil(lat_max) + 1, GRIDLINE_STEP),
            )
        else:
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_xticks(np.arange(math.floor(lon_min), math.ceil(lon_max) + 1, GRIDLINE_STEP))
            ax.set_yticks(np.arange(math.floor(lat_min), math.ceil(lat_max) + 1, GRIDLINE_STEP))
            ax.grid(True, alpha=0.30, linewidth=0.5)

        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")

    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.83, wspace=0.10)
    return fig, axes


def map_scatter(ax, x, y, **kwargs):
    if HAVE_CARTOPY:
        ax.scatter(x, y, transform=ccrs.PlateCarree(), **kwargs)
    else:
        ax.scatter(x, y, **kwargs)


def map_plot(ax, x, y, **kwargs):
    if HAVE_CARTOPY:
        ax.plot(x, y, transform=ccrs.PlateCarree(), **kwargs)
    else:
        ax.plot(x, y, **kwargs)


def map_quiver(ax, x, y, u, v, **kwargs):
    if HAVE_CARTOPY:
        ax.quiver(x, y, u, v, transform=ccrs.PlateCarree(), **kwargs)
    else:
        ax.quiver(x, y, u, v, **kwargs)


def draw_panel(ax, day_df: pd.DataFrame, track: List[Tuple[float, float]], title: str, color: str):
    if SHOW_HW_CELLS and not day_df.empty:
        map_scatter(
            ax,
            day_df["lon1"].values,
            day_df["lat1"].values,
            s=CELL_PT_SIZE,
            c="gray",
            alpha=CELL_ALPHA,
            linewidths=0,
            zorder=1
        )

    if len(track) >= 2:
        x = np.array([p[0] for p in track], dtype=float)
        y = np.array([p[1] for p in track], dtype=float)

        map_plot(ax, x, y, color=color, linewidth=TRACK_LW, alpha=0.95, zorder=3)

        u = x[1:] - x[:-1]
        v = y[1:] - y[:-1]
        good = np.isfinite(x[:-1]) & np.isfinite(y[:-1]) & np.isfinite(u) & np.isfinite(v)

        if np.any(good):
            map_quiver(
                ax,
                x[:-1][good],
                y[:-1][good],
                u[good],
                v[good],
                angles="xy",
                scale_units="xy",
                scale=1,
                width=ARROW_WIDTH,
                headwidth=ARROW_HEADWIDTH,
                headlength=ARROW_HEADLENGTH,
                headaxislength=ARROW_HEADAxisLENGTH,
                color=color,
                alpha=0.95,
                zorder=4
            )

    if len(track) >= 1:
        xc, yc = track[-1]
        map_scatter(
            ax, [xc], [yc],
            s=CURR_PT_SIZE, c=color, alpha=1.0,
            edgecolors="k", linewidths=0.7, zorder=6
        )

    ax.set_title(title, fontsize=12)

def make_weights(day_df: pd.DataFrame) -> pd.DataFrame:
    out = day_df.copy()

    T_max_i = pd.to_numeric(out["val"], errors="coerce")
    T95_i = pd.to_numeric(out["thr_p95"], errors="coerce")

    # Tmax
    out["w_tmax"] = np.maximum(T_max_i, 0.0)

    # True exceedance in °C
    out["w_exceed_true"] = np.maximum(T_max_i - T95_i, 0.0)

    # Percent above P95, clipped to [0, 5].
    out["w_pct_above_p95"] = np.maximum(100.0 * (T_max_i - T95_i) / T95_i, 0.0)
    out["w_pct_above_p95"] = np.clip(out["w_pct_above_p95"], 0.0, 5.0)

    # Hybrid
    out["w_hybrid"] = out["w_tmax"] * out["w_exceed_true"]

    return out


# Generate the comparisons
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(MASTER_CSV, usecols=["date", "lat1", "lon1", "val", "thr_p95", "heatwave_id"])
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["heatwave_id"] = pd.to_numeric(df["heatwave_id"], errors="coerce").fillna(0).astype(int)
    df["val"] = pd.to_numeric(df["val"], errors="coerce")
    df["thr_p95"] = pd.to_numeric(df["thr_p95"], errors="coerce")

    df = df.dropna(subset=["date", "lat1", "lon1", "val", "thr_p95"])
    df = df[df["date"].dt.month.isin(MONTHS_IN_SCOPE)].copy()

    hw = df[df["heatwave_id"] > 0].copy()
    if hw.empty:
        raise RuntimeError("No heatwave cells found in the requested CSV.")

    start_dt = pd.to_datetime(START_DATE)
    end_dt = pd.to_datetime(END_DATE)
    hw = hw[(hw["date"] >= start_dt) & (hw["date"] <= end_dt)].copy()

    days = pd.date_range(start_dt, end_dt, freq="D")

    lon_min = float(hw["lon1"].min() - LON_PAD)
    lon_max = float(hw["lon1"].max() + LON_PAD)
    lat_min = float(hw["lat1"].min() - LAT_PAD)
    lat_max = float(hw["lat1"].max() + LAT_PAD)

    methods: Dict[str, Dict[str, str]] = {
        "Tmax": {
            "weight_col": "w_tmax",
            "formula": (
                r"$w_i=T_{\max,i}$" + "\n" +
                r"$\bar{x}_w=\frac{\sum x_i\,w_i}{\sum w_i},\ \bar{y}_w=\frac{\sum y_i\,w_i}{\sum w_i}$"
            ),
            "subtitle": "Weight: absolute Tmax"
        },
        "Exceed": {
            "weight_col": "w_exceed_true",
            "formula": (
                r"$w_i=T_{\max,i}-T_{95,i}$" + "\n" +
                r"$\bar{x}_w=\frac{\sum x_i\,w_i}{\sum w_i},\ \bar{y}_w=\frac{\sum y_i\,w_i}{\sum w_i}$"
            ),
            "subtitle": "Weight: exceedance above P95 (°C)"
        },
        "PercentAboveP95": {
            "weight_col": "w_pct_above_p95",
            "formula": (
                r"$w_i=100\cdot\frac{T_{\max,i}-T_{95,i}}{T_{95,i}}$" + "\n" +
                r"$\mathrm{then\ capped\ at}\ 5$" + "\n" +
                r"$\bar{x}_w=\frac{\sum x_i\,w_i}{\sum w_i},\ \bar{y}_w=\frac{\sum y_i\,w_i}{\sum w_i}$"
            ),
            "subtitle": "Weight: percent above P95, capped at 5"
        },
        "Hybrid": {
            "weight_col": "w_hybrid",
            "formula": (
                r"$w_i=T_{\max,i}\cdot(T_{\max,i}-T_{95,i})$" + "\n" +
                r"$\bar{x}_w=\frac{\sum x_i\,w_i}{\sum w_i},\ \bar{y}_w=\frac{\sum y_i\,w_i}{\sum w_i}$"
            ),
            "subtitle": "Weight: Tmax × exceedance"
        }
    }

    summary_rows = []

    for method_name, meta in methods.items():
        weight_col = meta["weight_col"]
        formula = meta["formula"]
        subtitle = meta["subtitle"]

        frames_dir = os.path.join(OUT_DIR, f"{method_name}_frames")
        ensure_frames_dir(frames_dir)

        gif_path = os.path.join(OUT_DIR, f"{method_name}_vs_Original_{START_DATE}_to_{END_DATE}.gif")

        track_orig: List[Tuple[float, float]] = []
        track_meth: List[Tuple[float, float]] = []
        shifts: List[float] = []
        frame_paths: List[str] = []

        frame_idx = 0

        for day in days:
            day_df = hw[hw["date"] == day].copy()
            if day_df.empty:
                continue

            day_df = make_weights(day_df)

            lon_o = float(day_df["lon1"].mean())
            lat_o = float(day_df["lat1"].mean())

            lon_m, lat_m = weighted_centroid(day_df, weight_col)

            track_orig.append((lon_o, lat_o))
            track_meth.append((lon_m, lat_m))

            shift_today = distance_km(lon_o, lat_o, lon_m, lat_m)
            shifts.append(shift_today)
            shift_mean = float(np.nanmean(shifts))
            shift_max = float(np.nanmax(shifts))

            fig, axes = setup_axes(lon_min, lon_max, lat_min, lat_max)

            draw_panel(axes[0], day_df, track_orig, "Original centroid", color="black")
            draw_panel(axes[1], day_df, track_meth, f"{method_name}-weighted centroid", color="red")

            fig.suptitle(
                f"Centroid Comparison | {day.strftime('%Y-%m-%d')}\n{subtitle}",
                fontsize=14, y=0.96
            )

            # Weighting equation below the heading.
            fig.text(0.50, 0.835, formula, ha="center", va="center", fontsize=10.5)

            fig.text(
                0.50, 0.045,
                f"Position difference today: {shift_today:.1f} km    |    Running mean difference: {shift_mean:.1f} km    |    Running max difference: {shift_max:.1f} km",
                ha="center", va="center", fontsize=11
            )

            axes[0].text(
                0.5, 0.97,
                f"({lon_o:.3f}, {lat_o:.3f})",
                transform=axes[0].transAxes,
                ha="center", va="top", fontsize=10,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2)
            )

            axes[1].text(
                0.5, 0.97,
                f"({lon_m:.3f}, {lat_m:.3f})",
                transform=axes[1].transAxes,
                ha="center", va="top", fontsize=10,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2)
            )

            frame_path = os.path.join(frames_dir, f"frame_{frame_idx:03d}_{day.strftime('%Y-%m-%d')}.png")
            fig.savefig(frame_path, dpi=FRAME_DPI, bbox_inches="tight")
            plt.close(fig)

            frame_paths.append(frame_path)
            frame_idx += 1

            summary_rows.append({
                "method": method_name,
                "date": day.strftime("%Y-%m-%d"),
                "original_lon": lon_o,
                "original_lat": lat_o,
                "method_lon": lon_m,
                "method_lat": lat_m,
                "shift_km": shift_today,
                "weight_sum": float(np.nansum(day_df[weight_col].to_numpy(dtype=float)))
            })

        if not frame_paths:
            print(f"[SKIP] No frames created for {method_name}")
            continue

        duration_seconds = 1.0 / GIF_FPS
        images = [imageio.imread(fp) for fp in frame_paths]
        imageio.mimsave(gif_path, images, duration=duration_seconds, loop=GIF_LOOP)
        print(f"[SAVED] {gif_path}")

    summary_csv = os.path.join(OUT_DIR, f"centroid_method_comparison_summary_{START_DATE}_to_{END_DATE}.csv")
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)
    print(f"[SAVED] {summary_csv}")
    print(f"[DONE] All outputs are in one folder:\n{OUT_DIR}")


if __name__ == "__main__":
    main()