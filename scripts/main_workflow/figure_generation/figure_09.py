"""Figure 9: mean ellipse geometry for the 51 compound events.

Area and axis ratio use arithmetic means; orientation uses an axial mean
because directions 180 degrees apart describe the same axis. Panels show
interpolated empirical CDFs of mean axis ratio and area, and kernel densities
of mean orientation. All three panels share a Type 1–4 legend.

Writes the full PNG and its panels to output/figures/figure_9/.
Run: python scripts/main_workflow/figure_generation/figure_09.py
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

import math
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

import figure_style as Q
from figure_exports import save_figure  # noqa: E402


warnings.filterwarnings("ignore", category=RuntimeWarning)

TYPE_ORDER = ["Type 1", "Type 2", "Type 3", "Type 4"]
TYPE_COLORS = {
    "Type 1": "#d7191c",
    "Type 2": "#f57c00",
    "Type 3": "#d9a300",
    "Type 4": "#6a3d9a",
}

LINE_WIDTH = 2.15
LABEL_FS = 11.5
TICK_FS = 9.5
LEGEND_FS = 11
SPINE_LW = 0.8
GRID_KW = dict(alpha=0.22, linewidth=0.45)

# Layout in figure fractions: three data axes of identical size.
FIG_W, FIG_H = 15.0, 4.9
AX_W, AX_H = 0.28, 0.655
AX_BOTTOM = 0.225            # the region below is kept for the legend
AX_LEFT0, AX_GAP = 0.050, 0.0455
LEGEND_Y = 0.035


def orientation_from_north_deg(pc1_vec_x, pc1_vec_y):
    """Major-axis direction in degrees clockwise from north, in [-90, 90]."""
    vx = float(pc1_vec_x)
    vy = float(pc1_vec_y)
    if not np.isfinite(vx) or not np.isfinite(vy):
        return np.nan
    n = math.hypot(vx, vy)
    if n <= 0:
        return np.nan
    vx /= n
    vy /= n
    theta = np.degrees(np.arctan2(vx, vy))
    if theta > 90:
        theta -= 180
    if theta < -90:
        theta += 180
    return theta


def prepare_ellipse_dataframe(ellipses, compound_events):
    """One row per ellipse with its event, type and the three metrics."""
    events = compound_events[["compound_event_id", "event_type",
                              "event_type_name"]]
    df = ellipses.merge(events, on="compound_event_id", how="left",
                        validate="many_to_one")
    df = pd.DataFrame({
        "date": pd.to_datetime(df["date"]),
        "event_id": df["compound_event_id"],
        "event_type_name": df["event_type_name"],
        "type": "Type " + df["event_type"].astype(int).astype(str),
        "axis1_len_km": df["major_axis_km"].astype(float),
        "axis2_len_km": df["minor_axis_km"].astype(float),
        "pc1_vec_x": df["major_axis_east"].astype(float),
        "pc1_vec_y": df["major_axis_north"].astype(float),
    })
    df = df[(df["axis1_len_km"] > 0) & (df["axis2_len_km"] > 0)].copy()
    df = df[df["type"].isin(TYPE_ORDER)].copy()
    df["orientation_deg"] = df.apply(
        lambda r: orientation_from_north_deg(r["pc1_vec_x"], r["pc1_vec_y"]),
        axis=1)
    df["ellipse_area_km2"] = (
        np.pi * (df["axis1_len_km"] / 2.0) * (df["axis2_len_km"] / 2.0))
    df["ellipse_area_million_km2"] = df["ellipse_area_km2"] / 1e6
    df["ratio_L2_L1"] = df["axis2_len_km"] / df["axis1_len_km"]
    df = df[np.isfinite(df["orientation_deg"])
            & np.isfinite(df["ellipse_area_million_km2"])
            & np.isfinite(df["ratio_L2_L1"])].copy()
    return df


def convert_to_event_average_dataframe(df):
    """One row per compound event.

    Area and axis ratio: arithmetic means. Orientation: axial mean
    (functions.axial_mean_and_resultant), in the same clockwise-from-north
    convention, with its resultant length R (0 = no preferred axis,
    1 = all ellipses parallel).
    """
    keys = ["event_id", "event_type_name", "type"]
    event_df = (
        df.groupby(keys, as_index=False)
        .agg(
            start_date=("date", "min"),
            end_date=("date", "max"),
            n_ellipse_rows=("date", "size"),
            n_active_dates=("date", "nunique"),
            ellipse_area_million_km2=("ellipse_area_million_km2", "mean"),
            ratio_L2_L1=("ratio_L2_L1", "mean"),
        )
    )

    def _axial(group):
        values = group.to_numpy(dtype=float)
        mu, r = f.axial_mean_and_resultant(values)
        return pd.Series({"orientation_deg": mu,
                          "axial_resultant_R": r,
                          "n_orientation_values": int(np.isfinite(values).sum())})

    axial_stats = (
        df.groupby(keys)["orientation_deg"].apply(_axial).unstack()
        .reset_index())
    event_df = event_df.merge(axial_stats, on=keys, how="left",
                              validate="1:1")
    assert np.isfinite(event_df["orientation_deg"]).all()
    event_df["n_orientation_values"] = (
        event_df["n_orientation_values"].astype(int))
    event_df = event_df[
        ["event_id", "event_type_name", "type", "start_date", "end_date",
         "n_ellipse_rows", "n_active_dates", "orientation_deg",
         "ellipse_area_million_km2", "ratio_L2_L1",
         "axial_resultant_R", "n_orientation_values"]]
    event_df["event_id"] = event_df["event_id"].astype(int)
    return event_df


def smooth_empirical_cdf(values, x_min, x_max, n_grid=700):
    """Empirical CDF drawn as straight lines between its steps, anchored at
    0 at x_min and 1 at x_max."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    vals = np.sort(vals)
    if len(vals) == 0:
        return None, None
    y_emp = np.arange(1, len(vals) + 1) / len(vals)
    x_anchor = np.concatenate(([x_min], vals, [x_max]))
    y_anchor = np.concatenate(([0.0], y_emp, [1.0]))
    x_grid = np.linspace(x_min, x_max, n_grid)
    y_grid = np.interp(x_grid, x_anchor, y_anchor)
    return x_grid, y_grid


def get_dynamic_xlim(event_df, metric_col):
    """Data range plus 3.5 % on each side."""
    vals = event_df[metric_col].replace([np.inf, -np.inf], np.nan).dropna().values
    xmin = np.nanmin(vals)
    xmax = np.nanmax(vals)
    if xmin == xmax:
        pad = 0.05 if xmin == 0 else abs(xmin) * 0.05
    else:
        pad = 0.035 * (xmax - xmin)
    return xmin - pad, xmax + pad


def kde_curve(values, x_grid):
    """Gaussian kernel density (scipy default bandwidth, Scott's rule)."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 2 or np.nanstd(vals) == 0:
        return None
    return gaussian_kde(vals)(x_grid)


def _style_axes(ax):
    ax.grid(True, **GRID_KW)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(SPINE_LW)
    ax.spines["bottom"].set_linewidth(SPINE_LW)
    ax.tick_params(axis="both", labelsize=TICK_FS, length=3,
                   width=SPINE_LW, direction="out")


def draw_cdf(ax, event_df, metric_col, xlabel):
    x_min, x_max = get_dynamic_xlim(event_df, metric_col)
    for typ in TYPE_ORDER:
        sub = event_df[event_df["type"] == typ]
        values = sub[metric_col].replace([np.inf, -np.inf],
                                         np.nan).dropna().values
        if len(values) == 0:
            continue
        x_curve, y_curve = smooth_empirical_cdf(values, x_min, x_max)
        ax.plot(x_curve, y_curve, color=TYPE_COLORS[typ],
                linewidth=LINE_WIDTH, solid_capstyle="round",
                solid_joinstyle="round")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(xlabel, fontsize=LABEL_FS)
    ax.set_ylabel("Cumulative probability", fontsize=LABEL_FS)
    _style_axes(ax)


def draw_orientation(ax, event_df):
    x_grid = np.linspace(-90, 90, 600)
    for typ in TYPE_ORDER:
        vals = event_df.loc[event_df["type"] == typ,
                            "orientation_deg"].dropna().values
        density = kde_curve(vals, x_grid)
        if density is None:
            # Too few events for a density: the values are shown as points.
            ax.scatter(vals, np.zeros_like(vals), color=TYPE_COLORS[typ],
                       s=18, edgecolor="white", linewidth=0.35, zorder=4)
            continue
        ax.plot(x_grid, density, color=TYPE_COLORS[typ],
                linewidth=LINE_WIDTH, zorder=3)
    ax.axvline(0, color="0.35", linewidth=0.75, linestyle="--", zorder=1)
    ax.set_xlim(-90, 90)
    ax.set_xticks(np.arange(-90, 91, 30))
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Orientation (°)", fontsize=LABEL_FS)
    ax.set_ylabel("Probability Density", fontsize=LABEL_FS)
    _style_axes(ax)


# Panel order as in Figure 3: (a) axis ratio, (b) area, (c) orientation.
PANELS = [
    ("a", lambda ax, e: draw_cdf(ax, e, "ratio_L2_L1", r"Ratio ($L_2/L_1$)")),
    ("b", lambda ax, e: draw_cdf(ax, e, "ellipse_area_million_km2",
                                 "Area (km² × 10⁶)")),
    ("c", lambda ax, e: draw_orientation(ax, e)),
]


def make_figure_09(ellipses, compound_events, out_stem):
    """Draw Figure 9 and save its full PNG and panels."""
    plt.rcParams.update({
        "font.family": "Arial",
        # Mathtext in Arial so "Ratio ($L_2/L_1$)" reads as one label.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "it",
        "font.size": TICK_FS,
        "axes.labelsize": LABEL_FS,
        "xtick.labelsize": TICK_FS,
        "ytick.labelsize": TICK_FS,
        "axes.linewidth": SPINE_LW,
    })

    df_ellipse = prepare_ellipse_dataframe(ellipses, compound_events)
    event_df = convert_to_event_average_dataframe(df_ellipse)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    axes = []
    for i, (letter, drawer) in enumerate(PANELS):
        ax = fig.add_axes([AX_LEFT0 + i * (AX_W + AX_GAP), AX_BOTTOM,
                           AX_W, AX_H])
        drawer(ax, event_df)
        axes.append(ax)
    for (letter, _), ax in zip(PANELS, axes):
        Q.label_above(fig, ax, f"({letter})")
    fig.legend(handles=[Line2D([0], [0], color=TYPE_COLORS[t], lw=LINE_WIDTH,
                               label=t) for t in TYPE_ORDER],
               loc="lower center", bbox_to_anchor=(0.5, LEGEND_Y), ncol=4,
               frameon=False, fontsize=LEGEND_FS, handlelength=2.6,
               columnspacing=1.8)

    # Layout check: identical axes; the (c) label sits above orientation 0.
    fig.canvas.draw()
    sizes = {(round(a.get_position().width * FIG_W, 9),
              round(a.get_position().height * FIG_H, 9)) for a in axes}
    assert len(sizes) == 1, sizes
    ax_c = axes[2]
    x0_disp = ax_c.transData.transform((0.0, 0.0))[0]
    cx_disp = ((ax_c.get_position().x0 + ax_c.get_position().x1) / 2.0
               * FIG_W * fig.dpi)
    assert abs(x0_disp - cx_disp) < 0.51

    written = save_figure(fig, out_stem, dpi=600)
    plt.close(fig)

    return written

def main():
    config = f.load_config()
    inputs = {
        "ellipses": f.catalog_dir(config, "04_extreme_heatwave_ellipses") / "ellipses.csv",
        "events": f.catalog_dir(config, "05_compound_heatwave_events") / "compound_events.csv",
        "event_days": f.catalog_dir(config, "05_compound_heatwave_events") / "event_days.csv",
    }
    ellipses = f.attach_compound_events(pd.read_csv(inputs["ellipses"]),
                                       pd.read_csv(inputs["event_days"]))
    events = pd.read_csv(inputs["events"])
    out_stem = f.output_subdir(config, "figures") / "figure_9"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_09(ellipses, events, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
