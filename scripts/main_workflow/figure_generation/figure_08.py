"""Figure 8: ellipse geometry by compound-event type.

Histograms and box plots show axis ratio, area (million km2) and major-axis
orientation (degrees clockwise from north) for the 753 ellipses. Each column
is one event type. Histograms show the sample size and median.

The crop is measured with a plain-text ratio label before drawing its mathtext
version, to keep the manuscript panel positions.

Writes the full PNG and its panels to output/figures/figure_8/.
Run: python scripts/main_workflow/figure_generation/figure_08.py
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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from figure_exports import save_figure

TYPE_ORDER = ["Type 1", "Type 2", "Type 3", "Type 4"]

TYPE_COLORS = {
    "Type 1": "#d7191c",
    "Type 2": "#f57c00",
    "Type 3": "#d9a300",
    "Type 4": "#6a3d9a",
}

FIG_W = 12.2
FIG_H = 9.8
DPI = 500

LABEL_FS = 9.5
TICK_FS = 7.5
TITLE_FS = 11
ANNOT_FS = 6.8

HIST_ALPHA = 0.90
EDGE_LW = 0.35

BOX_FACE = "#E6E6E6"
BOX_EDGE = "#404040"
MEDIAN_COLOR = "#8B1E1E"

ORIENTATION_BINS = np.arange(-90, 90 + 15, 15)
RATIO_BINS = np.arange(0.0, 1.0001, 0.05)
AREA_BIN_COUNT = 12

RATIO_LABEL = r"Ratio ($L_2/L_1$)"
RATIO_LABEL_FOR_CROP = "Ratio (L2/L1)"


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


def summarise(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return {"n": 0, "mean": np.nan, "std": np.nan, "min": np.nan,
                "q1": np.nan, "median": np.nan, "q3": np.nan, "max": np.nan}

    return {
        "n": len(values),
        "mean": np.mean(values),
        "std": np.std(values, ddof=1) if len(values) > 1 else np.nan,
        "min": np.min(values),
        "q1": np.percentile(values, 25),
        "median": np.median(values),
        "q3": np.percentile(values, 75),
        "max": np.max(values),
    }


def prepare_dataframe(ellipses, compound_events):
    """One row per ellipse with its event, type and the three metrics."""
    types = compound_events[["compound_event_id", "event_type"]]
    df = ellipses.merge(types, on="compound_event_id", how="left",
                        validate="many_to_one")
    df = pd.DataFrame({
        "ellipse_id": df["ellipse_id"],
        "date": pd.to_datetime(df["date"]),
        "event_id": df["compound_event_id"],
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


def nice_spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)


def draw_hist_box(ax_hist, ax_box, values, bins, color, xlim=None):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        ax_hist.text(0.5, 0.5, "No data", transform=ax_hist.transAxes,
                     ha="center", va="center", fontsize=8)
        ax_box.axis("off")
        return

    ax_hist.hist(vals, bins=bins, color=color, alpha=HIST_ALPHA,
                 edgecolor="black", linewidth=EDGE_LW)

    ax_hist.grid(axis="y", alpha=0.22, linewidth=0.45)
    ax_hist.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_hist.tick_params(axis="both", labelsize=TICK_FS, length=2.5)
    nice_spines(ax_hist)

    s = summarise(vals)
    txt = f"n={s['n']}\nMed={s['median']:.2f}"

    ax_hist.text(0.96, 0.92, txt, transform=ax_hist.transAxes,
                 ha="right", va="top", fontsize=ANNOT_FS,
                 bbox=dict(boxstyle="round,pad=0.20", facecolor="white",
                           edgecolor="0.75", alpha=0.95))

    ax_box.boxplot(
        vals,
        vert=False,
        widths=0.55,
        patch_artist=True,
        boxprops=dict(facecolor=BOX_FACE, edgecolor=BOX_EDGE, linewidth=0.75),
        medianprops=dict(color=MEDIAN_COLOR, linewidth=1.25),
        whiskerprops=dict(color=BOX_EDGE, linewidth=0.75),
        capprops=dict(color=BOX_EDGE, linewidth=0.75),
        flierprops=dict(marker="o", markerfacecolor="white",
                        markeredgecolor=BOX_EDGE, markersize=2.2,
                        alpha=0.85),
    )

    ax_box.set_yticks([])
    ax_box.grid(axis="x", alpha=0.18, linewidth=0.45)
    ax_box.tick_params(axis="x", labelsize=TICK_FS, length=2.5)
    nice_spines(ax_box)

    if xlim is not None:
        ax_hist.set_xlim(*xlim)
        ax_box.set_xlim(*xlim)

    plt.setp(ax_hist.get_xticklabels(), visible=False)


def make_figure_08(ellipses, compound_events, out_stem):
    """Draw Figure 8 and save its full PNG and panels."""
    plt.rcParams.update({
        "font.family": "serif",
        # Mathtext in the serif body face so "Ratio ($L_2/L_1$)" reads as
        # one label.
        "mathtext.fontset": "dejavuserif",
        "font.size": 9,
        "axes.labelsize": LABEL_FS,
        "axes.titlesize": TITLE_FS,
        "xtick.labelsize": TICK_FS,
        "ytick.labelsize": TICK_FS,
    })

    df = prepare_dataframe(ellipses, compound_events)

    area_min = np.nanmin(df["ellipse_area_million_km2"])
    area_max = np.nanmax(df["ellipse_area_million_km2"])
    if area_min == area_max:
        area_bins = np.linspace(area_min - 0.1, area_max + 0.1,
                                AREA_BIN_COUNT + 1)
    else:
        area_bins = np.linspace(area_min, area_max, AREA_BIN_COUNT + 1)

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
    outer = fig.add_gridspec(3, 4, hspace=0.46, wspace=0.28, left=0.075,
                             right=0.995, bottom=0.075, top=0.965)

    # Rows in the order of Figure 3: (a) axis ratio, (b) area,
    # (c) orientation.
    row_info = [
        {"metric": "ratio_L2_L1", "xlabel": RATIO_LABEL,
         "bins": RATIO_BINS, "xlim": (0, 1)},
        {"metric": "ellipse_area_million_km2", "xlabel": "Area (km² × 10⁶)",
         "bins": area_bins, "xlim": None},
        {"metric": "orientation_deg", "xlabel": "Orientation (°)",
         "bins": ORIENTATION_BINS, "xlim": (-90, 90)},
    ]

    for r, info in enumerate(row_info):
        for c, typ in enumerate(TYPE_ORDER):
            subgs = outer[r, c].subgridspec(2, 1, height_ratios=[4.1, 1.0],
                                            hspace=0.03)
            ax_hist = fig.add_subplot(subgs[0, 0])
            ax_box = fig.add_subplot(subgs[1, 0], sharex=ax_hist)

            values = df.loc[df["type"] == typ, info["metric"]].values
            draw_hist_box(ax_hist=ax_hist, ax_box=ax_box, values=values,
                          bins=info["bins"], color=TYPE_COLORS[typ],
                          xlim=info["xlim"])
            ax_box.set_xlabel(info["xlabel"], fontsize=LABEL_FS, labelpad=4)

            if r == 0:
                ax_hist.set_title(typ, fontsize=TITLE_FS, fontweight="bold",
                                  pad=7, color=TYPE_COLORS[typ])
            if c == 0:
                ax_hist.set_ylabel("Number of ellipses", fontsize=LABEL_FS)
            else:
                ax_hist.set_ylabel("")
            ax_hist.tick_params(axis="y", labelleft=True)

    # Row letters (a)-(c): Arial bold 17 pt, one shared x position, each
    # centred on the vertical span of its row (rows found from the two
    # largest vertical gaps between the axes).
    boxes = [a.get_position() for a in fig.axes]
    ycs = sorted((((b.y0 + b.y1) / 2.0) for b in boxes), reverse=True)
    gaps = sorted(range(len(ycs) - 1), key=lambda i: ycs[i] - ycs[i + 1],
                  reverse=True)[:2]
    cuts = sorted(((ycs[i] + ycs[i + 1]) / 2.0 for i in gaps), reverse=True)
    rows = [[], [], []]
    for b in boxes:
        yc = (b.y0 + b.y1) / 2.0
        rows[0 if yc > cuts[0] else (1 if yc > cuts[1] else 2)].append(b)
    x_left = min(b.x0 for b in boxes) - 0.052   # clear of the y-axis titles
    for letter, row in zip(("(a)", "(b)", "(c)"), rows):
        yc = (min(b.y0 for b in row) + max(b.y1 for b in row)) / 2.0
        fig.text(x_left, yc, letter, ha="right", va="center", fontsize=17,
                 fontweight="bold", color="black", family="Arial")

    # Measure the crop box with the plain-text ratio label (see the module
    # description), then draw the mathtext label.
    ratio_axes = [a for a in fig.axes if a.get_xlabel() == RATIO_LABEL]
    assert len(ratio_axes) == len(TYPE_ORDER)

    def set_ratio_labels(text):
        for a in ratio_axes:
            a.set_xlabel(text, fontsize=LABEL_FS, labelpad=4)

    original_get_tightbbox = fig.get_tightbbox

    def get_tightbbox_with_plain_label(renderer=None, *args, **kwargs):
        set_ratio_labels(RATIO_LABEL_FOR_CROP)
        try:
            return original_get_tightbbox(renderer, *args, **kwargs)
        finally:
            set_ratio_labels(RATIO_LABEL)

    fig.get_tightbbox = get_tightbbox_with_plain_label
    written = save_figure(fig, out_stem, dpi=DPI, bbox_inches="tight")
    fig.get_tightbbox = original_get_tightbbox
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
    out_stem = f.output_subdir(config, "figures") / "figure_8"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_08(ellipses, events, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
