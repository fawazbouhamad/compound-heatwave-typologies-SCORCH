"""Figure 10: annual event counts and mean duration.

Panel (a) shows event counts by type, excluding years without events.
Panel (b) shows mean Type 3 and Type 4 duration by event start year.
Both use step-5 tables. Duration trend tests are saved separately in
supporting_files/additional_results/main_workflow/duration_trend_statistics.csv.

Writes the full PNG and its panels to output/figures/figure_10/.
Run: python scripts/main_workflow/figure_generation/figure_10.py
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
from figure_exports import save_figure

TYPE_ORDER = ["Type 1", "Type 2", "Type 3", "Type 4"]
DURATION_TYPES = ["Type 3", "Type 4"]
COLORS = {
    "Type 1": "#d7191c",
    "Type 2": "#f57c00",
    "Type 3": "#d9a300",
    "Type 4": "#6a3d9a",
}
FIRST_YEAR, LAST_YEAR = 1983, 2025      # years shown on the time axis
YEAR_TICKS = [1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025]


def make_figure_10(annual_counts_table, annual_duration_table,
                   compound_events, out_stem):
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.75,
        "ytick.major.width": 0.75,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 600,
    })

    # Panel (a): events per year and type; years without events left out.
    years = np.arange(FIRST_YEAR, LAST_YEAR + 1)
    annual_counts = (annual_counts_table.set_index("year")
                     [["type_1_events", "type_2_events", "type_3_events",
                       "type_4_events"]]
                     .set_axis(TYPE_ORDER, axis=1)
                     .reindex(index=years, fill_value=0))
    annual_plot = annual_counts[annual_counts.sum(axis=1) > 0]

    # Panel (b): mean duration per start year, Type 3 and Type 4.
    duration_summary = annual_duration_table.assign(
        type="Type " + annual_duration_table["event_type"].astype(str))
    duration_summary = duration_summary[
        duration_summary["type"].isin(DURATION_TYPES)]
    longest_duration = compound_events.loc[
        compound_events["event_type"].isin([3, 4]), "duration_days"].max()

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9.4, 6.45),
        gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.33})

    bottom = np.zeros(len(annual_plot))
    for typ in TYPE_ORDER:
        vals = annual_plot[typ].values
        ax1.bar(annual_plot.index, vals, bottom=bottom, width=0.72,
                color=COLORS[typ], edgecolor="white", linewidth=0.45,
                alpha=0.95, label=f"{typ}", zorder=3)
        bottom += vals

    ax1.set_xlim(1982.5, 2025.5)
    ax1.set_ylim(0, annual_plot.sum(axis=1).max() + 1)
    ax1.set_ylabel("Number of events")
    ax1.set_xticks(YEAR_TICKS)
    ax1.set_yticks(np.arange(0, annual_plot.sum(axis=1).max() + 2, 1))
    ax1.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.22)
    ax1.set_axisbelow(True)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.tick_params(axis="both", labelsize=8.5)
    ax1.legend(loc="upper left", frameon=False, fontsize=7.5, ncol=2,
               handlelength=1.2, columnspacing=1.1, handletextpad=0.45)
    # Panel letter to the left of the y axis, vertically centred; the same
    # offset for both panels.
    ax1.text(-0.135, 0.5, "(a)", transform=ax1.transAxes, ha="right",
             va="center", fontsize=17, fontweight="bold", color="black",
             family="Arial")

    for typ in DURATION_TYPES:
        sub = duration_summary[duration_summary["type"] == typ].sort_values(
            "year")
        ax2.plot(sub["year"], sub["mean_duration_days"], color=COLORS[typ],
                 linewidth=2.0, marker="o", markersize=3.6,
                 markeredgewidth=0.0, solid_capstyle="round",
                 label=f"{typ}", zorder=3)

    ax2.set_xlim(1983, 2025)
    ax2.set_ylim(0, longest_duration + 4)
    ax2.set_ylabel("Mean duration (days)")
    ax2.set_xlabel("Year")
    ax2.set_xticks(YEAR_TICKS)
    ax2.set_yticks(np.arange(0, 51, 5))
    ax2.grid(axis="y", linestyle="--", linewidth=0.40, alpha=0.18)
    ax2.set_axisbelow(True)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.tick_params(axis="both", labelsize=8.5, length=3.2, width=0.75)
    ax2.legend(loc="upper left", frameon=False, fontsize=7.5, ncol=2,
               handlelength=2.0, columnspacing=1.2, handletextpad=0.45)
    ax2.text(-0.135, 0.5, "(b)", transform=ax2.transAxes, ha="right",
             va="center", fontsize=17, fontweight="bold", color="black",
             family="Arial")

    fig.align_ylabels([ax1, ax2])

    written = save_figure(fig, out_stem, bbox_inches="tight")
    plt.close(fig)

    return written

def main():
    config = f.load_config()
    inputs = {
        "counts": f.supporting_results_dir(config) / "main_workflow/annual_event_counts.csv",
        "durations": f.supporting_results_dir(config) / "main_workflow/annual_mean_duration.csv",
        "events": f.catalog_dir(config, "05_compound_heatwave_events") / "compound_events.csv",
    }
    counts = pd.read_csv(inputs["counts"])
    durations = pd.read_csv(inputs["durations"])
    events = pd.read_csv(inputs["events"])
    out_stem = f.output_subdir(config, "figures") / "figure_10"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = make_figure_10(counts, durations, events, out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
