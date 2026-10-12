# SPDX-License-Identifier: GPL-3.0-only
"""The overview figure of ``scorch demo`` (matplotlib only, no map projections).

Top: the number of heatwave cells on every day of the synthetic record, the
regional threshold and the selected days coloured by the type of the
compound event they belong to. Bottom: one day of one event of each type,
showing the heatwave cells, their DBSCAN clusters, the fitted ellipses and
the temperature-weighted centroids.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from .. import config as cfg  # noqa: E402
from .. import ellipses as ell  # noqa: E402
from .. import events as ev  # noqa: E402

TYPE_COLORS = {1: "#D7191C", 2: "#F17C1E", 3: "#E6B400", 4: "#7030A0"}
CLUSTER_COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]


def representative_day(event, event_days):
    """The event day drawn in the figure: the first day with most clusters."""
    days = event_days[event_days["compound_event_id"] == event["compound_event_id"]]
    counts = [int(c) for c in str(event["daily_cluster_counts"]).split("-")]
    dates = days.sort_values("day_in_event")["date"].tolist()
    return dates[int(np.argmax(counts))]


def draw_day(ax, date, cluster_cells, ellipses, lat, lon):
    day_cells = cluster_cells[cluster_cells["date"] == date]
    day_ellipses = ellipses[ellipses["date"] == date]
    ax.set_xlim(lon[0] - 0.5, lon[-1] + 0.5)
    ax.set_ylim(lat[0] - 0.5, lat[-1] + 0.5)
    ax.set_xticks(np.arange(lon[0] - 0.5, lon[-1] + 1.0, 5.0))
    ax.set_yticks(np.arange(lat[0] - 0.5, lat[-1] + 1.0, 5.0))
    ax.grid(True, color="0.85", linewidth=0.5)
    ax.set_aspect(1.0 / np.cos(np.radians(lat.mean())))
    # Heatwave cells, coloured by cluster (grey for noise cells).
    for row in day_cells.itertuples():
        if pd.isna(row.cluster_id):
            colour = "0.6"
        else:
            number = int(day_ellipses.loc[day_ellipses["cluster_id"] == int(row.cluster_id),
                                          "cluster_number_in_day"].iloc[0])
            colour = CLUSTER_COLORS[number % len(CLUSTER_COLORS)]
        ax.add_patch(plt.Rectangle((row.lon - 0.5, row.lat - 0.5), 1.0, 1.0,
                                   facecolor=colour, edgecolor="white", linewidth=0.4, alpha=0.85))
    # Ellipses and weighted centroids.
    for row in day_ellipses.itertuples():
        geometry = {"centroid_lon_unweighted": row.centroid_lon_unweighted,
                    "centroid_lat_unweighted": row.centroid_lat_unweighted,
                    "major_axis_km": row.major_axis_km, "minor_axis_km": row.minor_axis_km,
                    "orientation_east_ccw_deg": row.orientation_east_ccw_deg}
        x, y = ell.ellipse_outline_lonlat(geometry)
        ax.plot(x, y, color="black", linewidth=1.6)
        ax.plot(row.centroid_lon_weighted, row.centroid_lat_weighted, "o", color="black",
                markersize=5, markeredgecolor="white")
    ax.set_xlabel("longitude (synthetic)", fontsize=8)
    ax.set_ylabel("latitude (synthetic)", fontsize=8)
    ax.tick_params(labelsize=7)


def draw_overview(config, directory, events, daily, selected, path):
    """Draw and save the overview figure; return its path."""
    catalogs = cfg.output_dir(config) / "catalogs"
    event_days = pd.read_csv(catalogs / "05_compound_heatwave_events" / "event_days.csv")
    cluster_cells = pd.read_csv(catalogs / "03_extreme_heatwave_clusters" / "cluster_cells.csv")
    ellipses = pd.read_csv(catalogs / "04_extreme_heatwave_ellipses" / "ellipses.csv")
    grid = pd.read_csv(cfg.supporting_input_dir(config) / "maps" / "grid_cells.csv")
    lat = np.sort(grid["lat"].unique())
    lon = np.sort(grid["lon"].unique())
    threshold = float(daily["regional_threshold_cells"].iloc[0])

    fig = plt.figure(figsize=(13.5, 9.5))
    grid_spec = GridSpec(2, 4, figure=fig, height_ratios=[1.0, 1.45], hspace=0.42, wspace=0.32,
                         left=0.06, right=0.98, top=0.90, bottom=0.10)

    # Top panel: daily heatwave-cell counts with the selected days.
    ax = fig.add_subplot(grid_spec[0, :])
    x = np.arange(len(daily))
    ax.fill_between(x, daily["n_heatwave_cells"], color="0.75", step="mid", linewidth=0)
    ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0)
    ax.text(len(daily) - 1, threshold + 2, f"regional threshold: {threshold:g} cells",
            ha="right", va="bottom", fontsize=8)
    type_of_date = (event_days.merge(events[["compound_event_id", "event_type"]],
                                     on="compound_event_id").set_index("date")["event_type"])
    selected_index = daily.index[daily["selected"] == 1]
    for t in (1, 2, 3, 4):
        dates_t = [i for i in selected_index if type_of_date.get(daily.loc[i, "date"]) == t]
        if dates_t:
            ax.plot(dates_t, daily.loc[dates_t, "n_heatwave_cells"], "o", color=TYPE_COLORS[t],
                    markersize=6, markeredgecolor="black", linewidth=0)
    years = pd.to_datetime(daily["date"]).dt.year
    top = float(daily["n_heatwave_cells"].max()) * 1.25
    ax.set_ylim(0, top)
    centres, labels = [], []
    for year in sorted(years.unique()):
        rows = np.where(years == year)[0]
        if rows[0] > 0:
            ax.axvline(rows[0] - 0.5, color="0.5", linewidth=0.8)
        centres.append(rows.mean())
        labels.append(f"{year} (Apr-Sep)")
    ax.set_xlim(0, len(daily) - 1)
    ax.set_xticks(centres)
    ax.set_xticklabels(labels, fontsize=9)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("cells in a heatwave")
    ax.set_title("Daily number of grid cells in a local heatwave; dots mark the selected days",
                 fontsize=10)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=TYPE_COLORS[t],
                      markeredgecolor="black", markersize=7, label=f"Type {t}: {ev.TYPE_NAMES[t]}")
               for t in (1, 2, 3, 4)]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True, ncol=2)

    # Bottom panels: one representative day per type.
    for column, t in enumerate((1, 2, 3, 4)):
        ax = fig.add_subplot(grid_spec[1, column])
        of_type = events[events["event_type"] == t]
        if of_type.empty:
            ax.text(0.5, 0.5, f"no Type {t} event", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        event = of_type.iloc[0]
        date = representative_day(event, event_days)
        draw_day(ax, date, cluster_cells, ellipses, lat, lon)
        span = (event["start_date"] if event["start_date"] == event["end_date"]
                else f"{event['start_date']} to {event['end_date']}")
        ax.set_title(f"Type {t}: {ev.TYPE_NAMES[t]}\nevent {int(event['compound_event_id'])}, "
                     f"{span}\nshown: {date}", fontsize=9, color=TYPE_COLORS[t], fontweight="bold")
    legend_handles = [Patch(facecolor=CLUSTER_COLORS[0], label="heatwave cells, cluster 1"),
                      Patch(facecolor=CLUSTER_COLORS[1], label="heatwave cells, cluster 2"),
                      Patch(facecolor="0.6", label="heatwave cells not in a cluster"),
                      Line2D([0], [0], color="black", linewidth=1.6, label="PCA ellipse"),
                      Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                             markeredgecolor="white", markersize=7,
                             label="temperature-weighted centroid")]
    fig.legend(handles=legend_handles, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 0.035))
    fig.suptitle("SCORCH demo on synthetic temperatures (not observations)", fontsize=14,
                 fontweight="bold")
    fig.text(0.5, 0.005, "Synthetic data generated by scorch.demo.synthetic; the coordinates do "
             "not describe a real place. Each map shows one day of one compound event.",
             ha="center", fontsize=8.5, color="0.35")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
