"""Draw one day of a compound event for Figures 5–7 and S1.

Heatwave cells are orange for DBSCAN clusters and grey for noise. Each cluster
has an unweighted PCA ellipse (scale 1.25), its axes and a Tmax-weighted centroid.
The shape is recomputed with the step-4 functions and centred on the weighted
centroid from the ellipse catalog.

The caller supplies rows from cluster_cells.csv and ellipses.csv. This module
draws on the supplied axes and does not write files.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import functions as f  # noqa: E402

PROJ = ccrs.PlateCarree()

# Drawing style of the snapshot maps.
COL_CLUSTERED = "orange"
COL_NOISE = "0.55"
COL_ELLIPSE = "black"
SCALE_FACTOR = 1.25            # manuscript ellipse scale, also recorded in functions.py
EIGENVALUE_FLOOR_KM2 = 1.0     # PCA eigenvalue floor, also recorded in functions.py
N_ELLIPSE_POINTS = 240
# Fixed map extent (the study area plus 0.5 degree) so that all maps have the
# same frame.
STUDY_EXTENT = (19.5, 70.5, 9.5, 46.5)   # lon_min, lon_max, lat_min, lat_max


def setup_ax(ax, extent=STUDY_EXTENT):
    lon_min, lon_max, lat_min, lat_max = extent
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=PROJ)
    ax.add_feature(cfeature.LAND.with_scale("50m"), alpha=0.12, zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.55,
                   zorder=8)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.45,
                   zorder=8)
    gl = ax.gridlines(crs=PROJ, draw_labels=True, linewidth=0.28,
                      alpha=0.18, linestyle="-")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 7}
    gl.ylabel_style = {"size": 7}


def _draw_cells(ax, cells, facecolor, alpha, lw, zorder):
    for lon, lat in cells:
        ax.add_patch(Rectangle((lon - 0.5, lat - 0.5), 1.0, 1.0,
                               facecolor=facecolor, edgecolor="k",
                               linewidth=lw, alpha=alpha, zorder=zorder,
                               transform=PROJ))


def _km_to_lonlat(center_lon, center_lat, x_km, y_km):
    """Inverse of the local projection of functions.project_lonlat_to_km."""
    km_per_deg_lon = f.KM_PER_DEG_LON_AT_EQUATOR * math.cos(
        math.radians(center_lat))
    return (center_lon + x_km / km_per_deg_lon,
            center_lat + y_km / f.KM_PER_DEG_LAT)


def ellipse_outline_lonlat(center_lon, center_lat, eigvecs, eigvals):
    """240 points of the ellipse outline around (center_lon, center_lat)."""
    a, b = f.ellipse_semi_axes_km(eigvals, SCALE_FACTOR, EIGENVALUE_FLOOR_KM2)
    t = np.linspace(0.0, 2.0 * np.pi, N_ELLIPSE_POINTS)
    circle = np.vstack([a * np.cos(t), b * np.sin(t)])
    ell = np.zeros(2).reshape(2, 1) + eigvecs @ circle
    return _km_to_lonlat(center_lon, center_lat, ell[0, :], ell[1, :])


def axis_endpoints_lonlat(center_lon, center_lat, eigvecs, eigvals):
    """End points of the major axis (L1) and the minor axis (L2)."""
    a, b = f.ellipse_semi_axes_km(eigvals, SCALE_FACTOR, EIGENVALUE_FLOOR_KM2)
    v1 = eigvecs[:, 0] * a
    v2 = eigvecs[:, 1] * b
    c = np.zeros(2)
    return tuple(_km_to_lonlat(center_lon, center_lat, p[0], p[1])
                 for p in (c - v1, c + v1, c - v2, c + v2))


def render_day(ax, day_cells, day_ellipses):
    """Draw one day on a map axes.

    day_cells: the rows of cluster_cells.csv for this date, in file order
    (columns lon, lat, dbscan_label; -1 = noise).
    day_ellipses: the rows of ellipses.csv for this date.
    Returns one dict per cluster with the drawn centroid and ellipse size.
    """
    setup_ax(ax)
    lon = day_cells["lon"].to_numpy(float)
    lat = day_cells["lat"].to_numpy(float)
    labels = day_cells["dbscan_label"].to_numpy(int)
    clustered = [(round(float(lon[i]), 6), round(float(lat[i]), 6))
                 for i in range(lon.size) if labels[i] >= 0]
    noise = [(round(float(lon[i]), 6), round(float(lat[i]), 6))
             for i in range(lon.size) if labels[i] < 0]
    if noise:
        _draw_cells(ax, noise, COL_NOISE, 0.70, 0.20, 3)
    if clustered:
        _draw_cells(ax, clustered, COL_CLUSTERED, 0.92, 0.24, 5)

    drawn = []
    for label in sorted(set(labels.tolist()) - {f.NOISE}):
        member = labels == label
        geometry = f.ellipse_geometry(lon[member], lat[member], SCALE_FACTOR,
                                      EIGENVALUE_FLOOR_KM2)
        # The catalog numbers the clusters of a day by their DBSCAN label.
        row = day_ellipses[day_ellipses["cluster_number_in_day"] == label]
        assert len(row) == 1
        # Same cluster and same ellipse as in the catalog. The member cells
        # are used in file order, as in the manuscript maps; step 4 sorts them
        # by longitude first, so the PCA sums run in a different order and
        # the axis lengths can differ in the last binary digit (relative
        # difference at most 2.4e-16 over the 760 ellipses of the historical analysis).
        assert int(row["n_member_cells"].iloc[0]) == int(member.sum())
        assert np.isclose(float(row["major_axis_km"].iloc[0]),
                          geometry["major_axis_km"], rtol=1e-12, atol=0.0)
        assert np.isclose(float(row["minor_axis_km"].iloc[0]),
                          geometry["minor_axis_km"], rtol=1e-12, atol=0.0)
        wlon = float(row["centroid_lon_weighted"].iloc[0])
        wlat = float(row["centroid_lat_weighted"].iloc[0])
        eigvecs, eigvals = geometry["_eigvecs"], geometry["_eigvals"]

        elon, elat = ellipse_outline_lonlat(wlon, wlat, eigvecs, eigvals)
        ax.plot(elon, elat, lw=2.0, color=COL_ELLIPSE, zorder=7,
                transform=PROJ)
        p1a, p1b, p2a, p2b = axis_endpoints_lonlat(wlon, wlat, eigvecs, eigvals)
        ax.plot([p1a[0], p1b[0]], [p1a[1], p1b[1]], lw=2.4, color="black",
                zorder=7, transform=PROJ)
        ax.plot([p2a[0], p2b[0]], [p2a[1], p2b[1]], lw=1.5, color="black",
                zorder=7, transform=PROJ)
        ax.scatter([wlon], [wlat], s=60, marker="o", c="black",
                   edgecolors="k", linewidths=0.7, zorder=8, transform=PROJ)
        drawn.append(dict(ellipse_id=int(row["ellipse_id"].iloc[0]),
                          dbscan_label=int(label),
                          n_member_cells=int(member.sum()),
                          centroid_lon_unweighted=geometry[
                              "centroid_lon_unweighted"],
                          centroid_lat_unweighted=geometry[
                              "centroid_lat_unweighted"],
                          centroid_lon_weighted=wlon,
                          centroid_lat_weighted=wlat,
                          major_axis_km=geometry["major_axis_km"],
                          minor_axis_km=geometry["minor_axis_km"],
                          area_km2=geometry["area_km2"]))
    return drawn
