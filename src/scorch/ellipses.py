# SPDX-License-Identifier: GPL-3.0-only
"""Stage 4: PCA ellipse geometry and temperature-weighted centroid.

Each cluster's member cells are projected to kilometres around the cluster's
mean position. The sample covariance of the projected coordinates gives the
ellipse: semi-axes ``s * sqrt(max(eigenvalue, floor))`` along the two
eigenvectors, with ``s`` the scale factor. The temperature-weighted centroid
is the Tmax-weighted mean position of the member cells.
"""
from __future__ import annotations

import math

import numpy as np

KM_PER_DEG_LON_AT_EQUATOR = 111.320
KM_PER_DEG_LAT = 110.574


def project_lonlat_to_km(lon, lat, lon0, lat0):
    """Local equirectangular projection to km around (lon0, lat0).

    One degree of latitude is 110.574 km and one degree of longitude is
    111.320 km times cos(lat0).
    """
    km_per_deg_lon = KM_PER_DEG_LON_AT_EQUATOR * math.cos(math.radians(lat0))
    x = (np.asarray(lon, dtype=float) - lon0) * km_per_deg_lon
    y = (np.asarray(lat, dtype=float) - lat0) * KM_PER_DEG_LAT
    return np.column_stack((x, y))


def fit_pca(points_xy):
    """Sample-covariance PCA of km coordinates.

    Returns (centre, eigenvectors, eigenvalues) with eigenvalues in
    descending order (clipped at zero); eigenvector column 0 is the major
    axis. A single-cell cluster gets unit eigenvalues.
    """
    points_xy = np.asarray(points_xy, dtype=float)
    if points_xy.shape[0] == 1:
        return points_xy[0].copy(), np.eye(2), np.array([1.0, 1.0])
    center = points_xy.mean(axis=0)
    cov = np.cov(points_xy - center, rowvar=False, bias=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    return center, eigvecs[:, order], np.clip(eigvals[order], 0.0, None)


def ellipse_semi_axes_km(eigvals, scale_factor, eigenvalue_floor_km2):
    """Semi-axes a = s*sqrt(max(lambda1, floor)), b = s*sqrt(max(lambda2, floor))."""
    a = scale_factor * math.sqrt(max(float(eigvals[0]), eigenvalue_floor_km2))
    b = scale_factor * math.sqrt(max(float(eigvals[1]), eigenvalue_floor_km2))
    return a, b


def major_axis_orientation_east_ccw(eigvecs):
    """Major-axis angle counterclockwise from east, in [0, 180) degrees."""
    return math.degrees(math.atan2(float(eigvecs[1, 0]), float(eigvecs[0, 0]))) % 180.0


def azimuth_from_north(orientation_east_ccw):
    """Convert to the manuscript convention: clockwise from north, (-90, 90].

    theta_north = 90 - theta_east, wrapped so that an exact east-west axis
    is +90.
    """
    out = np.mod(90.0 - np.asarray(orientation_east_ccw, dtype=float) + 90.0, 180.0) - 90.0
    return np.where(out <= -90.0, out + 180.0, out)


def ellipse_geometry(lon, lat, scale_factor, eigenvalue_floor_km2):
    """PCA ellipse of one cluster (unweighted, centred on the cell mean).

    L1 = 2 s sqrt(max(lambda1, floor)), L2 = 2 s sqrt(max(lambda2, floor)),
    area = pi L1 L2 / 4, ratio = L2 / L1.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    if lon.size == 0 or lon.size != lat.size:
        raise ValueError("an ellipse needs at least one member cell with matching lon and lat")
    lon0, lat0 = float(lon.mean()), float(lat.mean())
    center_xy, eigvecs, eigvals = fit_pca(project_lonlat_to_km(lon, lat, lon0, lat0))
    a, b = ellipse_semi_axes_km(eigvals, scale_factor, eigenvalue_floor_km2)
    orientation = major_axis_orientation_east_ccw(eigvecs)
    return {
        "n_member_cells": int(lon.size),
        "centroid_lon_unweighted": lon0,
        "centroid_lat_unweighted": lat0,
        "eigenvalue_major_km2": float(eigvals[0]),
        "eigenvalue_minor_km2": float(eigvals[1]),
        "major_axis_km": 2.0 * a,
        "minor_axis_km": 2.0 * b,
        "area_km2": math.pi * a * b,
        "axis_ratio": b / a,
        "orientation_east_ccw_deg": orientation,
        "azimuth_north_cw_deg": float(azimuth_from_north(orientation)),
        "major_axis_east": math.cos(math.radians(orientation)),
        "major_axis_north": math.sin(math.radians(orientation)),
        "_center_xy": center_xy, "_eigvecs": eigvecs, "_eigvals": eigvals,
    }


def points_inside_ellipse(lon, lat, center_lon, center_lat, center_xy, eigvecs,
                          eigvals, scale_factor, eigenvalue_floor_km2,
                          scale_factor_minor=None):
    """True for points inside the ellipse (boundary included, 1e-12 slack).

    A different scale factor can be given for the minor axis, as used in
    the Appendix B sensitivity analysis.
    """
    if scale_factor_minor is None:
        scale_factor_minor = scale_factor
    rel = project_lonlat_to_km(lon, lat, center_lon, center_lat) - center_xy
    pc = rel @ eigvecs
    a = scale_factor * math.sqrt(max(float(eigvals[0]), eigenvalue_floor_km2))
    b = scale_factor_minor * math.sqrt(max(float(eigvals[1]), eigenvalue_floor_km2))
    return ((pc[:, 0] / a) ** 2 + (pc[:, 1] / b) ** 2) <= 1.0 + 1e-12


def ellipse_outline_lonlat(geometry, n_points=120):
    """Points on the ellipse boundary in degrees (for drawing maps).

    Uses the same local projection as the fit, so the outline is exact for
    the geometry returned by ``ellipse_geometry``.
    """
    lon0, lat0 = geometry["centroid_lon_unweighted"], geometry["centroid_lat_unweighted"]
    a, b = geometry["major_axis_km"] / 2.0, geometry["minor_axis_km"] / 2.0
    theta = math.radians(geometry["orientation_east_ccw_deg"])
    t = np.linspace(0.0, 2.0 * math.pi, n_points)
    x = a * np.cos(t) * math.cos(theta) - b * np.sin(t) * math.sin(theta)
    y = a * np.cos(t) * math.sin(theta) + b * np.sin(t) * math.cos(theta)
    km_per_deg_lon = KM_PER_DEG_LON_AT_EQUATOR * math.cos(math.radians(lat0))
    return lon0 + x / km_per_deg_lon, lat0 + y / KM_PER_DEG_LAT


def weighted_centroid(lon, lat, tmax_c):
    """Temperature-weighted centroid: sum(lon*T)/sum(T), sum(lat*T)/sum(T).

    T is the raw daily Tmax in degrees Celsius of the cluster's own member
    cells. Stops if any weight is missing, zero or negative (there is no
    silent fallback to the unweighted centroid).
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    w = np.asarray(tmax_c, dtype=float)
    if w.size == 0 or not np.all(np.isfinite(w)) or np.any(w <= 0.0):
        raise ValueError("weighted centroid needs finite, positive Tmax weights")
    wsum = float(np.sum(w))
    return float(np.sum(lon * w) / wsum), float(np.sum(lat * w) / wsum), wsum
