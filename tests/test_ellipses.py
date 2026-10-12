# SPDX-License-Identifier: GPL-3.0-only
"""Stage 4 rules: ellipse geometry, orientation conventions and weighted centroids."""
from __future__ import annotations

import math

import numpy as np
import pytest

from scorch import ellipses as ell

S, FLOOR = 1.25, 1.0


def test_east_west_line_of_cells():
    lon = np.array([20.5, 21.5, 22.5, 23.5, 24.5])
    lat = np.full(5, 30.5)
    g = ell.ellipse_geometry(lon, lat, S, FLOOR)
    km_per_deg = ell.KM_PER_DEG_LON_AT_EQUATOR * math.cos(math.radians(30.5))
    major_eigenvalue = 2.5 * km_per_deg ** 2                 # sample variance of -2..2 degrees
    assert g["n_member_cells"] == 5
    assert g["eigenvalue_major_km2"] == pytest.approx(major_eigenvalue)
    assert g["eigenvalue_minor_km2"] == pytest.approx(0.0)
    assert g["major_axis_km"] == pytest.approx(2 * S * math.sqrt(major_eigenvalue))
    assert g["minor_axis_km"] == pytest.approx(2 * S * math.sqrt(FLOOR))   # floor applies
    assert g["area_km2"] == pytest.approx(math.pi * g["major_axis_km"] * g["minor_axis_km"] / 4)
    assert g["axis_ratio"] == pytest.approx(g["minor_axis_km"] / g["major_axis_km"])
    assert g["orientation_east_ccw_deg"] == pytest.approx(0.0)
    assert g["azimuth_north_cw_deg"] == pytest.approx(90.0)     # an east-west axis is +90


def test_single_cell_cluster_has_unit_eigenvalues():
    g = ell.ellipse_geometry([45.5], [28.5], S, FLOOR)
    assert g["major_axis_km"] == g["minor_axis_km"] == pytest.approx(2 * S)
    assert g["area_km2"] == pytest.approx(math.pi * S ** 2) and g["axis_ratio"] == 1.0


def test_square_block_floor_not_applied():
    lon, lat = np.meshgrid(np.arange(30.5, 35.5), np.arange(30.5, 35.5))
    g = ell.ellipse_geometry(lon.ravel(), lat.ravel(), S, FLOOR)
    assert g["eigenvalue_minor_km2"] > FLOOR
    assert g["axis_ratio"] == pytest.approx(
        math.sqrt(g["eigenvalue_minor_km2"] / g["eigenvalue_major_km2"]))


@pytest.mark.parametrize("east_ccw, north_cw", [(0.0, 90.0), (90.0, 0.0), (45.0, 45.0),
                                                (135.0, -45.0), (179.0, -89.0)])
def test_azimuth_convention(east_ccw, north_cw):
    assert float(ell.azimuth_from_north(east_ccw)) == pytest.approx(north_cw)


def test_points_inside_ellipse_includes_the_boundary():
    lon = np.array([20.5, 21.5, 22.5, 23.5, 24.5])
    lat = np.full(5, 30.5)
    g = ell.ellipse_geometry(lon, lat, S, FLOOR)
    outline_lon, outline_lat = ell.ellipse_outline_lonlat(g, n_points=16)
    inside = ell.points_inside_ellipse(outline_lon, outline_lat, g["centroid_lon_unweighted"],
                                       g["centroid_lat_unweighted"], g["_center_xy"], g["_eigvecs"],
                                       g["_eigvals"], S, FLOOR)
    assert inside.all()
    far = ell.points_inside_ellipse([30.0], [30.5], g["centroid_lon_unweighted"],
                                    g["centroid_lat_unweighted"], g["_center_xy"], g["_eigvecs"],
                                    g["_eigvals"], S, FLOOR)
    assert not far.any()


def test_weighted_centroid_uses_temperature_weights():
    lon_w, lat_w, total = ell.weighted_centroid([0.0, 1.0], [0.0, 2.0], [1.0, 3.0])
    assert (lon_w, lat_w, total) == (0.75, 1.5, 4.0)
    for bad in ([0.0, 3.0], [-1.0, 3.0], [np.nan, 3.0], []):
        with pytest.raises(ValueError, match="finite, positive"):
            ell.weighted_centroid([0.0, 1.0][:len(bad)], [0.0, 2.0][:len(bad)], bad)


def test_fit_pca_orders_eigenvalues_descending():
    points = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 1.0], [10.0, 1.0]])
    _, vectors, values = ell.fit_pca(points)
    assert values[0] >= values[1] and abs(vectors[0, 0]) == pytest.approx(1.0)


def test_empty_cluster_is_rejected():
    with pytest.raises(ValueError, match="at least one member cell"):
        ell.ellipse_geometry([], [], S, FLOOR)
