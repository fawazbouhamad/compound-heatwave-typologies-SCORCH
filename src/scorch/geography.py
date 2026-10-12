# SPDX-License-Identifier: GPL-3.0-only
"""Geographic classification of points (country of each ellipse centroid)."""
from __future__ import annotations

import numpy as np


def country_of_points(lon, lat, boundaries_path):
    """Name of the country (Natural Earth ADMIN field) containing each point.

    Uses the even-odd ray-casting test on every polygon ring, the rule used
    for the country counts reported in Section 4 (points outside every
    country, for example over the sea, get "Unassigned").
    """
    import shapefile
    points = np.column_stack((np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)))
    names = np.array(["Unassigned"] * len(points), dtype=object)
    reader = shapefile.Reader(str(boundaries_path))
    admin_field = [field[0] for field in reader.fields[1:]].index("ADMIN")
    for record in reader.iterShapeRecords():
        vertices = np.array(record.shape.points)
        parts = list(record.shape.parts) + [len(vertices)]
        inside = np.zeros(len(points), dtype=bool)
        for a, b in zip(parts[:-1], parts[1:]):
            inside ^= _points_in_ring(vertices[a:b], points)
        names[inside] = record.record[admin_field]
    return names


def _points_in_ring(ring, points):
    out = np.zeros(len(points), dtype=bool)
    if len(ring) < 3:
        return out
    xmin, ymin = ring.min(axis=0)
    xmax, ymax = ring.max(axis=0)
    subset = np.flatnonzero((points[:, 0] >= xmin) & (points[:, 0] <= xmax)
                            & (points[:, 1] >= ymin) & (points[:, 1] <= ymax))
    if len(subset) == 0:
        return out
    x, y = points[subset].T
    flags = np.zeros(len(subset), dtype=bool)
    xj, yj = ring[-1]
    for xi, yi in ring:
        if yi != yj:
            flags ^= ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        xj, yj = xi, yi
    out[subset] = flags
    return out
