# SPDX-License-Identifier: GPL-3.0-only
"""Stage 3 rules: parameter search size, modal-count ties, rounding, cell order."""
from __future__ import annotations

import numpy as np
import pytest

from scorch import clustering, config as cfg


def disc(center_lon, center_lat, radius):
    lon, lat = np.meshgrid(np.arange(0.5, 40.5), np.arange(0.5, 30.5))
    inside = (lon - center_lon) ** 2 + (lat - center_lat) ** 2 <= radius ** 2
    return lon[inside], lat[inside]


def test_parameter_search_tries_every_pair():
    lon, lat = disc(10.5, 10.5, 3)
    distances = clustering.grid_distance_matrix(lon, lat)
    rows = clustering.parameter_search(distances, [1.0, 2.0, 3.0], [4, 5, 6, 7])
    assert len(rows) == 12
    config = cfg.load_config(overrides={"clustering": {"eps_values": [1.0, 2.0, 3.0],
                                                       "min_samples_values": [4, 5, 6, 7]}},
                             environment={})
    assert cfg.number_of_parameter_pairs(config) == 12       # not the manuscript's 63
    # One disc: never more than one cluster; with eps >= 2 every min_samples finds it,
    # with eps = 1 (four neighbours) the larger min_samples values find none.
    assert all(r["n_clusters"] in (0, 1) for r in rows)
    assert all(r["n_clusters"] == 1 for r in rows if r["eps"] >= 2.0)


def test_two_distant_discs_give_two_clusters_in_west_to_east_order():
    lon_a, lat_a = disc(5.5, 10.5, 3)
    lon_b, lat_b = disc(30.5, 20.5, 3)
    lon = np.concatenate([lon_b, lon_a])          # deliberately east first
    lat = np.concatenate([lat_b, lat_a])
    order = clustering.longitude_major_order(lon, lat)
    lon, lat = lon[order], lat[order]
    assert lon[0] == lon.min()                    # west to east ...
    assert lat[0] == lat[lon == lon[0]].min()     # ... and south to north within a column
    labels, is_core = clustering.run_dbscan(clustering.grid_distance_matrix(lon, lat), 1.5, 4)
    assert clustering.count_clusters(labels) == 2
    assert labels[0] == 0                         # the western disc is found first
    assert is_core.sum() > 0 and (labels == clustering.NOISE).sum() == 0


def test_round_half_up_differs_from_bankers_rounding():
    assert clustering.round_half_up(6.5) == 7 and round(6.5) == 6
    assert clustering.round_half_up(7.5) == 8 and clustering.round_half_up(7.49) == 7
    assert clustering.round_half_up(-0.5) == 0


def test_select_daily_parameters_keeps_tied_modal_counts():
    rows = [
        {"eps": 1.0, "min_samples": 4, "n_clusters": 1, "n_noise_cells": 0},
        {"eps": 2.0, "min_samples": 5, "n_clusters": 1, "n_noise_cells": 0},
        {"eps": 3.0, "min_samples": 8, "n_clusters": 2, "n_noise_cells": 0},
        {"eps": 4.0, "min_samples": 9, "n_clusters": 2, "n_noise_cells": 0},
        {"eps": 4.0, "min_samples": 12, "n_clusters": 0, "n_noise_cells": 50},   # zeros never count
    ]
    choice = clustering.select_daily_parameters(rows)
    assert choice["modal_cluster_counts"] == [1, 2] and choice["modal_frequency"] == 2
    assert choice["n_pairs_retained"] == 4
    assert choice["daily_eps"] == pytest.approx(2.5)
    assert choice["daily_min_samples_mean"] == pytest.approx(6.5)
    assert choice["daily_min_samples"] == 7                    # half rounded up
    assert choice["retained_pairs"] == {(1.0, 4), (2.0, 5), (3.0, 8), (4.0, 9)}


def test_select_daily_parameters_when_every_pair_gives_zero_clusters():
    rows = [{"eps": e, "min_samples": m, "n_clusters": 0, "n_noise_cells": 3}
            for e in (1.0, 2.0) for m in (4, 5)]
    choice = clustering.select_daily_parameters(rows)
    assert choice["modal_cluster_counts"] == [0] and choice["n_pairs_retained"] == 4
    assert choice["daily_min_samples"] == 5                    # mean 4.5 rounded half up
    with pytest.raises(ValueError):
        clustering.select_daily_parameters([])


def test_box_statistics_handle_empty_groups():
    stats = clustering.box_statistics(np.array([]))
    assert stats["n"] == 0 and np.isnan(stats["median"])
    stats = clustering.box_statistics(np.array([1.0, 2.0, 3.0, 4.0]))
    assert stats == {"n": 4, "mean": 2.5, "median": 2.5, "min": 1.0, "max": 4.0, "q1": 1.75, "q3": 3.25}
