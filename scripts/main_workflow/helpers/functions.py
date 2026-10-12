# SPDX-License-Identifier: GPL-3.0-only
"""Calculations shared by the analysis, models and supplementary scripts.

The calculations live in the ``scorch`` package under ``src/scorch``:

    scorch.config      study settings and folder locations
    scorch.io          reading the NetCDF input, grid table and run records
    scorch.heatwaves   stage 1: thresholds and local heatwaves
    scorch.regional    stage 2: regional days and consecutive-day sequences
    scorch.clustering  stage 3: DBSCAN and its parameter selection
    scorch.ellipses    stage 4: PCA ellipses and weighted centroids
    scorch.events      stage 5: event types, Table 1, axial and trend statistics
    scorch.geography   country of a point (Section 4 statistics)

This module keeps the names the scripts have always imported
(``import functions as f``). It works whether or not the package has been
installed with pip: when it has not, the ``src`` folder of the repository is
added to the module search path.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parents[3]

try:
    import scorch  # noqa: F401
except ImportError:  # running from a checkout without ``pip install -e .``
    sys.path.insert(0, str(REPO_DIR / "src"))
    import scorch  # noqa: F401

from scorch.config import (  # noqa: E402,F401
    MANUSCRIPT_SETTINGS, SettingsError, base_dir, catalog_dir, country_boundaries_path,
    input_dir, input_file, load_config, number_of_parameter_pairs, output_dir,
    output_subdir, project_path, relative_to_repo, run_records_dir,
    supporting_input_dir, supporting_results_dir,
)
from scorch.io import (  # noqa: E402,F401
    InputError, centroid_temperature_weights, check_grid_matches, check_study_calendar,
    expected_study_dates, grid_cells_table, now_utc, read_grid_cells, read_netcdf_tmax,
    read_tmax_input, sha256_of, software_versions, write_csv, write_json, write_run_record,
)
from scorch.heatwaves import (  # noqa: E402,F401
    cell_thresholds, detect_all_cells, label_cell_by_season, label_heatwaves,
)
from scorch.regional import (  # noqa: E402,F401
    attach_compound_events, group_consecutive_days, regional_threshold,
)
from scorch.clustering import (  # noqa: E402,F401
    DBSCAN_GROUPS, NOISE, box_statistics, count_clusters, figure_a2_panels,
    grid_distance_matrix, longitude_major_order, parameter_search, round_half_up,
    run_dbscan, select_daily_parameters, statistics_table, values_by_group,
)
from scorch.ellipses import (  # noqa: E402,F401
    KM_PER_DEG_LAT, KM_PER_DEG_LON_AT_EQUATOR, azimuth_from_north, ellipse_geometry,
    ellipse_outline_lonlat, ellipse_semi_axes_km, fit_pca, major_axis_orientation_east_ccw,
    points_inside_ellipse, project_lonlat_to_km, weighted_centroid,
)
from scorch.events import (  # noqa: E402,F401
    TYPE_DESCRIPTIONS, TYPE_NAMES, axial_mean_and_resultant, classify_event,
    mann_kendall_test, sens_slope, table_1, trend_statistics,
)
from scorch.geography import country_of_points  # noqa: E402,F401
