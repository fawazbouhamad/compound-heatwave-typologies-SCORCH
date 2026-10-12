# SPDX-License-Identifier: GPL-3.0-only
"""Study settings and folder locations.

``MANUSCRIPT_SETTINGS`` records every setting of the published analysis. The
values are also written out, with explanations, in
``config/manuscript_settings.yml``; a test keeps the two copies identical.

``load_config()`` returns a fresh copy of these settings. A different study
can be run without editing any Python file: give ``load_config`` a settings
file (YAML, TOML or JSON) that overrides some of the values, and point the
folder settings at your own input and output folders. The ``scorch`` command
does this through its options; the stage scripts read the same information
from environment variables (``SCORCH_CONFIG``, ``SCORCH_INPUT_DIR``,
``SCORCH_OUTPUT_DIR``, ``SCORCH_SUPPORTING_INPUT_DIR``,
``SCORCH_SUPPORTING_RESULTS_DIR``, ``SCORCH_RUN_RECORDS_DIR`` and
``SCORCH_BASE_DIR``).
"""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path

# The settings of the published study. Values marked informational describe
# the study and are not read by the calculations.
MANUSCRIPT_SETTINGS = {
    # Folders. Relative paths are resolved against ``base_dir`` (the
    # repository folder when SCORCH runs from a checkout, otherwise the
    # current working directory). ``run_records_dir`` defaults to
    # ``.scorch/logs/run-records`` under ``base_dir``.
    "paths": {
        "base_dir": None,
        "input_dir": "input",
        "input_file": "tmax_processed.nc",
        "supporting_input_dir": "supporting_files/additional_inputs",
        "supporting_results_dir": "supporting_files/additional_results",
        "output_dir": "output",
        "run_records_dir": None,
    },
    # Bounds are informational. The input supplies the grid; counts are checked.
    # ``expected_dates`` and ``expected_cells`` may be null, in which case the
    # calendar is derived from the years and months and the cell count is
    # taken from the input file.
    "study_area": {
        "lat_min": 10.0,
        "lat_max": 46.0,
        "lon_min": 20.0,
        "lon_max": 70.0,
        "grid_step_deg": 1.0,
        "warm_season_months": [4, 5, 6, 7, 8, 9],
        "first_year": 1940,
        "last_year": 2025,
        "expected_dates": 15738,
        "expected_cells": 1800,
    },
    # Linear percentile; Tmax >= threshold. Isolated gaps are bridged.
    "local_heatwaves": {
        "threshold_quantile": 0.95,
        "min_duration_days": 3,
        "min_exceeding_days": 3,
    },
    # Includes zero-count days; the higher method selects an observed value.
    "regional_days": {
        "quantile": 0.975,
        "quantile_method": "higher",
    },
    # Neighbourhood distances are in grid-cell units.
    "clustering": {
        "eps_values": [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
        "min_samples_values": [4, 5, 6, 7, 8, 9, 10, 11, 12],
    },
    # Only centroid weights use float32; other steps retain input precision.
    "ellipses": {
        "scale_factor": 1.25,
        "eigenvalue_floor_km2": 1.0,
        "weighted_centroid_tmax_precision": "float32",
    },
    # Stage 5 summaries. The in-text statistics of Section 4 need the Natural
    # Earth country boundaries and describe the study region only; switch
    # them off for other regions. Trend statistics are written when a type
    # has at least three years with events, otherwise the row is left empty.
    "summaries": {
        "in_text_statistics": True,
        "country_boundaries": "maps/countries_10m/ne_10m_admin_0_countries.shp",
    },
    "representative_events": {
        "figure_5_type_1": ["2002-07-31"],
        "figure_5_type_2": ["1983-07-13"],
        "figure_6a_type_3": ["2017-07-17", "2017-07-18"],
        "figure_6b_type_3": ["1991-06-04", "1991-06-05"],
        "figure_7_type_4": [
            "2001-08-07", "2001-08-08", "2001-08-09",
            "2001-08-10", "2001-08-11", "2001-08-12",
        ],
        "figure_s1_type_3": ["2016-06-07", "2016-06-08"],
    },
    # Clauset fit: all unique cutoffs except the largest; smaller cutoff on ties.
    # Uses the MLE exponent and full two-sided KS distance.
    # The significance level is informational.
    "power_law": {
        "ks_mode": "full",
        "n_bootstrap": 5000,
        "seed_daily_maxima": 20261020,
        "seed_event_maxima": 20261120,
        "significance_level": 0.1,
    },
    # Projection is informational; the LGCP preparation derives it from the data.
    # Lambert azimuthal equal-area, WGS84, centred at 45 E / 28 N.
    "lgcp": {
        "projection": (
            "+proj=laea +lat_0=28 +lon_0=45 +x_0=0 +y_0=0 "
            "+ellps=WGS84 +units=m +no_defs"
        ),
        "raster_resolution_km": 25,
        "trend_terms": ["lon", "lat", "mean_tmax", "std_tmax"],
        "seed": 20260704,
    },
    # Scale factors are tested on the two axes independently.
    "appendix_b": {
        "scale_factors": [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5],
    },
    "appendix_c": {
        "n_folds": 5,
        "fold_seed": 20260704,
        "zone_area_fractions": [0.1, 0.2, 0.3, 0.5],
    },
    "supplementary_power_law": {
        "n_bootstrap": 5000,
        "base_seed": 20260918,
    },
}

# Environment variables that override folder settings (used by the scripts
# and set by the ``scorch`` command).
PATH_ENVIRONMENT_VARIABLES = {
    "base_dir": "SCORCH_BASE_DIR",
    "input_dir": "SCORCH_INPUT_DIR",
    "input_file": "SCORCH_INPUT_FILE",
    "supporting_input_dir": "SCORCH_SUPPORTING_INPUT_DIR",
    "supporting_results_dir": "SCORCH_SUPPORTING_RESULTS_DIR",
    "output_dir": "SCORCH_OUTPUT_DIR",
    "run_records_dir": "SCORCH_RUN_RECORDS_DIR",
}
SETTINGS_FILE_VARIABLE = "SCORCH_CONFIG"


class SettingsError(ValueError):
    """A settings file or override is not usable."""


def repository_root():
    """The repository folder when this package runs from a checkout, else None."""
    candidate = Path(__file__).resolve().parents[2]
    if (candidate / "scripts" / "main_workflow").is_dir() and (candidate / "pyproject.toml").is_file():
        return candidate
    return None


def default_base_dir(environment=None):
    """Folder against which relative study paths are resolved."""
    environment = os.environ if environment is None else environment
    explicit = environment.get(PATH_ENVIRONMENT_VARIABLES["base_dir"], "").strip()
    if explicit:
        return Path(explicit).expanduser().resolve()
    return repository_root() or Path.cwd().resolve()


def read_settings_file(path):
    """Read a YAML, TOML or JSON settings file into a dict."""
    path = Path(path)
    if not path.is_file():
        raise SettingsError(f"settings file not found: {path}")
    suffix = path.suffix.lower()
    text = path.read_text(encoding="utf-8")
    if suffix in (".yml", ".yaml"):
        try:
            import yaml
        except ImportError as error:  # pragma: no cover - depends on the environment
            raise SettingsError("reading YAML settings needs the PyYAML package "
                                "(pip install pyyaml), or use a .toml or .json file") from error
        data = yaml.safe_load(text)
    elif suffix == ".toml":
        import tomllib
        data = tomllib.loads(text)
    elif suffix == ".json":
        data = json.loads(text)
    else:
        raise SettingsError(f"settings files must end in .yml, .yaml, .toml or .json: {path}")
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise SettingsError(f"the settings file must contain a mapping of sections: {path}")
    return data


def merge_settings(base, overrides, _prefix=""):
    """Update ``base`` in place with ``overrides``; unknown keys are errors."""
    for key, value in overrides.items():
        where = f"{_prefix}{key}"
        if key not in base:
            raise SettingsError(f"unknown setting '{where}'. Known settings in this section: "
                                + ", ".join(sorted(base)))
        if isinstance(base[key], dict):
            if not isinstance(value, dict):
                raise SettingsError(f"setting '{where}' must be a section (mapping)")
            merge_settings(base[key], value, _prefix=where + ".")
        else:
            base[key] = copy.deepcopy(value)
    return base


def load_config(settings_file=None, overrides=None, environment=None):
    """Return the study settings as a fresh dict.

    Parameters
    ----------
    settings_file : path, optional
        A YAML, TOML or JSON file whose values replace the manuscript defaults.
        When omitted, the ``SCORCH_CONFIG`` environment variable is used if set.
    overrides : dict, optional
        Values applied after the settings file (the same nested structure).
    environment : mapping, optional
        Environment variables to read (defaults to ``os.environ``).
    """
    environment = os.environ if environment is None else environment
    config = copy.deepcopy(MANUSCRIPT_SETTINGS)
    if settings_file is None:
        settings_file = environment.get(SETTINGS_FILE_VARIABLE, "").strip() or None
    if settings_file:
        merge_settings(config, read_settings_file(settings_file))
    if overrides:
        merge_settings(config, overrides)
    for key, variable in PATH_ENVIRONMENT_VARIABLES.items():
        value = environment.get(variable, "").strip()
        if value:
            config["paths"][key] = value
    if not config["paths"].get("base_dir"):
        config["paths"]["base_dir"] = str(default_base_dir(environment))
    validate_settings(config)
    return config


def validate_settings(config):
    """Stop early when a setting cannot work."""
    area = config["study_area"]
    months = area["warm_season_months"]
    if not months or any(int(m) not in range(1, 13) for m in months):
        raise SettingsError("study_area.warm_season_months must list months between 1 and 12")
    if int(area["first_year"]) > int(area["last_year"]):
        raise SettingsError("study_area.first_year must not be after last_year")
    local = config["local_heatwaves"]
    if not 0.0 < float(local["threshold_quantile"]) < 1.0:
        raise SettingsError("local_heatwaves.threshold_quantile must lie strictly between 0 and 1")
    if int(local["min_duration_days"]) < 1 or int(local["min_exceeding_days"]) < 1:
        raise SettingsError("local_heatwaves durations must be at least 1 day")
    regional = config["regional_days"]
    if not 0.0 < float(regional["quantile"]) < 1.0:
        raise SettingsError("regional_days.quantile must lie strictly between 0 and 1")
    if regional["quantile_method"] not in ("linear", "lower", "higher", "midpoint", "nearest"):
        raise SettingsError("regional_days.quantile_method must be a pandas interpolation method")
    clustering = config["clustering"]
    if not clustering["eps_values"] or not clustering["min_samples_values"]:
        raise SettingsError("clustering.eps_values and min_samples_values must not be empty")
    if any(float(e) <= 0 for e in clustering["eps_values"]):
        raise SettingsError("clustering.eps_values must be positive distances in grid-cell units")
    if any(int(m) < 1 for m in clustering["min_samples_values"]):
        raise SettingsError("clustering.min_samples_values must be positive integers")
    ellipses = config["ellipses"]
    if float(ellipses["scale_factor"]) <= 0 or float(ellipses["eigenvalue_floor_km2"]) < 0:
        raise SettingsError("ellipses.scale_factor must be positive and eigenvalue_floor_km2 non-negative")
    if ellipses["weighted_centroid_tmax_precision"] not in ("float32", "float64"):
        raise SettingsError("ellipses.weighted_centroid_tmax_precision must be float32 or float64")
    return config


def number_of_parameter_pairs(config):
    """Number of (eps, min_samples) pairs tried on every selected day."""
    clustering = config["clustering"]
    return len(clustering["eps_values"]) * len(clustering["min_samples_values"])


# Folder helpers
def base_dir(config):
    return Path(config["paths"]["base_dir"]).expanduser().resolve()


def project_path(relative_or_absolute, config=None):
    """Turn a configured path into an absolute path."""
    path = Path(relative_or_absolute).expanduser()
    if path.is_absolute():
        return path
    root = base_dir(config) if config is not None else default_base_dir()
    return (root / path).resolve()


def input_dir(config):
    return project_path(config["paths"]["input_dir"], config)


def input_file(config):
    """The processed temperature file (``paths.input_file`` inside ``paths.input_dir``)."""
    name = Path(config["paths"]["input_file"])
    return name if name.is_absolute() else input_dir(config) / name


def supporting_input_dir(config):
    return project_path(config["paths"]["supporting_input_dir"], config)


def supporting_results_dir(config):
    return project_path(config["paths"]["supporting_results_dir"], config)


def output_dir(config):
    return project_path(config["paths"]["output_dir"], config)


def catalog_dir(config, stage_folder):
    """Folder of one catalog stage, created if needed."""
    folder = output_dir(config) / "catalogs" / stage_folder
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def output_subdir(config, name):
    folder = output_dir(config) / name
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def run_records_dir(config):
    configured = config["paths"].get("run_records_dir")
    if configured:
        return project_path(configured, config)
    return base_dir(config) / ".scorch" / "logs" / "run-records"


def country_boundaries_path(config):
    """The Natural Earth country shapefile used by the in-text statistics."""
    path = Path(config["summaries"]["country_boundaries"])
    return path if path.is_absolute() else supporting_input_dir(config) / path


def relative_to_repo(path, config=None):
    """A short path for printing: relative to the base folder when inside it."""
    path = Path(path).resolve()
    root = base_dir(config) if config is not None else default_base_dir()
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def settings_as_yaml_text(config):
    """The effective settings as YAML text (for ``scorch settings``)."""
    import yaml
    return yaml.safe_dump(config, sort_keys=False, default_flow_style=None, allow_unicode=True)
