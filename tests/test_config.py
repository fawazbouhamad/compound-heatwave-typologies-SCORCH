# SPDX-License-Identifier: GPL-3.0-only
"""Settings: defaults, the YAML template, overrides, validation and stage selection."""
from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from scorch import config as cfg
from scorch.stages import parse_stage_selection

REPO = Path(__file__).resolve().parents[1]


def test_defaults_are_the_manuscript_settings():
    config = cfg.load_config(environment={})
    expected = copy.deepcopy(cfg.MANUSCRIPT_SETTINGS)
    expected["paths"]["base_dir"] = config["paths"]["base_dir"]
    assert config == expected
    assert cfg.number_of_parameter_pairs(config) == 63


def test_yaml_template_matches_the_code():
    template = yaml.safe_load((REPO / "config" / "manuscript_settings.yml").read_text(encoding="utf-8"))
    assert template == cfg.MANUSCRIPT_SETTINGS


def test_settings_file_overrides_values(tmp_path):
    settings = tmp_path / "study.yml"
    settings.write_text("local_heatwaves:\n  threshold_quantile: 0.9\nclustering:\n"
                        "  eps_values: [1.0, 2.0]\n  min_samples_values: [3, 4, 5]\n", encoding="utf-8")
    config = cfg.load_config(settings_file=settings, environment={})
    assert config["local_heatwaves"]["threshold_quantile"] == 0.9
    assert config["local_heatwaves"]["min_duration_days"] == 3      # unchanged default
    assert cfg.number_of_parameter_pairs(config) == 6


def test_toml_and_json_settings_files(tmp_path):
    toml = tmp_path / "s.toml"
    toml.write_text('[regional_days]\nquantile = 0.99\n', encoding="utf-8")
    assert cfg.load_config(settings_file=toml, environment={})["regional_days"]["quantile"] == 0.99
    json_file = tmp_path / "s.json"
    json_file.write_text('{"study_area": {"first_year": 1950}}', encoding="utf-8")
    assert cfg.load_config(settings_file=json_file, environment={})["study_area"]["first_year"] == 1950


def test_unknown_setting_is_an_error(tmp_path):
    settings = tmp_path / "bad.yml"
    settings.write_text("local_heatwaves:\n  threshold: 0.9\n", encoding="utf-8")
    with pytest.raises(cfg.SettingsError, match="unknown setting 'local_heatwaves.threshold'"):
        cfg.load_config(settings_file=settings, environment={})
    with pytest.raises(cfg.SettingsError, match="not found"):
        cfg.load_config(settings_file=tmp_path / "missing.yml", environment={})


@pytest.mark.parametrize("overrides, message", [
    ({"local_heatwaves": {"threshold_quantile": 1.5}}, "threshold_quantile"),
    ({"regional_days": {"quantile_method": "nearest_rank"}}, "quantile_method"),
    ({"clustering": {"eps_values": []}}, "must not be empty"),
    ({"clustering": {"min_samples_values": [0, 4]}}, "positive integers"),
    ({"study_area": {"warm_season_months": [4, 13]}}, "months between 1 and 12"),
    ({"ellipses": {"weighted_centroid_tmax_precision": "float16"}}, "float32 or float64"),
])
def test_invalid_settings_are_rejected(overrides, message):
    with pytest.raises(cfg.SettingsError, match=message):
        cfg.load_config(overrides=overrides, environment={})


def test_environment_variables_override_folders(tmp_path):
    environment = {"SCORCH_INPUT_DIR": str(tmp_path / "in"), "SCORCH_OUTPUT_DIR": "results",
                   "SCORCH_BASE_DIR": str(tmp_path)}
    config = cfg.load_config(environment=environment)
    assert cfg.input_file(config) == tmp_path / "in" / "tmax_processed.nc"
    assert cfg.output_dir(config) == (tmp_path / "results").resolve()
    assert cfg.run_records_dir(config) == (tmp_path / ".scorch" / "logs" / "run-records").resolve()


def test_stage_selection():
    assert parse_stage_selection("all") == [1, 2, 3, 4, 5]
    assert parse_stage_selection("2-4") == [2, 3, 4]
    assert parse_stage_selection("5,1") == [1, 5]
    with pytest.raises(ValueError, match="no stage 6"):
        parse_stage_selection("6")
