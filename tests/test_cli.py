# SPDX-License-Identifier: GPL-3.0-only
"""The scorch command: help, info, settings, check-input and error messages."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scorch import __version__
from scorch.cli import main
from conftest import write_small_netcdf


def test_version_and_info(capsys):
    with pytest.raises(SystemExit) as stop:
        main(["--version"])
    assert stop.value.code == 0
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert __version__ in out and "pandas" in out


def test_settings_are_printed_and_written(tmp_path, capsys):
    assert main(["settings"]) == 0
    assert "threshold_quantile: 0.95" in capsys.readouterr().out
    target = tmp_path / "s.yml"
    assert main(["settings", "--write", str(target)]) == 0
    assert "eps_values" in target.read_text(encoding="utf-8")


def test_missing_input_gives_a_helpful_error(tmp_path, capsys):
    code = main(["run", "--input-dir", str(tmp_path / "nowhere"), "--output-dir", str(tmp_path / "out")])
    assert code == 1
    err = capsys.readouterr().err
    assert "input file not found" in err and "zenodo" in err and "scorch demo" in err


def test_check_input(tmp_path, capsys):
    dates = pd.date_range("2000-04-01", "2000-09-30", freq="D")
    tmax = np.random.default_rng(2).normal(30, 2, size=(len(dates), 2, 3))
    path = write_small_netcdf(tmp_path / "t.nc", dates, [10.5, 11.5], [20.5, 21.5, 22.5], tmax)
    assert main(["check-input", str(path)]) == 0
    out = capsys.readouterr().out
    assert "SCORCH can read this file" in out and "183" in out
    # Against the manuscript calendar the small file is rejected with a clear message.
    assert main(["check-input", str(path), "--strict"]) == 1
    assert "warm-season calendar" in capsys.readouterr().err


def test_bad_settings_file_is_reported(tmp_path, capsys):
    settings = tmp_path / "bad.yml"
    settings.write_text("clustering:\n  epsilon: [1.0]\n", encoding="utf-8")
    assert main(["run", "--settings", str(settings)]) == 1
    assert "unknown setting 'clustering.epsilon'" in capsys.readouterr().err
