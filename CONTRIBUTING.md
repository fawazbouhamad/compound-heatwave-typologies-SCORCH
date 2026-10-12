# Contributing

Questions, bug reports and suggestions are welcome in the
[issue tracker](https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH/issues).

## Setting up for development

```bash
git clone https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH.git
cd compound-heatwave-typologies-SCORCH
python -m venv .venv
# Windows: .venv\Scripts\activate      Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements/core.txt -r requirements/figures.txt -r requirements/dev.txt
python -m pip install -e .
python -m pytest
```

`python -m pytest` runs the unit tests and the synthetic demo end to end in
about a minute. The regression test against the published study results runs
only when `SCORCH_STUDY_INPUT` points to the archived `tmax_processed.nc`
(see [docs/reproduction_record.md](docs/reproduction_record.md)).

## What to keep stable

The settings in `src/scorch/config.py` (mirrored in
`config/manuscript_settings.yml`), the order of operations, the tie rules and
the random seeds define the published results. A change to any of them is a
new analysis, not a bug fix: keep the manuscript defaults, and add new
behaviour as a setting with a different value. The test
`tests/test_study_regression.py` guards the published catalogs; run it before
proposing a change to the calculations.

## Style

Plain, readable Python and R with short docstrings that say what a function
computes and which rule it follows. Line length up to 100 characters. New
source files start with `# SPDX-License-Identifier: GPL-3.0-only`.
