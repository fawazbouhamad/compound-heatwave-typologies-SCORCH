# Changelog

All notable changes to the SCORCH software. The processed temperature
dataset has its own version (1.0.0) and DOI (10.5281/zenodo.21717874); the
entries below are about the software only.

## 1.0.0 - 2026-10-11

First citable software release. The scientific results of the study are
unchanged: the five stages reproduce every published catalog, Table 1 and
the supporting results byte for byte, the power-law point estimates and the
seven LGCP parameters are identical, and the Appendix C numbers are
reproduced (see `docs/reproduction_record.md`).

### Added

- Python package `scorch` in `src/` (settings, input reading, the five stage
  calculations, stage runners) and the `scorch` command: `run`, `demo`,
  `check-input`, `settings`, `info`, with explicit input and output folders,
  settings files (YAML, TOML or JSON) and plain-language error messages.
- A synthetic demo (`scorch demo`) with six designed events covering all
  four event types, an overview figure labelled as synthetic, a summary and
  recorded expected results.
- Tests (`tests/`): threshold equality, bridged and split heatwaves, season
  boundaries, the regional threshold, DBSCAN parameter counts and ties,
  rounding, ellipse geometry and orientation conventions, weighted
  centroids, event classification, Table 1 with a missing type, trend
  statistics with few years, input errors, settings, the command line, the
  demo end to end, and a regression test against the published study
  results (`SCORCH_STUDY_INPUT`). R function checks in `tests/r/`.
- Continuous integration (GitHub Actions) on Ubuntu, Windows and macOS with
  the pinned package versions, an informative job with current versions on
  Python 3.13, and an R job.
- Installation files: `pyproject.toml`, pinned `requirements/core.txt`,
  `requirements/figures.txt` and `requirements/dev.txt`, `environment/renv.lock`
  and `environment/install_r_packages.R`; `config/manuscript_settings.yml`
  with every study setting explained.
- Documentation in `docs/`: installation, usage, input/output reference,
  demo, reproduction record, SoftwareX checklist, release checklist; a
  rewritten README; `CONTRIBUTING.md`; `LICENSE-DATA-AND-ARTWORK.md`.
- Citation and archiving metadata: `CITATION.cff`, `codemeta.json`,
  `.zenodo.json`; SPDX licence identifiers in every Python and R source file.

### Changed

- `LICENSE` renamed to `LICENSE.txt` (same GPL-3.0-only text), as the
  SoftwareX guide asks.
- `scripts/main_workflow/helpers/functions.py` now re-exports the package
  functions, and the five stage scripts call the package; every script keeps
  its name and command line.
- The LGCP preparation and Figure C1 scripts compare counts with the
  catalogs instead of the study's numbers (753 ellipses, 395 days, 51 events);
  the appendix scripts write their run records to the configured folder.
- Stage 1 checks that `grid_cells.csv` matches the NetCDF coordinates; stage
  2 stops with an explanation when the regional threshold would be zero.

### Fixed

- An empty event or an event with a day without clusters is now rejected
  with an error instead of being classified as Type 2.
- Table 1 handles a dataset with no events of one type (zero events, empty
  statistics) instead of failing; trend statistics handle types with fewer
  than three years.
- The stage 3 completeness check uses the configured numbers of `eps` and
  `min_samples` values instead of the hard-coded 63.

### Development notes

The packaging, tests, demo and documentation of this release were prepared
with the assistance of an AI coding tool (Claude Code by Anthropic) operated
and reviewed by Fawaz Bouhamad; the scientific method, settings and results
are those of the authors and were verified against the published outputs as
recorded in `docs/reproduction_record.md`.

## Before 1.0.0

The public repository started on 8 October 2026 with a single commit holding
the research record of the study (scripts, catalogs, figures, supporting
files). Earlier development took place in a private repository; earlier
public revisions of the SCORCH code were released under the MIT License and
keep that licence.
