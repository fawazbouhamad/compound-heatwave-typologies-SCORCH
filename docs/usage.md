# Running SCORCH

## The `scorch` command

`pip install -e .` installs the `scorch` command. Run it from the repository
folder.

| Command | What it does |
|---|---|
| `scorch run` | runs the five stages in order on `input/tmax_processed.nc` and writes the catalogs to `output/` (the manuscript reproduction) |
| `scorch run --stages 3-5` | runs only some stages; each stage reads what the earlier ones wrote |
| `scorch run --input-dir D --output-dir E` | reads the temperature file from folder `D` and writes catalogs, tables, supporting results and run records under `E` |
| `scorch run --settings my.yml ...` | applies your own study settings (see below) |
| `scorch demo` | builds the synthetic example and runs it end to end ([demo.md](demo.md)) |
| `scorch check-input FILE [--strict]` | describes a NetCDF file and checks that SCORCH can read it (`--strict` also checks the study calendar and grid) |
| `scorch settings [--write FILE]` | prints the settings in force, or writes them to a file to edit |
| `scorch info` | versions of SCORCH and its packages |

Every command accepts `--debug` before the subcommand to show a full
traceback when something fails. Errors are reported in plain language, for
example a missing input file names the dataset record to download, and an
unknown setting lists the settings of its section.

## Folders

| Option | Default | Content |
|---|---|---|
| `--input-dir` | `input` | the processed temperature file (`--input-file`, default `tmax_processed.nc`) |
| `--supporting-input-dir` | `supporting_files/additional_inputs` | `maps/grid_cells.csv` and, for the study region, the country boundaries |
| `--output-dir` | `output` | `catalogs/`, `tables/` and `figures/` |
| `--supporting-results-dir` | `supporting_files/additional_results` | intermediate results of every stage (thresholds, daily summary, DBSCAN searches, model inputs and outputs) |
| `--run-records-dir` | `.scorch/logs/run-records` | one JSON file per stage with the checksums of what it read and wrote and the package versions |

When `--output-dir` is given without `--supporting-results-dir` or
`--run-records-dir`, the supporting results and run records go under the
output folder (`E/supporting_results`, `E/run_records`), so one folder holds
a complete run.

## Your own study: settings without editing code

All study settings live in [`config/manuscript_settings.yml`](../config/manuscript_settings.yml)
with explanations. SCORCH always starts from the manuscript values; a
settings file only needs the keys you change:

```yaml
# my_study.yml
study_area:
  first_year: 1991
  last_year: 2020
  warm_season_months: [5, 6, 7, 8, 9]
  expected_dates: null        # derived from the years and months
  expected_cells: null        # taken from the input file
local_heatwaves:
  threshold_quantile: 0.9
summaries:
  in_text_statistics: false   # the Section 4 statistics describe the Middle East only
```

```bash
scorch run --settings my_study.yml --input-dir my_data --supporting-input-dir my_data --output-dir my_results
```

Your input folder needs the NetCDF file and `maps/grid_cells.csv`
(see [io_reference.md](io_reference.md); `scorch.io.grid_cells_table`
writes the table for any regular grid). The stage scripts read the same
settings from the environment variables `SCORCH_CONFIG`, `SCORCH_INPUT_DIR`,
`SCORCH_OUTPUT_DIR`, `SCORCH_SUPPORTING_INPUT_DIR`,
`SCORCH_SUPPORTING_RESULTS_DIR`, `SCORCH_RUN_RECORDS_DIR` and `SCORCH_BASE_DIR`.

## The five stages in words

1. **Local heatwaves.** For every grid cell, the warm-season 95th percentile
   of daily maximum temperature (Tmax) over all study days is the cell's
   threshold. A day exceeds when Tmax is at or above it. A local heatwave is
   a run of at least 3 days with at least 3 exceeding days; a single
   non-exceeding day inside the run is bridged, two consecutive
   non-exceeding days end it. 30 September and the following 1 April are
   not consecutive days.
2. **Regionally extensive days.** The daily number of cells in a heatwave is
   counted over all study days (including days with none). Days with at
   least the 97.5th percentile of that count (the next observed value at or
   above the quantile position) are selected.
3. **Clusters.** On every selected day the heatwave cells are clustered with
   DBSCAN (distance in grid cells) for every pair of `eps` (1.0 to 4.0) and
   `min_samples` (4 to 12). The most frequent non-zero number of clusters
   decides which pairs are kept; ties keep the pairs of all tied counts. The
   day's settings are the mean `eps` and the mean `min_samples` (rounded
   half up) of the kept pairs. Within a run of consecutive selected days the
   largest daily `eps` and `min_samples` are applied to every day.
4. **Ellipses.** Each cluster's cells are projected to kilometres around
   their mean position; the sample covariance gives a PCA ellipse with
   semi-axes `1.25 * sqrt(eigenvalue)` (eigenvalues below 1 km² are raised
   to 1). The temperature-weighted centroid is the Tmax-weighted mean
   position of the member cells.
5. **Compound events.** Consecutive selected days form one compound event.
   The number of clusters on each day gives the type (Table 1):

| Type | Name | Rule |
|---|---|---|
| 1 | Widespread (Isolated) | one day, one cluster |
| 2 | Spatially Clustered | one day, two or more clusters |
| 3 | Temporally Clustered | two or more days, all with one cluster or all with several |
| 4 | Compound Clustering (Multi-Type) | two or more days mixing one-cluster and several-cluster days |

An event without a day, or with a day without clusters, cannot occur in a
consistent catalog; SCORCH stops with an error instead of classifying it.

## The manuscript scripts

The original stage scripts still work and give exactly the same files:

```bash
python scripts/main_workflow/01_detect_heatwaves.py
python scripts/main_workflow/02_select_regional_days.py
python scripts/main_workflow/03_cluster_heatwaves.py
python scripts/main_workflow/04_fit_ellipses.py
python scripts/main_workflow/05_summarize_compound_events.py
```

Each is a short wrapper around the matching module in `src/scorch/stages/`.
The figure scripts, the R model scripts and the appendix scripts are run as
described in the [reproduction guide](../supporting_files/instructions/reproduction_guide.md).

## Reading the code

| Module | Content |
|---|---|
| `src/scorch/config.py` | the manuscript settings, settings files, folder helpers |
| `src/scorch/io.py` | reading the NetCDF input and the grid table, run records, checksums |
| `src/scorch/heatwaves.py` | stage 1 calculations |
| `src/scorch/regional.py` | stage 2 calculations |
| `src/scorch/clustering.py` | stage 3 calculations and the DBSCAN summaries of Appendix A |
| `src/scorch/ellipses.py` | stage 4 calculations |
| `src/scorch/events.py` | stage 5: event types, Table 1, axial statistics, trend tests |
| `src/scorch/geography.py` | the country of a point (Section 4 statistics) |
| `src/scorch/stages/` | the five stage runners (what each stage reads and writes) |
| `src/scorch/demo/` | the synthetic example and its figure |
| `src/scorch/cli.py` | the `scorch` command |
| `scripts/` | the manuscript scripts: stage wrappers, figures, R models, appendices |
