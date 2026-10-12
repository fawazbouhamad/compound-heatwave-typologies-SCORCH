# Reproduction record

A verification of the published study results with this software, done on
11 October 2026 before the release commit of version 1.0.0. The supplied
catalogs, tables, model results and figures in the repository were **not
changed**; every regenerated file was compared with the committed one.

## Environment

| Item | Value |
|---|---|
| Computer | Windows 11 Home 10.0.26200 (the study computer) |
| Python | CPython 3.12.10 in a fresh virtual environment (`py -3.12 -m venv .venv`) built from `requirements/core.txt` and `requirements/figures.txt` |
| Python packages | numpy 2.4.6, pandas 3.0.3, scipy 1.17.1, scikit-learn 1.9.0, netCDF4 1.7.4, cftime 1.6.5, pyshp 3.0.12, PyYAML 6.0.3, matplotlib 3.10.9, Cartopy 0.25.0, shapely 2.1.2, pillow 12.2.0, pyproj 3.7.2 |
| R | R 4.5.1 with spatstat.geom 3.8-1, spatstat.model 3.7-1, spatstat.explore 3.8-1, jsonlite 2.0.0, poweRlaw 1.0.0, pracma 2.4.6 (as in `environment/renv.lock`) |
| Repository | the working tree of the release commit tagged `v1.0.0`, whose parent is the initial commit `fdfd5619d1587b09a1459c7ddb8dceece60a1e9e` (8 October 2026) |
| Settings | the manuscript settings (`config/manuscript_settings.yml`), no overrides |

## Input checksums

| File | Checksum |
|---|---|
| `SCORCH_data.zip` downloaded from https://zenodo.org/records/21717874 (84,435,803 bytes) | MD5 `80772df9410306cbc419191d0dcb824c` (matches the Zenodo record) |
| `scorch_data/tmax_processed.nc` extracted from it (84,709,777 bytes) | SHA-256 `fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90` (matches the README and the reproduction guide) |

`scorch check-input input/tmax_processed.nc --strict` reported 15,738 dates
(1940-04-01 to 2025-09-30), 36 x 50 = 1,800 cells, 10.5-45.5 N, 20.5-69.5 E,
Tmax from -13.24 to 52.51 degC, and a match with the study calendar.

## Commands and results

### The five stages (Catalogs 1-5 and Table 1)

```
scorch run --output-dir <separate folder>       # 157 s in total; stage 1 97 s, stage 3 58 s
```

| Stage | Result | Expected |
|---|---|---|
| 1 | 176,061 local heatwaves | 176,061 |
| 2 | regional threshold 371 cells; 395 selected days | 395 |
| 3 | 753 clusters on 395 days (1-5 per day); 24,885 parameter-search rows | 753 |
| 4 | 753 ellipses | 753 |
| 5 | 51 compound events: Types 1-4 = 3, 3, 20, 25; longest 40 days | 3, 3, 20, 25 |

All nine catalog and table CSVs and the nine supporting CSVs of the stages
(`cell_thresholds.csv`, `daily_summary.csv`, `parameter_selection.csv`,
`daily_settings.csv`, `sequence_settings.csv`, `annual_event_counts.csv`,
`annual_mean_duration.csv`, `duration_trend_statistics.csv`,
`in_text_statistics.csv`) are **byte-identical** to the committed files, and
so is `heatwave_flags.nc`. The same files were produced by the five original
scripts (`scripts/main_workflow/0*.py`), which now call the package.

### Models and appendices (run in place with the committed catalogs)

```
python scripts/modeling/power_law/01_prepare_input.py            #  1 s
Rscript scripts/modeling/power_law/02_fit_power_law.R           # 31 s (2 x 5,000 bootstrap repetitions)
python scripts/modeling/lgcp/01_prepare_input.py                 # 17 s
Rscript scripts/modeling/lgcp/02_fit_lgcp.R                     #  3 s
python scripts/modeling/lgcp/03_centroid_concentration.py        #  1 s
python scripts/modeling/lgcp/validation/01_assign_folds.py       #  1 s
Rscript scripts/modeling/lgcp/validation/02_validate_lgcp.R     # 11 s
python scripts/modeling/lgcp/validation/03_score_validation.py   #  1 s
python scripts/main_workflow/dbscan_parameters/parameter_summary.py   #  2 s
python scripts/main_workflow/ellipse_sensitivity/pca_sensitivity.py   # 14 s
```

| Result | Reproduced value | Supplied value |
|---|---|---|
| Power law, daily maxima (n = 395, tail 43) | A_min = 4,854,155.1583398106 km², alpha = 10.520776128248148, p = 0.2128 | identical |
| Power law, event maxima (n = 51, tail 19) | A_min = 3,585,652.3950546649 km², alpha = 5.3480923823530473, p = 0.3360 | identical |
| LGCP (all 753 centroids) | intercept -14.225999669504311, lon -0.016387583199287341, lat 0.065120560820585205, mean_tmax 0.091388157345756879, std_tmax -0.00050547977736141512, var_sigma2 1.645391601721077, scale_km 277.57792596278443; minimum contrast converged (code 0) | identical (17 significant digits) |
| Five-fold validation | the all-data fit is bitwise identical to the modeling fit; all Appendix C numbers reproduced as printed | identical |
| Appendix B | largest purity sum on the diagonal at scale factor 1.25 (sum 205.056013) | identical |
| Appendix A | 395 days, 51 events, Types 3/3/20/25 | identical |

Every regenerated model and appendix file was byte-identical to the
committed file, including the two `.rds` model objects, with one exception:
`fold_lgcp_parameters.csv` differs only in its `fit_seconds` column (the run
time of each fold's fit). The committed file was kept. No numerical
tolerance was needed anywhere; the comparison is exact.

### Figures and Supporting Information scripts

All sixteen figure scripts (Figures 2, 3, 5-12, A1, A2, B1, C1 and the
PowerPoint exports of Figures 1 and 4) and the station comparison script
(Figure S1) ran without errors through the new package layout. Every figure
and result file they wrote was byte-identical to the committed one, except
the PNGs of Figures 5, 6, 7 and S1: their typology headings use the Abadi
font, which was not supplied in this run (`SCORCH_SLIDE_FONT_DIR` unset), so
they were rendered with the fallback font and differ in pixels only. The
committed PNGs were kept. The six-configuration power-law comparison of the
Supporting Information (about 45 minutes in its separate environment) was
not rerun for this release; its results are the supplied ones.

### Tests

`python -m pytest` with `SCORCH_STUDY_INPUT` set to the archived file:
101 passed in 162 s (81 unit, demo and command-line tests plus the 20 study
regression checks, which run the five stages into a temporary folder and
compare 18 CSV files and the flags NetCDF byte for byte). `Rscript --vanilla
tests/r/test_r_functions.R` passed (11 checks).

## Notes on the validation design

The LGCP validation of Appendix C assigns the 753 centroids to five random
folds (seed 20260704), so the centroids of one compound event usually fall
into several folds. It checks the stability of the fitted concentration
surface; it does not estimate how the model would do on an unseen event or
period. No event-grouped or temporal validation is included in this release;
if one is added later it must be kept separate from the Appendix C results.

## Boundary of the reproduction

The hourly ERA5 preprocessing (daily maxima at 0.25 degrees and the
1-degree cell means) was done outside this repository; the archived NetCDF
file is the starting point of everything reproduced here.
