# Input and output reference

## The processed temperature input (`tmax_processed.nc`)

One NetCDF-4 file with the daily maximum 2-m air temperature of every grid
cell on every warm-season date.

| Item | Requirement |
|---|---|
| Dimensions | `time`, `lat`, `lon`, in this order for `tmax` |
| `tmax(time, lat, lon)` | float (float64 in the study file), attribute `units = "degC"`; every value present and finite (SCORCH does not fill gaps) |
| `time(time)` | CF time coordinate: integer or float with `units` such as `days since 1940-01-01 00:00:00` and `calendar = "standard"` (the default when the attribute is missing); one value per date, no repeats |
| `lat(lat)`, `lon(lon)` | cell-centre coordinates in degrees, strictly increasing |
| Calendar | exactly the warm-season days of the study: every day of `study_area.warm_season_months` from `first_year` to `last_year`, in order (the study: 1 April to 30 September of 1940-2025, 15,738 dates) |
| Grid | a regular grid; the study uses 1 degree x 1 degree cells with centres at .5 degrees (36 x 50 = 1,800 cells, 10.5-45.5 N, 20.5-69.5 E) |
| Missing data | none allowed; days outside the warm season are simply absent from `time` |

The study file (Zenodo record 10.5281/zenodo.21717874, version 1.0.0, file
`SCORCH_data.zip`, MD5 `80772df9410306cbc419191d0dcb824c`; inside it
`scorch_data/tmax_processed.nc`, SHA-256
`fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90`) was
prepared from hourly ERA5 2-m temperature: the daily maximum of the 24 hourly
values of each UTC calendar day at 0.25 degrees, then the unweighted mean of
the 16 daily maxima of the 0.25-degree points in each 1-degree cell. That
preprocessing was done outside this repository (its global attributes
describe it); SCORCH starts from the prepared file. `scorch check-input
FILE --strict` checks a file against these requirements.

## The grid-cell table (`maps/grid_cells.csv`)

Lists every cell of the input grid; stage 1 checks that it matches the
NetCDF coordinates.

| Column | Meaning |
|---|---|
| `cell_id` | 1, 2, ... by latitude row (south to north), then longitude (west to east) |
| `lat`, `lon` | cell-centre coordinates (degrees) |
| `lat_index`, `lon_index` | 0-based position of the cell in the `lat` and `lon` coordinates |
| `lat_south`, `lat_north`, `lon_west`, `lon_east` | cell edges |

`scorch.io.grid_cells_table(lat, lon)` writes this table for any regular
grid; `scorch demo` writes one for the synthetic grid.

## Dates and identifiers

All dates are ISO strings (`YYYY-MM-DD`). Identifiers are integers starting
at 1 and are assigned in a fixed order so that two runs on the same input
give the same files:

| Identifier | Assigned by | Order |
|---|---|---|
| `local_event_id` | stage 1 | by `cell_id`, then start date |
| `cluster_id` = `ellipse_id` | stages 3-4 | by date, then DBSCAN label (0, 1, 2, ... in the order the clusters were found, west to east) |
| `compound_event_id` | stage 5 | by the first date of the event |

## Catalog 1: local heatwaves (`output/catalogs/01_heatwave_events/heatwave_events.csv`)

One row per local heatwave at one grid cell (the study: 176,061 rows).

| Column | Meaning |
|---|---|
| `local_event_id`, `cell_id`, `lat`, `lon` | identifiers and cell centre |
| `event_number_in_cell` | 1, 2, ... within the cell through the record |
| `start_date`, `end_date`, `duration_days` | first and last day (both exceeding) and the span |
| `exceeding_days` | days with Tmax at or above the cell threshold |
| `bridged_days` | single non-exceeding days inside the heatwave |
| `peak_date`, `peak_tmax_c` | day of the highest Tmax (the earlier day on a tie) and its value in degC |

Supporting results: `heatwave_detection/cell_thresholds.csv` (threshold of
every cell in degC and its numbers of exceeding days, heatwave days and local
events) and `heatwave_detection/heatwave_flags.nc` (daily `exceedance`,
`heatwave` and `local_event_id` on the grid; 0/1 flags, no missing values).

## Catalog 2: regionally extensive days (`output/catalogs/02_regionally_extensive_heatwaves/`)

`selected_days.csv`, one row per selected day (the study: 395):

| Column | Meaning |
|---|---|
| `date` | the selected day |
| `n_heatwave_cells` | cells in a local heatwave that day |
| `heatwave_fraction_of_domain` | that number divided by the number of cells |
| `n_exceeding_cells` | cells at or above their threshold that day |
| `n_heatwave_and_exceeding_cells` | cells that are both |
| `heatwave_share_of_exceeding_cells` | the ratio of the last two (empty when no cell exceeds) |

`selected_day_cells.csv` lists the heatwave cells of every selected day
(`date`, `cell_id`, `lat`, `lon`, `local_event_id`), west to east and south
to north within each day. The supporting file
`regional_selection/daily_summary.csv` has the same counts for every study
day plus `regional_threshold_cells` and `selected` (0/1).

## Catalog 3: clusters (`output/catalogs/03_extreme_heatwave_clusters/`)

`clusters.csv`, one row per DBSCAN cluster (the study: 753):

| Column | Meaning |
|---|---|
| `cluster_id`, `date`, `cluster_number_in_day`, `dbscan_label` | identifiers; clusters are numbered 0, 1, ... within the day |
| `n_member_cells`, `n_core_cells` | cells in the cluster and DBSCAN core cells among them |
| `eps_applied`, `min_samples_applied` | the settings applied to the day (those of its run of consecutive selected days) |
| `n_clusters_on_day`, `n_noise_cells_on_day` | clusters found on the day and heatwave cells that belong to none |

`cluster_cells.csv` lists every heatwave cell of every selected day with its
`dbscan_label`, `cluster_id` (empty for noise cells), `is_noise`,
`is_core_cell` and `n_cells_within_eps`. Supporting results in `clustering/`:
`parameter_selection.csv` (the number of clusters for every day and every
`eps`/`min_samples` pair, with `retained` = 1 for the kept pairs),
`daily_settings.csv` (the modal counts and daily settings) and
`sequence_settings.csv` (the settings of every run of consecutive days).

## Catalog 4: ellipses (`output/catalogs/04_extreme_heatwave_ellipses/ellipses.csv`)

One row per cluster (`ellipse_id` = `cluster_id`; the study: 753).

| Column | Meaning |
|---|---|
| `n_member_cells` | cells in the cluster |
| `centroid_lon_unweighted`, `centroid_lat_unweighted` | mean position of the member cells (ellipse centre) |
| `eigenvalue_major_km2`, `eigenvalue_minor_km2` | eigenvalues of the sample covariance of the projected positions |
| `major_axis_km`, `minor_axis_km` | full axes `2 * 1.25 * sqrt(max(eigenvalue, 1 km2))` |
| `area_km2` | `pi * major * minor / 4` |
| `axis_ratio` | minor divided by major, in (0, 1] |
| `orientation_east_ccw_deg` | major-axis angle counterclockwise from east, in [0, 180) |
| `azimuth_north_cw_deg` | the same axis clockwise from north, in (-90, 90] (an east-west axis is +90) |
| `major_axis_east`, `major_axis_north` | unit vector of the major axis |
| `centroid_lon_weighted`, `centroid_lat_weighted` | Tmax-weighted mean position of the member cells |
| `tmax_weight_sum_c`, `tmax_min_c`, `tmax_max_c` | sum, minimum and maximum of the member-cell Tmax (degC; single precision, as in the study) |
| `scale_factor`, `eigenvalue_floor_km2` | the settings used |

Positions are projected with a local equirectangular projection (one degree
of latitude = 110.574 km, one degree of longitude = 111.320 km times the
cosine of the mean latitude).

## Catalog 5: compound events (`output/catalogs/05_compound_heatwave_events/`)

`compound_events.csv`, one row per compound event (the study: 51):

| Column | Meaning |
|---|---|
| `compound_event_id`, `start_date`, `end_date`, `duration_days`, `n_selected_days` | identifiers and span (every day of an event is a selected day) |
| `event_type`, `event_type_name` | 1-4 and its name (see below) |
| `daily_cluster_counts` | clusters on each day, for example `1-2-2-2` |
| `total_clusters`, `min_daily_clusters`, `max_daily_clusters`, `n_single_cluster_days`, `n_multi_cluster_days` | counts behind the type |
| `event_eps`, `event_min_samples` | the DBSCAN settings applied to the event's days |
| `mean_area_km2`, `mean_axis_ratio` | means over the event's ellipses |
| `axial_mean_azimuth_deg`, `axial_resultant_length` | doubled-angle mean orientation of the ellipses and its resultant length (1 = all parallel) |
| `largest_ellipse_id`, `largest_ellipse_date`, `largest_ellipse_area_km2` | the event's largest ellipse (the earlier one on a tie) |

`event_days.csv` maps every selected day to its event (`date`,
`compound_event_id`, `day_in_event`).

### The four event types

| Type | Name | Days | Clusters per day |
|---|---|---|---|
| 1 | Widespread (Isolated) | 1 | 1 |
| 2 | Spatially Clustered | 1 | 2 or more |
| 3 | Temporally Clustered | 2 or more | every day 1, or every day 2 or more |
| 4 | Compound Clustering (Multi-Type) | 2 or more | a mix of 1-cluster days and several-cluster days |

## Table 1 (`output/tables/table_1.csv`)

One row per type: `selected_event_days`, `number_of_compound_events`,
`longest_event_days`, `mean_duration_days_per_event`,
`mean_ellipses_per_event`, `minimum_ellipses_per_event`,
`maximum_ellipses_per_event`. A type without events has 0 events and empty
statistics.

## Other supporting results of stage 5

`annual_event_counts.csv` (events per year and type, every year from the
first event to the last study year), `annual_mean_duration.csv` and
`duration_trend_statistics.csv` (OLS, Sen's slope and Mann-Kendall tests on
the annual mean duration of Type 3 and Type 4 events; a type with fewer than
three years with events gets a row of missing values) and, for the study
region only, `in_text_statistics.csv` (the numbers quoted in Sections 3.1 and
4 of the manuscript; switched off with `summaries.in_text_statistics: false`).

## Run records

Every stage writes `<run_records_dir>/<stage>.json` with the time it
finished, the Python and package versions, and the SHA-256 checksum of every
file it read and wrote.

## Model inputs and outputs (`supporting_files/additional_results/modeling/`)

The R scripts read the catalogs and write: `power_law/daily_max_area.csv`
and `event_max_area.csv` (the largest ellipse area per selected day and per
event), `power_law/power_law_fits.csv` (cutoff, exponent, KS distance and
bootstrap p-value of the two fits) with the bootstrap replicates;
`lgcp/centroids.csv`, `temperature_covariates.csv`, `covariate_raster_25km.csv`
and `study_window.json` (the LGCP input), `lgcp/lgcp_parameters.csv` (the
seven fitted parameters), the predicted intensities, the centroid
concentration summaries, and the five-fold validation tables in
`lgcp/validation/`. The columns are named as in the manuscript's Appendix C.
