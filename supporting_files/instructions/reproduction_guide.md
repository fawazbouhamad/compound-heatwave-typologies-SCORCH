# Reproduction guide

See the [repository overview and dataset information](../../README.md#processed-temperature-dataset)
for the study background and temperature download.

## Reproduce the manuscript results

Follow Steps 1–6 in order to generate **Catalogs 1–5, Table 1, Figures 1–12
and Appendix Figures A1, A2, B1 and C1**. The following section covers the
manuscript's **Supporting Information, Figures S1–S3**.

Run the commands one at a time in **PowerShell**, from the repository's root
folder. The scripts overwrite the corresponding supplied results, so use a
fresh working copy and stop if a command reports an error.

### 1. Download the repository and place the temperature file

1. Download this repository using **Code → Download ZIP**, then extract it,
   or clone it with Git.
2. Download [SCORCH_data.zip](https://zenodo.org/records/21717874/files/SCORCH_data.zip?download=1)
   from the Zenodo dataset record (version 1.0.0,
   [DOI 10.5281/zenodo.21717874](https://doi.org/10.5281/zenodo.21717874)).
   Its MD5 checksum, shown on Zenodo, is `80772df9410306cbc419191d0dcb824c`
   (`Get-FileHash SCORCH_data.zip -Algorithm MD5`). Extract it and, inside its
   `scorch_data` folder, locate `tmax_processed.nc`.
3. In the repository's root folder, create a folder named `input`.
4. Copy `tmax_processed.nc` into `input`. The final path must be
   **`input/tmax_processed.nc`**; do not copy the enclosing `scorch_data` folder.
5. Open PowerShell in the repository's root folder, which contains
   `README.md`, `scripts/` and `output/`.

Check that the downloaded file is the current (corrected) study input:

```powershell
Get-FileHash input/tmax_processed.nc -Algorithm SHA256
```

The expected SHA-256 is:

```text
fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90
```

### 2. Prepare Python, R and the figure tools

The verified runs used **Windows, Python 3.12.10 and R 4.5.1**. Obtain the
main package-version and installation record from the authors, and use it to
prepare the required Python and R packages. This record is maintained
separately from the repository. Other operating systems have not been verified.

The commands below assume the main Python environment is named `.venv` in the
repository root. Set these interpreter paths in the same PowerShell session,
adjusting them to your installation:

```powershell
$py = ".\.venv\Scripts\python.exe"
$r = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
```

Study settings are already recorded in the [Python helpers](../../scripts/main_workflow/helpers/functions.py)
and [R model scripts](../../scripts/modeling/). No setting changes are needed to
reproduce the manuscript analysis.

For matching plot typography, supply the two Abadi font files listed in the
[font inventory](../additional_inputs/figure_assets/font_files.csv).
They are not included in the repository. Replace the example font-folder path
below with their actual location before running these commands:

```powershell
$env:SCORCH_SLIDE_FONT_DIR = "C:\path\to\Abadi"
$env:SCORCH_REQUIRE_SLIDE_FONTS = "1"
```

Without those fonts, omit these two settings; substitute fonts can change the
appearance of plots while the calculations remain the same.

**Figures 1 and 4 also require desktop Microsoft PowerPoint on Windows**, with
Aptos, Abadi and Cambria Math available to PowerPoint. Their editable slides
are included in [figure_assets](../additional_inputs/figure_assets/).
The font-folder setting above applies to Python plots; PowerPoint uses its
installed fonts.

### 3. Generate the five catalogs and Table 1

Run the five analysis scripts in this order. Each stage uses the results of
the earlier stages.

```powershell
& $py -B scripts/main_workflow/01_detect_heatwaves.py
& $py -B scripts/main_workflow/02_select_regional_days.py
& $py -B scripts/main_workflow/03_cluster_heatwaves.py
& $py -B scripts/main_workflow/04_fit_ellipses.py
& $py -B scripts/main_workflow/05_summarize_compound_events.py
```

| Stage | Catalog generated | Expected contents |
|---|---|---|
| 1 | [Local heatwave events](../../output/catalogs/01_heatwave_events/) | 176,061 heatwaves at individual grid cells |
| 2 | [Regionally extensive heatwaves](../../output/catalogs/02_regionally_extensive_heatwaves/) | 395 selected dates and their heatwave cells |
| 3 | [Extreme heatwave clusters](../../output/catalogs/03_extreme_heatwave_clusters/) | 753 spatial clusters and their cell membership |
| 4 | [Extreme heatwave ellipses](../../output/catalogs/04_extreme_heatwave_ellipses/) | Geometry and temperature-weighted centroids of 753 ellipses |
| 5 | [Compound heatwave events](../../output/catalogs/05_compound_heatwave_events/) | 51 classified events and their constituent dates |

The fifth script also writes [Table 1](../../output/tables/table_1.csv), summarizing
the four compound-heatwave types. Compound-event identifiers are assigned in
Catalog 5. Intermediate calculations are saved in
[additional_results](../additional_results/).

### 4. Generate Figures 1–10

Run the figure scripts after generating the catalogs:

```powershell
& $py -B scripts/main_workflow/figure_generation/figure_01.py
& $py -B scripts/main_workflow/figure_generation/figure_02.py
& $py -B scripts/main_workflow/figure_generation/figure_03.py
& $py -B scripts/main_workflow/figure_generation/figure_04.py
& $py -B scripts/main_workflow/figure_generation/figure_05.py
& $py -B scripts/main_workflow/figure_generation/figure_06.py
& $py -B scripts/main_workflow/figure_generation/figure_07.py
& $py -B scripts/main_workflow/figure_generation/figure_08.py
& $py -B scripts/main_workflow/figure_generation/figure_09.py
& $py -B scripts/main_workflow/figure_generation/figure_10.py
```

Each figure is saved in its own folder under [output/figures](../../output/figures/).
Complete figures retain their panel letters; separately exported panels omit
them. The scripts for Figures 1 and 4 read the supplied PowerPoints and export
PNGs directly to these output folders.

### 5. Fit the models and generate Figures 11–12

First, prepare the maximum-area samples, fit the power-law models and draw
Figure 11:

```powershell
& $py -B scripts/modeling/power_law/01_prepare_input.py
& $r --vanilla scripts/modeling/power_law/02_fit_power_law.R
& $py -B scripts/modeling/power_law/03_make_figure_11.py
```

Next, prepare and fit the log-Gaussian Cox process (LGCP), calculate centroid
concentration and draw Figure 12:

```powershell
& $py -B scripts/modeling/lgcp/01_prepare_input.py
& $r --vanilla scripts/modeling/lgcp/02_fit_lgcp.R
& $py -B scripts/modeling/lgcp/03_centroid_concentration.py
& $py -B scripts/modeling/lgcp/04_make_figure_12.py
```

Keep this order: LGCP preparation reads the power-law input tables to identify
the largest ellipse per selected day and per event. These identifiers label
comparison groups; **the LGCP fit uses all 753 centroids**.

### 6. Reproduce the appendix analyses and figures

**Appendix A — DBSCAN parameter summaries and Figures A1–A2:**

```powershell
& $py -B scripts/main_workflow/dbscan_parameters/parameter_summary.py
& $py -B scripts/main_workflow/dbscan_parameters/figure_A1.py
& $py -B scripts/main_workflow/dbscan_parameters/figure_A2.py
```

**Appendix B — ellipse sensitivity and Figure B1:**

```powershell
& $py -B scripts/main_workflow/ellipse_sensitivity/pca_sensitivity.py
& $py -B scripts/main_workflow/ellipse_sensitivity/figure_B1.py
```

**Appendix C — LGCP validation and Figure C1:**

```powershell
& $py -B scripts/modeling/lgcp/validation/01_assign_folds.py
& $r --vanilla scripts/modeling/lgcp/validation/02_validate_lgcp.R
& $py -B scripts/modeling/lgcp/validation/03_score_validation.py
& $py -B scripts/modeling/lgcp/validation/04_make_figure_C1.py
```

The appendix figures are saved under [output/figures](../../output/figures/).
Their calculations and model results are retained in
[additional_results](../additional_results/).

## Reproduce the Supporting Information

The manuscript's Supporting Information contains **Text S1** and **Figures
S1–S3**. Its scripts, station observations, results and figures are in
[supporting_information](../supporting_information/).

- **Figure S1:** comparison of ERA5 temperatures with GHCN-Daily observations
  at Aswan, Egypt.
- **Text S1 and Figures S2–S3:** comparison of six power-law fitting
  configurations, using the largest ellipse per selected day (S2) and per
  compound event (S3).

Continue after the manuscript steps, in the same PowerShell session. The
power-law comparison uses different Python package versions, so prepare a
separate environment using its included
[requirements file](../supporting_information/power_law_comparison/requirements.txt):

```powershell
py -3.12 -m venv .venv-powerlaw
$powerlaw = ".\.venv-powerlaw\Scripts\python.exe"
& $powerlaw -m pip install -r supporting_files/supporting_information/power_law_comparison/requirements.txt
```

Then run the station comparison, the R and Python power-law comparisons, and
the figure script:

```powershell
& $py -B supporting_files/supporting_information/station_comparison/station_comparison.py
& $r --vanilla supporting_files/supporting_information/power_law_comparison/01_compare_r.R
& $powerlaw -B supporting_files/supporting_information/power_law_comparison/02_compare_python.py
& $powerlaw -B supporting_files/supporting_information/power_law_comparison/03_make_figures.py
```

The resulting figures are saved in
[supporting_information/figures](../supporting_information/figures/).
Tables and fit results are saved in the respective comparison folders. These
scripts reproduce the analyses behind Text S1; they do not generate the
Supporting Information document itself.

## Check the reproduced results

Compare your outputs with the supplied catalogs, Table 1 and figure files.
The catalog counts above provide a first check. The recorded study runs took
about **5 minutes for the manuscript analyses and figures** and **45 minutes
for the Supporting Information comparisons** on the study computer. The
repeated bootstrap fits account for most of the latter; timings vary by machine.

The supplied results were generated on **6 October 2026** by a complete Windows
run of Steps 3–6 and the Supporting Information on the corrected input
(SHA-256 above). Their catalogs, Table 1, trend statistics, power-law fits,
LGCP parameters and centroid concentration agree exactly with an independent
implementation of the method (the scoRch.generator R package) run on the same
input. Figures 1 and 4 are schematic diagrams that do not depend on the data;
they were not regenerated and are exported from the editable PowerPoints.

The historical first analysis used an earlier version of the input whose
45.5°N cells averaged only the 0.25° points at 45.0°N. Its outputs, with
175,978 local heatwaves, 760 clusters and ellipses and 3, 4, 20 and 24
compound events of Types 1–4, are retained in the authors' archive.

Timing fields, stored script paths and text-file line endings can differ
between runs. Fonts and PowerPoint versions can also affect rendered pixels.

## Comparisons used to inform study decisions

The [extra_analysis](../extra_analysis/) folder preserves earlier work that
helped guide method choices but is not presented as figures or tables in the
final manuscript. These files are available for context; the manuscript
reproduction steps do not read or rerun them.

| Comparison | Decision examined |
|---|---|
| [Centroid weighting](../extra_analysis/centroid_weighting/) | Unweighted, Tmax, exceedance, percentile-above-95th and hybrid weighting; includes [four animations](../extra_analysis/centroid_weighting/gifs/) for 1–21 August 2010 |
| [Regional-day threshold](../extra_analysis/regional_day_threshold/regional_day_selection_sensitivity.csv) | How percentile cutoffs and quantile methods affect the selected dates |
| [DBSCAN distance choices](../extra_analysis/dbscan_eps_grid_comparison/) | Two candidate neighbourhood-distance lists within DBSCAN, tested on the same 395 days |
| [Clustering methods](../extra_analysis/clustering_method_comparison/) | Connected components, DBSCAN, OPTICS and HDBSCAN in an earlier workflow |
| [LGCP temperature covariates](../extra_analysis/lgcp_covariate_comparison/) | Four combinations of mean Tmax, coefficient of variation and standard deviation, with longitude and latitude included throughout |

## Other supporting files

- [additional_inputs](../additional_inputs/): study-grid
  definitions, geographic map layers, PowerPoint figure sources and the font
  inventory.
- [additional_results](../additional_results/): thresholds,
  daily summaries, clustering and ellipse diagnostics, model fits and
  validation results generated by the study scripts.
