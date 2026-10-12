# Installation

SCORCH has two parts: a Python package (the five catalog stages, the demo,
the tests and most figures) and a set of R scripts (the power-law and
log-Gaussian Cox process models and their validation). You need Python for
everything; you need R only for the model scripts.

## Supported systems

| | Verified | Expected to work |
|---|---|---|
| Operating system | Windows 11 (the study computer) | Linux and macOS: the continuous-integration workflow (`.github/workflows/ci.yml`) installs the same pinned versions on Ubuntu and macOS and runs the tests and the demo |
| Python | CPython 3.12.10 | 3.11 or newer (`tomllib` and pandas 3 need 3.11) |
| R | R 4.5.1 | R 4.3 or newer |

Figures 1 and 4 of the manuscript are exported from editable PowerPoint
files and need desktop Microsoft PowerPoint on Windows; nothing else does.
Some figure headings use the Abadi font; without it a fallback font is used
and only the look of the plots changes.

## 1. Python

Create a virtual environment and install the package. The `requirements`
folder pins the exact package versions of the verified study run; the
`pyproject.toml` lists the minimum versions.

```bash
git clone https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH.git
cd compound-heatwave-typologies-SCORCH
python -m venv .venv
```

Activate the environment (`.venv\Scripts\activate` on Windows,
`source .venv/bin/activate` on Linux and macOS), then install:

```bash
python -m pip install -r requirements/core.txt
python -m pip install -e .
```

That is enough for the five stages, the demo (without its figure) and the
tests. Three optional groups of packages add the rest:

| Group | Install | Needed for |
|---|---|---|
| figures | `python -m pip install -r requirements/figures.txt` | the manuscript figures (matplotlib, Cartopy, Shapely, Pillow), the demo overview figure, and the LGCP input preparation (pyproj) |
| development | `python -m pip install -r requirements/dev.txt` | running the tests with pytest |
| supplementary power-law comparison | `python -m pip install -r supporting_files/supporting_information/power_law_comparison/requirements.txt` into a **separate** environment (it pins other versions of numpy and pandas) | Text S1 and Figures S2-S3 of the Supporting Information |

The same groups are available as extras: `python -m pip install -e ".[figures,dev]"`
installs the minimum versions instead of the pinned ones.

Check the installation:

```bash
scorch info
scorch demo
```

`scorch demo` builds a small synthetic example and runs the whole workflow
in a few seconds (see [demo.md](demo.md)).

If you prefer conda, create an environment with Python 3.12
(`conda create -n scorch python=3.12`), activate it and run the same `pip`
commands inside it.

## 2. R (model scripts only)

Install R 4.5 from [CRAN](https://cran.r-project.org/), then install the
packages from the repository folder:

```bash
Rscript environment/install_r_packages.R
```

This installs the current CRAN versions of `spatstat.geom`, `spatstat.model`,
`spatstat.explore`, `jsonlite`, `poweRlaw` and `pracma`. To install exactly
the versions of the study (recorded in `environment/renv.lock`: spatstat.geom
3.8-1, spatstat.model 3.7-1, spatstat.explore 3.8-1, jsonlite 2.0.0, poweRlaw
1.0.0, pracma 2.4.6), run `Rscript environment/install_r_packages.R --exact`.
Users of [renv](https://rstudio.github.io/renv/) can also call
`renv::restore(lockfile = "environment/renv.lock")`.

Check the R side with:

```bash
Rscript --vanilla tests/r/test_r_functions.R
```

On Windows, `Rscript` is usually not on the `PATH`; use its full path (for
example `"C:\Program Files\R\R-4.5.1\bin\Rscript.exe"`) or set the
environment variable `SCORCH_RSCRIPT` to it. `scorch demo --with-r` finds
Rscript through `SCORCH_RSCRIPT`, the `PATH` or the usual Windows install
folders.

## 3. The study input (only for reproducing the study)

The processed temperature file is not in the repository. Download
`SCORCH_data.zip` from the Zenodo dataset record (version 1.0.0, DOI
[10.5281/zenodo.21717874](https://doi.org/10.5281/zenodo.21717874)), check
its MD5 (`80772df9410306cbc419191d0dcb824c`), extract it and copy
`scorch_data/tmax_processed.nc` to `input/tmax_processed.nc`. Then:

```bash
scorch check-input input/tmax_processed.nc --strict
```

The expected SHA-256 of the NetCDF file is
`fc080880b6f2065fb5eec9fe9dce71abdfdc44a72d1ae91e50abdd278fadaa90`.
See [reproduction_record.md](reproduction_record.md) for the full
reproduction and the verification results.

## Clean-environment check

The pinned Python environment was rebuilt from scratch on Windows 11 with
`py -3.12 -m venv` and the two `requirements` files on 11 October 2026, and
reproduced the published catalogs exactly. The GitHub Actions workflow repeats
the installation on Ubuntu, Windows and macOS runners for every change.
