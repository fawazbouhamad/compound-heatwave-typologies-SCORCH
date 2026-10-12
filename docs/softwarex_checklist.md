# SoftwareX requirements checked against the official sources

Checked on 11 October 2026 against the SoftwareX *Guide for Authors*
(https://www.sciencedirect.com/journal/softwarex/publish/guide-for-authors)
and the *Original Software Publication* template, version 6, March 2026
(https://legacyfileshare.elsevier.com/promis_misc/softwarex-osp-template.docx,
LaTeX version `softwarex-osp-template.tex`). Quotations are from those
documents.

## What the journal requires of the repository

| Requirement (source) | Where SCORCH meets it |
|---|---|
| "Store your code/software in an open-access GitHub repository where it is publicly available" (guide) | https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH, public |
| "README.md and LICENSE.txt files. Your README.md file should be well-written, explaining the use, installation, and purpose of your proposed software" (guide); "GitHub repositories must have a well-documented README.md file, and a Licence.txt file" (template) | `README.md` (purpose, installation, demo, use, reproduction, citation) and `LICENSE.txt` (GPL-3.0-only) |
| "Source code in a repo/src directory" (guide) | `src/scorch/` holds the package; `scripts/` holds the thin manuscript entry points that call it |
| "All software and code must be open source and released under a recognized legal license such as those listed ... and by the Open Source Initiative", examples "Apache License, GNU General Public License (GPL) or MIT License" (guide, template) | GPL-3.0-only, unchanged from the first commit |
| "We only accept GitHub repositories within Metadata Tables" (template) | the C2 link below is the GitHub repository |
| Research data: "Deposit your research data in a relevant data repository. Cite and link to this dataset in your article" (guide, Option C) | processed input on Zenodo, DOI 10.5281/zenodo.21717874, cited in `CITATION.cff` and the README |
| "Declaration of generative AI use ... Authors must declare the use of generative AI tools in the manuscript preparation process" (guide) | applies to the article; the repository records the AI assistance used for the software packaging in `CHANGELOG.md` and the README so the declaration can be written accurately |
| Software citation: "If the software repository you used supplied a DOI or another Persistent IDentifier (PID), please add a reference for your software" (template) | `CITATION.cff`, `codemeta.json` and `.zenodo.json` are ready; the software DOI is added after the Zenodo archive of the GitHub release (see `release_checklist.md`) |

The guide does not prescribe a build system, a test framework or a particular
documentation format beyond the README; the template limits the article to
4,000 words, six figures and the five mandatory sections (Motivation and
significance; Software description; Illustrative examples; Impact;
Conclusions).

## Code metadata table (template Table 1, "Code metadata (mandatory)")

| Nr | Code metadata description | Value for SCORCH |
|---|---|---|
| C1 | Current code version | v1.0.0 |
| C2 | Permanent link to code/repository used for this code version | https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH (release tag `v1.0.0`) |
| C3 | Legal code license | GNU General Public License v3.0 only (GPL-3.0-only) |
| C4 | Code versioning system used | git |
| C5 | Software code languages, tools and services used | Python 3.11+ (numpy, pandas, scipy, scikit-learn, netCDF4, pyshp, PyYAML; matplotlib, Cartopy, Shapely, pyproj for figures and model input), R 4.5 (spatstat.geom, spatstat.model, spatstat.explore, jsonlite, poweRlaw), GitHub Actions |
| C6 | Compilation requirements, operating environments and dependencies | no compilation; Windows 11 verified, Linux and macOS exercised by continuous integration; dependencies in `pyproject.toml`, `requirements/*.txt` and `environment/renv.lock` |
| C7 | If available, link to developer documentation/manual | https://github.com/fawazbouhamad/compound-heatwave-typologies-SCORCH/tree/main/docs |
| C8 | Support email for questions | the corresponding author's address (to be entered in the submission) |

## Illustrative example for the article

`scorch demo` (docs/demo.md) is the illustrative example: a synthetic
dataset with known events, the command, its printed comparison table and
the overview figure. The figure is clearly labelled synthetic; it is not a
climate result.

## Points to state honestly in the article

- The ERA5 preprocessing (hourly to daily maximum, 0.25 to 1 degree) was done
  outside this repository; the software starts from the archived processed
  file.
- The software has been used for one study so far (the manuscript in
  preparation); there are no download or user statistics to report.
- The LGCP validation of Appendix C splits the 753 centroids at random, so
  centroids of one event can fall into several folds; it is a check of the
  fitted surface, not of generalisation to unseen events.
