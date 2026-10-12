# SPDX-License-Identifier: GPL-3.0-only
# Install the R packages used by the SCORCH model scripts.
#
# Usage (from the repository folder):
#   Rscript environment/install_r_packages.R            # current CRAN versions
#   Rscript environment/install_r_packages.R --exact    # the versions of environment/renv.lock
#
# The study used R 4.5.1 with spatstat.geom 3.8-1, spatstat.model 3.7-1,
# spatstat.explore 3.8-1, jsonlite 2.0.0, poweRlaw 1.0.0 and pracma 2.4.6.
# "--exact" installs those versions with the remotes package (source
# compilation may be needed on Linux and macOS; Rtools is needed on Windows
# only when no binary is available).

args <- commandArgs(trailingOnly = TRUE)
exact <- "--exact" %in% args
repos <- c(CRAN = "https://cloud.r-project.org")
packages <- c("spatstat.geom", "spatstat.model", "spatstat.explore", "jsonlite", "poweRlaw", "pracma")

if (!exact) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = repos)
} else {
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", repos = repos)
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
  lock <- jsonlite::fromJSON("environment/renv.lock")
  wanted <- lock$Packages[packages]
  for (name in names(wanted)) {
    version <- wanted[[name]]$Version
    installed <- tryCatch(as.character(packageVersion(name)), error = function(e) NA_character_)
    if (identical(installed, version)) next
    message("installing ", name, " ", version)
    remotes::install_version(name, version = version, repos = repos, upgrade = "never")
  }
}

for (name in packages) {
  cat(sprintf("%-18s %s\n", name, as.character(packageVersion(name))))
}
cat("R packages for SCORCH are ready.\n")
