# SPDX-License-Identifier: GPL-3.0-only
# Fit the manuscript's power-law procedure to the demo ellipse areas.
# Usage: Rscript --vanilla demo/fit_power_law_demo.R <ellipses.csv> <output.csv>
# This only shows that the R side of SCORCH is installed and runs; the demo
# areas are synthetic and far too few for a meaningful tail fit.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript demo/fit_power_law_demo.R <ellipses.csv> <output.csv>")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source(file.path(repo_dir, "scripts", "modeling", "power_law", "functions.R"))

areas <- read.csv(args[1])$area_km2
fit <- pl_fit(areas, ks_mode = "full")
if (!identical(fit$status, "ok")) stop("the power-law fit did not succeed: ", fit$status)
boot <- pl_bootstrap(areas, fit, reps = 200L, ks_mode = "full", seed = 20261011L)
result <- data.frame(n = fit$n, n_tail = fit$n_tail, xmin_km2 = fit$xmin, alpha = fit$alpha,
                     ks_distance = fit$D, bootstrap_reps = boot$reps,
                     p_value = boot$p, n_undefined = boot$n_undefined)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
connection <- file(args[2], "wb")
write.csv(result, connection, row.names = FALSE)
close(connection)
cat(sprintf(paste0("R %s with %d ellipse areas: cutoff A_min = %.0f km2, tail of %d areas, ",
                   "exponent alpha = %.3f, KS distance = %.3f, bootstrap p = %.3f (%d repetitions)\n"),
            getRversion(), fit$n, fit$xmin, fit$n_tail, fit$alpha, fit$D, boot$p, boot$reps))
cat(sprintf("wrote %s\n", args[2]))
