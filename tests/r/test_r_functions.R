# SPDX-License-Identifier: GPL-3.0-only
# Checks of the R model functions with base R only (no test framework needed).
# Run from the repository folder:  Rscript --vanilla tests/r/test_r_functions.R

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/")
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAILED: ", what)
  cat("ok:", what, "\n")
}

# Power-law functions ---------------------------------------------------------
source(file.path(repo_dir, "scripts", "modeling", "power_law", "functions.R"))
set.seed(1)
alpha_true <- 2.5
x <- 100 * (1 - runif(2000))^(-1 / (alpha_true - 1))     # Pareto sample, x_min = 100
fit <- pl_fit(x, ks_mode = "full")
check(identical(fit$status, "ok"), "the power-law fit succeeds on a Pareto sample")
check(abs(fit$alpha - alpha_true) < 0.15, sprintf("alpha %.3f is close to %.1f", fit$alpha, alpha_true))
check(fit$xmin >= min(x) && fit$n_tail > 1, "the cutoff is an observed value with a tail behind it")
boot <- pl_bootstrap(x, fit, reps = 50L, ks_mode = "full", seed = 7L)
check(boot$reps == 50L && boot$n_undefined == 0L, "the bootstrap runs with defined refits")
check(boot$p >= 0 && boot$p <= 1, "the bootstrap p-value is a proportion")
check(identical(format_exact(c(0.1, NA)), c(sprintf("%.17g", 0.1), "NA")), "format_exact")
tied <- c(1, 1, 1, 2, 2, 3, 5, 8, 13)
check(pl_ks(sort(tied), 1, 1.5, "full") >= pl_ks(sort(tied), 1, 1.5, "lower"),
      "the full KS distance is at least the lower-step distance")

# LGCP functions --------------------------------------------------------------
source(file.path(repo_dir, "scripts", "modeling", "lgcp", "functions.R"))
check(requireNamespace("spatstat.model", quietly = TRUE), "spatstat.model is installed")
window <- spatstat.geom::owin(c(0, 200), c(0, 200))
set.seed(3)
points <- spatstat.random::rpoispp(0.005, win = window)
xy <- spatstat.geom::as.im(function(x, y) x / 200, W = window)
covariates <- list(x = xy, y = spatstat.geom::as.im(function(x, y) y / 200, W = window))
result <- fit_lgcp(points, covariates, c("x", "y"))
parameters <- lgcp_parameters(result$fit)
check(all(is.finite(parameters)) && length(parameters) == 5, "the LGCP fit gives five finite parameters")
check(all(names(parameters) == c("(Intercept)", "x", "y", "var_sigma2", "scale_km")), "parameter names")
cells <- data.frame(cell_id = 1:3, x = c(0.1, 0.5, 0.9), y = c(0.5, 0.5, 0.5))
intensity <- mean_intensity_at_cells(result$fit, cells, c("x", "y"))
check(all(intensity > 0), "mean intensities are positive")
check(identical(concentration_rank(c(3, 1, 2)), c(1, 1/3, 2/3)), "concentration_rank")
cat("All R checks passed.\n")
