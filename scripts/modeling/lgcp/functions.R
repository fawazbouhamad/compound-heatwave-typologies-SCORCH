# Shared LGCP fitting, prediction, and file-writing functions.
# Reads: Called by the full-data fit and Appendix C validation.
# Writes: Returns calculations to the calling script.
# Run: Sourced by the calculation scripts.

# Fixed manuscript settings shared by the full-data and validation fits.
lgcp_study_settings <- function() {
  list(lgcp = list(trend_terms = c("lon", "lat", "mean_tmax", "std_tmax"),
                   seed = 20260704L),
       appendix_c = list(n_folds = 5L))
}

find_repo_dir <- function() {
  script_file <- sub("^--file=", "",
                     grep("^--file=", commandArgs(FALSE), value = TRUE))
  folder <- normalizePath(dirname(script_file), winslash = "/")
  while (!file.exists(file.path(folder, "scripts", "main_workflow", "helpers", "functions.py"))) {
    parent <- dirname(folder)
    if (parent == folder) stop("Study repository not found above ", script_file)
    folder <- parent
  }
  folder
}

# Uses the R packages from the study environment.
suppressMessages({
  library(spatstat.geom)
  library(spatstat.model)
  library(spatstat.explore)
  library(jsonlite)
})

# A study path: relative to the repository unless absolute.
in_repo <- function(repo_dir, path) {
  if (grepl("^([A-Za-z]:|/)", path)) path else file.path(repo_dir, path)
}

# Read the prepared inputs and build the spatstat objects of the model.
# Returns a list with
#   centroids  : table of the 753 centroids (one row per ellipse)
#   cells      : table of the 1,800 grid cells with mean_tmax and std_tmax
#   raster     : table of the 25-km raster pixels
#   covariates : spatstat images lon, lat, mean_tmax, std_tmax (NA outside)
#   window     : the study window (spatstat owin)
#   points     : the 753 centroids as a spatstat point pattern (km)
load_lgcp_inputs <- function(repo_dir) {
  input_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "lgcp")
  centroids <- read.csv(file.path(input_dir, "centroids.csv"))
  cells <- read.csv(file.path(input_dir, "temperature_covariates.csv"))
  raster <- read.csv(file.path(input_dir, "covariate_raster_25km.csv"))
  window_info <- fromJSON(file.path(input_dir, "study_window.json"))

  # The raster rows run west to east within a row, rows south to north.
  nx <- window_info$nx
  ny <- window_info$ny
  pixel_x <- sort(unique(raster$x_km))
  pixel_y <- sort(unique(raster$y_km))
  stopifnot(length(pixel_x) == nx, length(pixel_y) == ny, nrow(raster) == nx * ny)
  as_matrix <- function(values) matrix(values, nrow = ny, ncol = nx, byrow = TRUE)
  inside <- as_matrix(raster$inside == 1)
  as_image <- function(values) {
    m <- as_matrix(values)
    m[!inside] <- NA
    im(m, xcol = pixel_x, yrow = pixel_y)
  }

  # Longitude and latitude images use each pixel's own coordinates; the
  # temperature images use the values of the nearest grid cell.
  covariates <- list(lon = as_image(raster$pixel_lon),
                     lat = as_image(raster$pixel_lat),
                     mean_tmax = as_image(raster$mean_tmax),
                     std_tmax = as_image(raster$std_tmax))
  window <- as.owin(as_image(rep(1, nrow(raster))))
  points <- make_point_pattern(centroids$x_km, centroids$y_km, window)
  stopifnot(points$n == nrow(centroids))

  list(centroids = centroids, cells = cells, raster = raster,
       covariates = covariates, window = window, points = points,
       window_info = window_info)
}

# Point pattern of centroid positions (km), allowing duplicate positions
# with checkdup = FALSE.
make_point_pattern <- function(x_km, y_km, window) {
  ppp(x_km, y_km, window = window, checkdup = FALSE)
}

# Fit the LGCP by minimum contrast and collect warnings.
fit_lgcp <- function(points, covariates, trend_terms) {
  model_formula <- as.formula(paste("points ~", paste(trend_terms, collapse = " + ")))
  warnings_seen <- character(0)
  fit <- withCallingHandlers(
    kppm(model_formula, clusters = "LGCP", method = "mincon", statistic = "K",
         model = "exponential", covariates = covariates),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = warnings_seen)
}

# Trend coefficients and Gaussian-field parameters of a fitted model:
# intercept and slopes, the field variance sigma2 and the exponential
# correlation scale (km).
lgcp_parameters <- function(fit) {
  coefficients <- coef(fit)
  field <- as.list(fit$clustpar)
  c(coefficients, var_sigma2 = as.numeric(field$var),
    scale_km = as.numeric(field$scale))
}

# Convergence information of the minimum-contrast step and of the trend fit.
fit_diagnostics <- function(fit) {
  optimisation <- fit$Fit$mcfit$opt
  single <- function(value, missing) {
    if (is.null(value) || length(value) != 1 || is.na(value)) missing else value
  }
  message_text <- if (is.null(optimisation$message)) "" else
    paste(as.character(optimisation$message), collapse = "; ")
  list(mincon_convergence_code = single(as.integer(optimisation$convergence), NA_integer_),
       mincon_objective = single(as.numeric(optimisation$value), NA_real_),
       mincon_n_objfun_evals = single(as.integer(optimisation$counts[["function"]]),
                                      NA_integer_),
       mincon_message = message_text,
       trend_glm_converged = isTRUE(fit$po$internal$glmfit$converged))
}

# Mean intensity (expected centroids per km^2) at the grid-cell centres.
# The fitted trend already contains the -sigma2/2 term of the LGCP, so the
# mean intensity is exp(intercept + sum of slope x covariate); nothing else
# is subtracted. Longitude and latitude are the cell-centre degrees.
mean_intensity_at_cells <- function(fit, cells, trend_terms) {
  coefficients <- coef(fit)
  linear <- coefficients[["(Intercept)"]]
  for (term in trend_terms) linear <- linear + coefficients[[term]] * cells[[term]]
  exp(linear)
}

# Mean intensity (per km^2) as a spatstat image on the 25-km raster.
mean_intensity_image <- function(fit, covariates, trend_terms) {
  coefficients <- coef(fit)
  linear <- covariates$lon * 0 + coefficients[["(Intercept)"]]
  for (term in trend_terms) linear <- linear + coefficients[[term]] * covariates[[term]]
  exp(linear)
}

# Relative concentration rank R(s) of each cell: its rank among all cells
# (ties share the average rank) divided by the number of cells, so the
# values run from 1/1800 to 1. This is pandas' rank(method = "average",
# pct = True), as used for Figure 12.
concentration_rank <- function(intensity) {
  rank(intensity, ties.method = "average") / length(intensity)
}

# Write a table as CSV with LF line endings (a binary connection stops
# Windows from writing CR LF), so file checksums are the same everywhere.
# write.csv keeps 15 significant digits.
write_csv_lf <- function(table, path) {
  connection <- file(path, "wb")
  on.exit(close(connection))
  write.csv(table, connection, row.names = FALSE)
}

# Write lines of text with LF line endings.
write_lines_lf <- function(lines, path) {
  connection <- file(path, "wb")
  on.exit(close(connection))
  writeLines(lines, connection, sep = "\n")
}
