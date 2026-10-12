# Fit the five validation models using the shared LGCP functions.
# Reads: Prepared LGCP inputs, fold assignments, and the full-data LGCP fit.
# Writes: supporting_files/additional_results/modeling/lgcp/validation/
# Run: Rscript scripts/modeling/lgcp/validation/02_validate_lgcp.R

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_dir <- normalizePath(file.path(dirname(script_file), "..", "..", "..", ".."), winslash = "/")
source(file.path(repo_dir, "scripts", "modeling", "lgcp", "functions.R"))

config <- lgcp_study_settings()
output_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "lgcp", "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trend_terms <- unlist(config$lgcp$trend_terms)
n_folds <- config$appendix_c$n_folds
# The minimum-contrast fit is deterministic. Set the configured seed for
# reproducibility if any package operation uses random numbers.
set.seed(config$lgcp$seed)

inputs <- load_lgcp_inputs(repo_dir)
cells <- inputs$cells
folds <- read.csv(file.path(output_dir, "fold_assignment.csv"))
stopifnot(nrow(folds) == nrow(inputs$centroids),
          identical(folds$ellipse_id, inputs$centroids$ellipse_id))

parameter_rows <- list()
full_precision_rows <- list()
surfaces <- list()
fits <- list()
log_lines <- character(0)

describe_fit <- function(label, fold, points, result, seconds) {
  parameters <- lgcp_parameters(result$fit)
  diagnostics <- fit_diagnostics(result$fit)
  data.frame(
    fit = label, fold = fold, n_training_centroids = points$n,
    n_withheld_centroids = nrow(folds) - points$n,
    intercept = parameters[["(Intercept)"]], beta_lon = parameters[["lon"]],
    beta_lat = parameters[["lat"]], beta_mean_tmax = parameters[["mean_tmax"]],
    beta_std_tmax = parameters[["std_tmax"]], var_sigma2 = parameters[["var_sigma2"]],
    scale_km = parameters[["scale_km"]],
    mincon_convergence_code = diagnostics$mincon_convergence_code,
    mincon_objective = diagnostics$mincon_objective,
    mincon_n_objfun_evals = diagnostics$mincon_n_objfun_evals,
    trend_glm_converged = diagnostics$trend_glm_converged,
    n_warnings = length(result$warnings), fit_seconds = seconds,
    stringsAsFactors = FALSE)
}

# The seven parameters with 17 significant digits (enough to give back every
# double exactly), so that they can be compared bit for bit without reading
# the model files.
full_precision_table <- function(label, fold, fit) {
  parameters <- lgcp_parameters(fit)
  data.frame(fit = label, fold = fold, parameter = names(parameters),
             value = sprintf("%.17g", parameters), stringsAsFactors = FALSE)
}

# one fit per fold, on the training centroids only
for (k in seq_len(n_folds)) {
  training <- folds$fold != k
  training_points <- make_point_pattern(folds$x_km[training], folds$y_km[training],
                                        inputs$window)
  stopifnot(training_points$n == sum(training))
  started <- Sys.time()
  result <- fit_lgcp(training_points, inputs$covariates, trend_terms)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  parameter_rows[[k]] <- describe_fit("fold", k, training_points, result, seconds)
  full_precision_rows[[k]] <- full_precision_table("fold", k, result$fit)
  surfaces[[k]] <- data.frame(fold = k, cell_id = cells$cell_id,
                              mean_intensity_per_km2 =
                                mean_intensity_at_cells(result$fit, cells, trend_terms))
  fits[[k]] <- list(fold = k, fit = result$fit, warnings = result$warnings,
                    training_ellipse_id = folds$ellipse_id[training],
                    withheld_ellipse_id = folds$ellipse_id[!training])
  log_lines <- c(log_lines, sprintf("fold %d: %d training centroids, %.2f s, %d warnings, convergence code %s",
                                    k, training_points$n, seconds, length(result$warnings),
                                    parameter_rows[[k]]$mincon_convergence_code))
  if (length(result$warnings)) log_lines <- c(log_lines, paste("   ", unique(result$warnings)))
}

# the model fitted to all 753 centroids
started <- Sys.time()
full_result <- fit_lgcp(inputs$points, inputs$covariates, trend_terms)
seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
parameter_rows[[n_folds + 1]] <- describe_fit("all_data", NA, inputs$points, full_result, seconds)
full_precision_rows[[n_folds + 1]] <- full_precision_table("all_data", NA, full_result$fit)
full_intensity <- mean_intensity_at_cells(full_result$fit, cells, trend_terms)
log_lines <- c(log_lines, sprintf("all data: %d centroids, %.2f s, %d warnings",
                                  inputs$points$n, seconds, length(full_result$warnings)))

parameter_table <- do.call(rbind, parameter_rows)
write_csv_lf(parameter_table, file.path(output_dir, "fold_lgcp_parameters.csv"))
write_csv_lf(do.call(rbind, full_precision_rows),
             file.path(output_dir, "fold_lgcp_parameters_full_precision.csv"))
write_csv_lf(do.call(rbind, surfaces), file.path(output_dir, "fold_mean_intensity_1deg.csv"))
write_csv_lf(data.frame(cell_id = cells$cell_id, lon = cells$lon, lat = cells$lat,
                        mean_intensity_per_km2 = full_intensity,
                        mean_intensity_per_million_km2 = full_intensity * 1e6),
             file.path(output_dir, "full_data_mean_intensity_1deg.csv"))
saveRDS(list(folds = fits, all_data = full_result$fit, trend_terms = trend_terms),
        file.path(output_dir, "fold_lgcp_fits.rds"), compress = "xz")
print(parameter_table[, c("fit", "fold", "intercept", "beta_lon", "beta_lat", "beta_mean_tmax",
                          "beta_std_tmax", "var_sigma2", "scale_km", "mincon_convergence_code",
                          "n_warnings")], digits = 10)

cat(log_lines, sep = "\n")

# check: the all-data fit is the fit of modeling/lgcp
full_model_file <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "lgcp", "lgcp_fit_variant3.rds")
same_as_modeling <- identical(unname(lgcp_parameters(full_result$fit)),
                              unname(readRDS(full_model_file)$parameters))
cat(sprintf("all-data fit identical to modeling/lgcp (bitwise): %s\n", same_as_modeling))
if (!same_as_modeling) quit(status = 1)
