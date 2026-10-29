# Fit the LGCP and predict spatial centroid intensity.
# Reads: supporting_files/additional_results/modeling/lgcp/ and shared study settings.
# Writes: supporting_files/additional_results/modeling/lgcp/
# Run: Rscript scripts/modeling/lgcp/02_fit_lgcp.R

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script_file), "functions.R"))

repo_dir <- find_repo_dir()
config <- lgcp_study_settings()
output_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "lgcp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trend_terms <- unlist(config$lgcp$trend_terms)
# The minimum-contrast fit is deterministic. Set the configured seed for
# reproducibility if any package operation uses random numbers.
set.seed(config$lgcp$seed)

inputs <- load_lgcp_inputs(repo_dir)
cat(sprintf("centroids = %d, grid cells = %d, raster pixels inside = %d, window area = %.1f km^2\n",
            inputs$points$n, nrow(inputs$cells), sum(inputs$raster$inside == 1),
            area(inputs$window)))

# fit
started <- Sys.time()
result <- fit_lgcp(inputs$points, inputs$covariates, trend_terms)
fit_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
fit <- result$fit
parameters <- lgcp_parameters(fit)
diagnostics <- fit_diagnostics(fit)
cat(sprintf("fit time %.1f s, warnings: %d\n", fit_seconds, length(result$warnings)))
print(parameters, digits = 15)

parameter_table <- data.frame(
  model = "variant3",
  formula = paste("~", paste(trend_terms, collapse = " + ")),
  parameter = names(parameters),
  value = sprintf("%.17g", parameters),
  value_6_decimals = sprintf("%.6f", parameters),
  stringsAsFactors = FALSE)
write_csv_lf(parameter_table, file.path(output_dir, "lgcp_parameters.csv"))
saveRDS(list(fit = fit, parameters = parameters, warnings = result$warnings,
             trend_terms = trend_terms),
        file.path(output_dir, "lgcp_fit_variant3.rds"), compress = "xz")

# predictions
cells <- inputs$cells
intensity <- mean_intensity_at_cells(fit, cells, trend_terms)
predictions <- data.frame(
  cell_id = cells$cell_id, lon = cells$lon, lat = cells$lat,
  mean_tmax = cells$mean_tmax, std_tmax = cells$std_tmax,
  mean_intensity_per_km2 = intensity,
  mean_intensity_per_million_km2 = intensity * 1e6,
  concentration_rank = concentration_rank(intensity * 1e6))
# write.csv keeps 15 significant digits.
write_csv_lf(predictions, file.path(output_dir, "predicted_intensity_1deg.csv"))

intensity_image <- mean_intensity_image(fit, inputs$covariates, trend_terms)
raster <- inputs$raster
pixel_intensity <- as.vector(t(intensity_image$v))   # rows south to north, west to east
stopifnot(identical(is.na(pixel_intensity), raster$inside == 0))
inside <- raster$inside == 1
write_csv_lf(data.frame(ix = raster$ix[inside], iy = raster$iy[inside],
                     x_km = raster$x_km[inside], y_km = raster$y_km[inside],
                     pixel_lon = raster$pixel_lon[inside], pixel_lat = raster$pixel_lat[inside],
                     mean_intensity_per_km2 = pixel_intensity[inside],
                     mean_intensity_per_million_km2 = pixel_intensity[inside] * 1e6),
          file.path(output_dir, "predicted_intensity_25km.csv"))

# log
log_lines <- c(
  sprintf("R: %s", R.version.string),
  sprintf("spatstat.geom %s, spatstat.model %s, spatstat.explore %s",
          packageVersion("spatstat.geom"), packageVersion("spatstat.model"),
          packageVersion("spatstat.explore")),
  sprintf("fit time: %.2f s", fit_seconds),
  sprintf("warnings during the fit: %d", length(result$warnings)),
  if (length(result$warnings)) paste("  ", unique(result$warnings)) else NULL,
  sprintf("minimum-contrast convergence code: %s (0 = converged)",
          diagnostics$mincon_convergence_code),
  sprintf("minimum-contrast objective: %.17g after %s evaluations",
          diagnostics$mincon_objective, diagnostics$mincon_n_objfun_evals),
  sprintf("trend fit converged: %s", diagnostics$trend_glm_converged))
cat(log_lines, sep = "\n")
