#!/usr/bin/env Rscript
# Draw maps from the historical LGCP covariate comparison.
#
# Uses saved coefficients from spatial_risk_lgcp_model_diagnostics.R and the
# covariate raster; models are not refitted. Images are written under
# model_diagnostics_and_variants/figures_final_for_review/.
# The default input path belongs to the original workspace. The first command
# argument can supply that input folder.

suppressMessages({
  library(spatstat.geom)
})

args <- commandArgs(trailingOnly = TRUE)
IN <- if (length(args) >= 1) args[1] else
  file.path("outputs", "event_global_max_algorithm",
            "spatial_risk_lgcp_tmax_centroids")
DIAG <- file.path(IN, "model_diagnostics_and_variants")
OUTF <- file.path(DIAG, "figures_final_for_review")
dir.create(OUTF, recursive = TRUE, showWarnings = FALSE)

ip  <- function(...) file.path(IN, ...)
out <- function(...) file.path(OUTF, ...)
say <- function(...) cat(sprintf(...), "\n")

# Read inputs and saved coefficients.
boxes  <- read.csv(ip("tmax_covariate_grid_all_boxes.csv"))
rast   <- read.csv(ip("covariate_raster_km.csv"))
cen    <- read.csv(ip("centroid_points_all_ellipses_validated.csv"))
meta   <- jsonlite::fromJSON(ip("validation_summary.json"))
params <- read.csv(file.path(DIAG, "tables", "lgcp_model_parameters.csv"),
                   stringsAsFactors = FALSE)
say("[load] raster pixels=%d  centroids=%d  param rows=%d",
    nrow(rast), nrow(cen), nrow(params))

# Map std_tmax and its standardized values from the grid to the raster.
pop_sd <- function(x) sqrt(mean((x - mean(x))^2))
mu_s <- mean(boxes$std_tmax); sd_s <- pop_sd(boxes$std_tmax)
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3), sep = "_")
std_lookup <- setNames(boxes$std_tmax, key(boxes$lon, boxes$lat))
rast$std_tmax   <- as.numeric(std_lookup[key(rast$lon, rast$lat)])
rast$std_tmax_z <- (rast$std_tmax - mu_s) / sd_s

# Build spatstat images from the covariate raster.
nx <- meta$raster$nx; ny <- meta$raster$ny
xs <- sort(unique(rast$x_km)); ys <- sort(unique(rast$y_km))
as_mat <- function(v) matrix(v, nrow = ny, ncol = nx, byrow = TRUE)
inside_mat <- as_mat(rast$inside == 1)
mk_im <- function(v) { m <- as_mat(v); m[!inside_mat] <- NA; im(m, xcol = xs, yrow = ys) }

lon_im <- mk_im(rast$pixel_lon); lat_im <- mk_im(rast$pixel_lat)
mtz_im <- mk_im(rast$mean_tmax_z); cvz_im <- mk_im(rast$cv_tmax_z)
stz_im <- mk_im(rast$std_tmax_z)
mt_im  <- mk_im(rast$mean_tmax); cv_im <- mk_im(rast$cv_tmax); st_im <- mk_im(rast$std_tmax)
covim  <- list(lon = lon_im, lat = lat_im, mean_tmax_z = mtz_im,
               cv_tmax_z = cvz_im, std_tmax_z = stz_im)

# Calculate intensity as exp(trend) from saved coefficients; scale by 1e6.
coef_of <- function(model) {
  d <- params[params$model == model & params$param_class == "fixed_effect", ]
  setNames(as.numeric(d$value), d$parameter)
}
intensity_image <- function(cf) {
  terms <- setdiff(names(cf), "(Intercept)")
  acc <- lon_im * 0 + cf[["(Intercept)"]]
  for (t in terms) acc <- acc + cf[[t]] * covim[[t]]
  exp(acc) * 1e6
}

MODELS <- list(
  list(key = "baseline", file = "risk_map_baseline_mean_cv.png",
       label = "Baseline: lon + lat + mean Tmax + Tmax CV"),
  list(key = "variant1", file = "risk_map_lon_lat_mean.png",
       label = "Mean Tmax model: lon + lat + mean Tmax"),
  list(key = "variant2", file = "risk_map_lon_lat_cv.png",
       label = "Tmax CV model: lon + lat + Tmax CV"),
  list(key = "variant3", file = "risk_map_lon_lat_mean_std.png",
       label = "Mean + std model: lon + lat + mean Tmax + Tmax std"))

lam_ims <- lapply(MODELS, function(m) intensity_image(coef_of(m$key)))
names(lam_ims) <- vapply(MODELS, function(m) m$key, "")

# Use one colour scale for all models, capped at the pooled 99th percentile.
pooled <- unlist(lapply(lam_ims, function(im) as.numeric(im$v)))
pooled <- pooled[is.finite(pooled)]
ZLO <- 0
ZHI <- as.numeric(quantile(pooled, 0.99))
say("[scale] shared intensity color range = [%.1f, %.1f] (pooled 99th pctile)", ZLO, ZHI)

# Values above the shared limit use the highest colour.
lam_plot <- lapply(lam_ims, function(im) eval.im(pmin(im, ZHI)))

RISK_COL <- hcl.colors(256, "YlOrRd", rev = TRUE)
COV_COL  <- hcl.colors(256, "viridis")
INT_TITLE <- expression("Fitted centroid intensity" %*% 10^6)

# Map layout helpers
add_axes <- function() {
  axis(1, cex.axis = 0.85); axis(2, cex.axis = 0.85); box()
  title(xlab = "Easting (km)", ylab = "Northing (km)",
        line = 2.2, cex.lab = 0.9)
}

# Draw a map with optional units, subtitle, centroids and colour limits.
clean_map <- function(image, main, file, col, bar = NULL, sub = NULL,
                      overlay = FALSE, zlim = NULL) {
  png(out(file), width = 7.6, height = 5.9, units = "in", res = 300)
  op <- par(mar = c(3.8, 3.8, if (is.null(sub)) 2.6 else 3.4, 5.4))
  if (is.null(zlim))
    plot(image, main = "", col = col, ribbon = TRUE, box = FALSE)
  else
    plot(image, main = "", col = col, ribbon = TRUE, box = FALSE, zlim = zlim)
  title(main = main, cex.main = 1.1, line = if (is.null(sub)) 1.1 else 1.9)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.4, cex = 0.85)
  if (overlay) points(cen$x_km, cen$y_km, pch = 16, cex = 0.28, col = "black")
  add_axes()
  if (!is.null(bar)) mtext(bar, side = 4, line = 3.6, cex = 0.95, las = 0)
  par(op); dev.off()
  say("[fig] %s", out(file))
}

# Covariate maps
say("Covariate maps")
clean_map(mt_im, "Historical warm-season mean Tmax", "map_tmax_mean.png",
          COV_COL, bar = expression(paste("Mean Tmax (", degree, "C)")))
clean_map(cv_im, "Historical warm-season Tmax coefficient of variation",
          "map_tmax_cv.png", COV_COL, bar = "Tmax CV")
clean_map(st_im, "Historical warm-season Tmax standard deviation",
          "map_tmax_std.png", COV_COL,
          bar = expression(paste("Tmax std (", degree, "C)")))

# LGCP intensity maps with a shared scale and centroid locations.
say("LGCP fitted-intensity maps (shared color scale)")
for (m in MODELS)
  clean_map(lam_plot[[m$key]], INT_TITLE, m$file, RISK_COL,
            bar = expression("intensity" %*% 10^6),
            sub = m$label, overlay = TRUE, zlim = c(ZLO, ZHI))

# Three-panel covariate figure
say("Combined panels")
png(out("panel_covariate_maps.png"), width = 15.0, height = 5.7, units = "in", res = 300)
op <- par(mfrow = c(1, 3), mar = c(3.6, 3.6, 2.6, 5.0), oma = c(0, 0, 2.8, 0))
cov_specs <- list(
  list(im = mt_im, t = "Mean Tmax", b = expression(paste("(", degree, "C)"))),
  list(im = cv_im, t = "Tmax coefficient of variation", b = "CV"),
  list(im = st_im, t = "Tmax standard deviation", b = expression(paste("(", degree, "C)"))))
for (s in cov_specs) {
  plot(s$im, main = "", col = COV_COL, ribbon = TRUE, box = FALSE)
  title(main = s$t, cex.main = 1.1, line = 0.6)
  add_axes(); mtext(s$b, side = 4, line = 3.0, cex = 0.9, las = 0)
}
mtext("Historical warm-season Tmax covariates over 1800 ERA5 boxes",
      side = 3, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
par(op); dev.off()
say("[fig] %s", out("panel_covariate_maps.png"))

# Four-panel model comparison with a shared scale and centroids.
png(out("panel_model_comparison.png"), width = 11.0, height = 9.8, units = "in", res = 300)
op <- par(mfrow = c(2, 2), mar = c(3.4, 3.6, 3.0, 4.6), oma = c(0, 0, 2.6, 0))
for (m in MODELS) {
  plot(lam_plot[[m$key]], main = "", col = RISK_COL, ribbon = TRUE,
       box = FALSE, zlim = c(ZLO, ZHI))
  title(main = m$label, cex.main = 1.0, line = 1.0)
  points(cen$x_km, cen$y_km, pch = 16, cex = 0.26, col = "black")
  add_axes()
}
mtext(INT_TITLE, side = 3, outer = TRUE, line = 0.4, cex = 1.25, font = 2)
par(op); dev.off()
say("[fig] %s", out("panel_model_comparison.png"))

# List saved figures
cat(strrep("-", 70), "\n")
say("Final-review figures written to: %s", OUTF)
files <- sort(list.files(OUTF, pattern = "\\.png$", full.names = TRUE))
for (f in files) say("  %s", f)
say("DONE — %d figures.", length(files))
