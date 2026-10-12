# Fit the power laws and run the bootstrap calculations.
# Reads: supporting_files/additional_results/modeling/power_law/ and shared study settings.
# Writes: supporting_files/additional_results/modeling/power_law/
# Run: Rscript scripts/modeling/power_law/02_fit_power_law.R

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_dir <- normalizePath(file.path(dirname(script_file), "..", "..", ".."),
                          winslash = "/")
source(file.path(repo_dir, "scripts", "modeling", "power_law", "functions.R"))
config <- power_law_study_settings()

input_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "power_law")
output_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "power_law")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

settings <- config$power_law
ks_mode <- settings$ks_mode
n_bootstrap <- as.integer(settings$n_bootstrap)

samples <- list(
  list(name = "daily_maxima", file = "daily_max_area.csv",
       seed = as.integer(settings$seed_daily_maxima)),
  list(name = "event_maxima", file = "event_max_area.csv",
       seed = as.integer(settings$seed_event_maxima)))

cat(sprintf("%s\nks_mode = %s, %d bootstrap repetitions\n\n",
            R.version.string, ks_mode, n_bootstrap))

fit_rows <- character(0)

for (sample in samples) {
  areas <- read.csv(file.path(input_dir, sample$file))$area_km2

  started <- proc.time()[["elapsed"]]
  fit <- pl_fit(areas, ks_mode = ks_mode)
  stopifnot(identical(fit$status, "ok"))
  boot <- pl_bootstrap(areas, fit, reps = n_bootstrap, ks_mode = ks_mode,
                       seed = sample$seed)
  seconds <- proc.time()[["elapsed"]] - started

  replicate_file <- file.path(output_dir, sprintf("bootstrap_replicates_%s.csv",
                                                  sample$name))
  write_bootstrap_replicates(boot, replicate_file)

  # Summaries of the exponents fitted to the simulated samples (R quantile
  # type 7).
  boot_alpha <- boot$alpha_boot[!is.na(boot$alpha_boot)]
  q <- quantile(boot_alpha, c(0.025, 0.25, 0.5, 0.75, 0.975), type = 7)

  fit_rows <- c(fit_rows, paste(
    sample$name, ks_mode, fit$n, fit$n_tail, format_exact(fit$xmin),
    format_exact(fit$alpha), format_exact(fit$loglik), format_exact(boot$gof),
    format_exact(boot$p), format_exact(boot$p_lower), format_exact(boot$p_upper),
    boot$n_exceed, boot$reps, boot$n_undefined, format_exact(boot$mc_se),
    sample$seed, format_exact(mean(boot_alpha)), format_exact(q[[3]]),
    format_exact(q[[2]]), format_exact(q[[4]]), format_exact(q[[1]]),
    format_exact(q[[5]]), format_exact(min(boot_alpha)),
    format_exact(max(boot_alpha)), sep = ","))

  cat(sprintf(paste0("%s: n = %d, n_tail = %d, A_min = %.17g km2, alpha = %.17g,\n",
                     "  D = %.17g, p = %.4f (%d of %d), undefined refits = %d,\n",
                     "  bootstrap alpha median %.4f, 2.5-97.5%% %.4f-%.4f (%.1f s)\n"),
              sample$name, fit$n, fit$n_tail, fit$xmin, fit$alpha, boot$gof,
              boot$p_lower, boot$n_exceed, boot$reps, boot$n_undefined,
              q[[3]], q[[1]], q[[5]], seconds))
}

fits_file <- file.path(output_dir, "power_law_fits.csv")
con <- file(fits_file, "wb")
writeLines(c(paste0("dataset,ks_mode,n,n_tail,xmin_km2,alpha,loglik,ks_distance,",
                    "p_value,p_lower,p_upper,n_exceed,reps,n_undefined,mc_se,seed,",
                    "alpha_boot_mean,alpha_boot_median,alpha_boot_q25,",
                    "alpha_boot_q75,alpha_boot_q2.5,alpha_boot_q97.5,",
                    "alpha_boot_min,alpha_boot_max"),
             fit_rows), con, sep = "\n")
close(con)

cat(sprintf("\nwrote %s\n", fits_file))
