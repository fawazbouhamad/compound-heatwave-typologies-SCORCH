# Run the R power-law comparison methods for Figures S2 and S3.
# Reads: Prepared maximum-area tables and shared study settings.
# Writes: supporting_files/supporting_information/power_law_comparison/results/
# Run: Rscript supporting_files/supporting_information/power_law_comparison/01_compare_r.R

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_dir <- normalizePath(file.path(dirname(script_file), "..", "..", ".."),
                          winslash = "/")
source(file.path(repo_dir, "scripts", "modeling", "power_law", "functions.R"))
config <- power_law_study_settings()

# bootstrap_p runs its simulations in separate R worker processes (even with
# threads = 1); pass the current library search path to those workers.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
suppressPackageStartupMessages(library(poweRlaw))
stopifnot(as.character(packageVersion("poweRlaw")) == "1.0.0")

base_seed <- as.integer(config$supplementary_power_law$base_seed)
n_bootstrap <- as.integer(config$supplementary_power_law$n_bootstrap)

input_dir <- file.path(repo_dir, "supporting_files", "additional_results", "modeling", "power_law")
output_dir <- file.path(repo_dir, "supporting_files", "supporting_information", "power_law_comparison", "results")
results_dir <- file.path(output_dir, "fits")
replicates_dir <- file.path(output_dir, "replicates")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(replicates_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  daily_maxima = list(file = "daily_max_area.csv", index_m1_m2 = 1L, index_m3_m4 = 1L),
  event_maxima = list(file = "event_max_area.csv", index_m1_m2 = 2L, index_m3_m4 = 6L))

configurations <- list(
  M1 = list(family = "clauset", ks_mode = "lower", number = 1L,
            label = "M1 reference-algorithm R (lower-rank KS)",
            distance_definition = "lower-rank KS"),
  M2 = list(family = "clauset", ks_mode = "full", number = 2L,
            label = "M2 full-KS variant (base R)",
            distance_definition = "full KS"),
  M3_default = list(family = "poweRlaw", id = "strict-default", xmax = 1e5,
                    number = 1L,
                    label = "M3 poweRlaw 1.0.0 (strict defaults; xmax=1e5)"),
  M3_adjusted = list(family = "poweRlaw", id = "adjusted-xmax", xmax = Inf,
                     number = 2L,
                     label = "M3 poweRlaw 1.0.0 (ADJUSTED: xmax=Inf)"))

arguments <- commandArgs(trailingOnly = TRUE)
selected_configurations <- intersect(arguments, names(configurations))
selected_datasets <- intersect(arguments, names(datasets))
if (length(selected_configurations) == 0) selected_configurations <- names(configurations)
if (length(selected_datasets) == 0) selected_datasets <- names(datasets)

write_lines_lf <- function(lines, path) {
  con <- file(path, "wb"); writeLines(lines, con, sep = "\n"); close(con)
}
number_or_na <- function(v) {
  if (is.null(v) || length(v) == 0 || !is.finite(v)) "NA" else sprintf("%.17g", v)
}

# M1 and M2: the Clauset procedure (functions.R)
run_clauset <- function(areas, dataset_name, configuration, seed) {
  started <- proc.time()[["elapsed"]]
  fit <- pl_fit(areas, ks_mode = configuration$ks_mode)
  header <- paste0("dataset,method,ks_mode,n,n_tail,xmin,alpha,loglik,",
                   "distance_definition,distance,p_value,p_lower,p_upper,",
                   "bootstrap_complete,n_exceed,reps,n_failed,n_undefined,",
                   "n_valid,mc_se,seed,seconds,r_version,status")
  if (!identical(fit$status, "ok")) {
    # An undefined fit is reported as such; nothing is substituted.
    seconds <- proc.time()[["elapsed"]] - started
    return(c(header, sprintf(
      "%s,%s,%s,%d,NA,NA,NA,NA,%s,NA,NA,NA,NA,NA,NA,%d,NA,NA,NA,NA,%d,%.3f,%s,%s",
      dataset_name, configuration$label, configuration$ks_mode, length(areas),
      configuration$distance_definition, n_bootstrap, seed, seconds,
      R.version.string, fit$status)))
  }
  boot <- pl_bootstrap(areas, fit, reps = n_bootstrap,
                       ks_mode = configuration$ks_mode, seed = seed)
  seconds <- proc.time()[["elapsed"]] - started
  write_bootstrap_replicates(boot, file.path(replicates_dir,
    sprintf("%s_%s.csv", dataset_name, configuration$name)))
  c(header, sprintf(
    "%s,%s,%s,%d,%d,%.17g,%.17g,%.17g,%s,%.17g,%s,%s,%s,%s,%d,%d,%d,%d,%d,%s,%d,%.3f,%s,%s",
    dataset_name, configuration$label, configuration$ks_mode, fit$n, fit$n_tail,
    fit$xmin, fit$alpha, fit$loglik, configuration$distance_definition, boot$gof,
    number_or_na(boot$p), number_or_na(boot$p_lower), number_or_na(boot$p_upper),
    if (boot$complete) "TRUE" else "FALSE", boot$n_exceed, boot$reps,
    boot$n_failed, boot$n_undefined, boot$n_valid, number_or_na(boot$mc_se),
    seed, seconds, R.version.string,
    if (boot$complete) "ok" else "incomplete_bootstrap"))
}

# M3: the poweRlaw package
run_powerRlaw <- function(areas, dataset_name, configuration, seed) {
  started <- proc.time()[["elapsed"]]
  status <- "ok"; note <- ""
  n_tail <- NA_integer_; xmin <- NA_real_; alpha <- NA_real_
  distance <- NA_real_; p_value <- NA_real_; n_failed <- NA_integer_
  n_valid <- NA_integer_; mc_se <- NA_real_

  tryCatch({
    model <- conpl$new(areas)                     # continuous power law
    estimate <- estimate_xmin(model, xmax = configuration$xmax)
    if (is.null(estimate$xmin) || !is.finite(estimate$xmin)) {
      status <- "no_usable_fit"
      note <- paste0("estimate_xmin returned no finite cutoff with xmax=",
                     format(configuration$xmax))
    } else {
      model$setXmin(estimate)
      xmin <- estimate$xmin
      alpha <- estimate$pars
      distance <- estimate$gof                    # the package's own KS distance
      n_tail <- sum(areas >= estimate$xmin)
      boot <- bootstrap_p(model, xmax = configuration$xmax,
                          no_of_sims = n_bootstrap, threads = 1, seed = seed)
      p_value <- boot$p
      simulated <- boot$bootstraps
      # poweRlaw leaves failed simulations out of the numerator but divides by
      # the requested simulation count. Record failed simulations separately.
      n_valid <- sum(!is.na(simulated$gof))
      n_failed <- n_bootstrap - n_valid
      mc_se <- sqrt(p_value * (1 - p_value) / n_bootstrap)
      write_lines_lf(c(
        "replicate,xmin,alpha,distance,exceeds_observed",
        sprintf("%d,%s,%s,%s,%s", seq_len(nrow(simulated)),
                sprintf("%.17g", simulated$xmin), sprintf("%.17g", simulated$pars),
                sprintf("%.17g", simulated$gof),
                ifelse(is.na(simulated$gof), "NA",
                       ifelse(simulated$gof >= distance, "1", "0")))),
        file.path(replicates_dir, sprintf("%s_%s.csv", dataset_name,
                                          configuration$name)))
    }
  }, error = function(e) {
    status <<- "error"
    note <<- gsub("[,\n]", " ", conditionMessage(e))
  })

  seconds <- proc.time()[["elapsed"]] - started
  c(paste0("dataset,method,config,n,n_tail,xmin,alpha,distance_definition,",
           "distance,p_value,reps,n_failed,n_valid,mc_se,seed,seconds,",
           "package_version,status,note"),
    sprintf("%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%d,%.3f,%s,%s,%s",
            dataset_name, configuration$label, configuration$id, length(areas),
            if (is.na(n_tail)) "NA" else as.character(n_tail),
            number_or_na(xmin), number_or_na(alpha), "poweRlaw native KS",
            number_or_na(distance), number_or_na(p_value), n_bootstrap,
            if (is.na(n_failed)) "NA" else as.character(n_failed),
            if (is.na(n_valid)) "NA" else as.character(n_valid),
            number_or_na(mc_se), seed, seconds,
            as.character(packageVersion("poweRlaw")), status, note))
}

# run
cat(sprintf("%s; poweRlaw %s; %d bootstrap repetitions\n", R.version.string,
            as.character(packageVersion("poweRlaw")), n_bootstrap))

for (dataset_name in selected_datasets) {
  dataset <- datasets[[dataset_name]]
  areas <- read.csv(file.path(input_dir, dataset$file))$area_km2   # file order

  for (configuration_name in selected_configurations) {
    configuration <- configurations[[configuration_name]]
    configuration$name <- configuration_name
    if (configuration$family == "clauset") {
      seed <- base_seed + 100L * dataset$index_m1_m2 + configuration$number
      row <- run_clauset(areas, dataset_name, configuration, seed)
    } else {
      seed <- base_seed + 100L * dataset$index_m3_m4 + 50L + configuration$number
      row <- run_powerRlaw(areas, dataset_name, configuration, seed)
    }
    write_lines_lf(row, file.path(results_dir, sprintf("%s_%s.csv", dataset_name,
                                                       configuration_name)))
    values <- read.csv(text = row, colClasses = "character")
    cat(sprintf("  %-13s %-12s status=%s n_tail=%s xmin=%s alpha=%s p=%s seed=%d (%s s)\n",
                dataset_name, configuration_name, values[["status"]],
                values[["n_tail"]], values[["xmin"]], values[["alpha"]],
                values[["p_value"]], seed, values[["seconds"]]))
    flush.console()
  }
}
