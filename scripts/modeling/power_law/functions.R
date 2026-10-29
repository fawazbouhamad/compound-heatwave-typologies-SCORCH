# Shared power-law fitting and bootstrap functions.
# Reads: Called by the modeling and supplementary scripts.
# Writes: Returns calculations to the calling script.
# Run: Sourced by the calculation scripts.

# Fixed manuscript settings for the main fit and supplementary comparisons.
power_law_study_settings <- function() {
  list(power_law = list(ks_mode = "full", n_bootstrap = 5000L,
                        seed_daily_maxima = 20261020L,
                        seed_event_maxima = 20261120L),
       supplementary_power_law = list(n_bootstrap = 5000L,
                                       base_seed = 20260918L))
}

pl_ks <- function(tail_sorted, xmin, a, ks_mode) {
  n <- length(tail_sorted)
  fitted_cdf <- 1 - (xmin / tail_sorted)^a
  lower_step <- (0:(n - 1)) / n
  if (identical(ks_mode, "lower")) {
    return(max(abs(fitted_cdf - lower_step)))
  }
  # Full distance. With tied observations the largest gap lies at the ends of
  # each tied run; taking both one-sided gaps at every index covers them.
  upper_step <- (1:n) / n
  max(max(fitted_cdf - lower_step), max(upper_step - fitted_cdf))
}

# ---- fit: try every candidate cutoff, keep the one with the smallest distance -
pl_fit <- function(x, ks_mode = c("lower", "full"), return_scan = FALSE) {
  ks_mode <- match.arg(ks_mode)
  x <- as.numeric(x)
  xs <- sort(x)
  N <- length(xs)

  distinct_values <- unique(xs)                 # increasing, because xs is sorted
  if (length(distinct_values) < 2L) {
    # No candidate is left after dropping the largest value: reported as an
    # undefined fit, never replaced by a fallback value.
    out <- list(xmin = NA_real_, alpha = NA_real_, D = NA_real_,
                n_tail = NA_integer_, loglik = NA_real_, n = N,
                ks_mode = ks_mode, n_candidates = 0L, n_tied_at_min = 0L,
                status = "no_candidates")
    if (return_scan) out$scan <- NULL
    return(out)
  }

  candidates <- distinct_values[-length(distinct_values)]  # drop only the largest
  first_index <- match(candidates, xs)          # first position of each candidate
  n_candidates <- length(candidates)

  distance <- numeric(n_candidates)
  a_scan <- numeric(n_candidates)
  n_scan <- integer(n_candidates)

  for (i in seq_len(n_candidates)) {
    xmin <- candidates[i]
    tail_values <- xs[first_index[i]:N]         # all copies of a tied cutoff kept
    n <- length(tail_values)
    # Written as sum(log(x / xmin)), not sum(log x) - n log(xmin): the two
    # differ in the last digits, which could decide an exact tie.
    a <- n / sum(log(tail_values / xmin))
    n_scan[i] <- n
    a_scan[i] <- a
    distance[i] <- pl_ks(tail_values, xmin, a, ks_mode)
  }

  D <- min(distance)
  selected <- which(distance <= D)[1L]          # smallest cutoff among exact ties
  xmin <- candidates[selected]

  tail_values <- xs[first_index[selected]:N]
  n <- length(tail_values)
  alpha <- 1 + n / sum(log(tail_values / xmin))
  loglik <- n * log((alpha - 1) / xmin) - alpha * sum(log(tail_values / xmin))

  status <- if (is.finite(alpha) && is.finite(D)) "ok" else "nonfinite_fit"

  out <- list(xmin = xmin, alpha = alpha, D = D, n_tail = n, loglik = loglik,
              n = N, ks_mode = ks_mode, n_candidates = n_candidates,
              n_tied_at_min = sum(distance <= D), status = status)
  if (return_scan) {
    out$scan <- data.frame(candidate_index = seq_len(n_candidates),
                           xmin = candidates, n_tail = n_scan,
                           alpha = 1 + a_scan, D = distance)
  }
  out
}

# goodness-of-fit bootstrap (the semiparametric scheme of plpva.m)
# Each repetition builds a simulated sample of the same size N:
#   - each of the N observations is independently assigned to the body (below
#     the cutoff) or to the tail, with the observed tail fraction;
#   - body values are drawn with replacement from the observed values below the
#     cutoff; tail values are drawn from the fitted power law;
#   - cutoff and exponent are then fitted again to the simulated sample, with
#     the same distance definition, and its distance is recorded.
# The random draws happen in this order in every repetition: N uniforms (body
# or tail), then one uniform per body value, then one per tail value.
pl_bootstrap <- function(x, fit, reps, ks_mode = c("lower", "full"),
                         seed = NULL, progress_every = 0L) {
  ks_mode <- match.arg(ks_mode)
  stopifnot(identical(fit$ks_mode, ks_mode))
  x <- as.numeric(x)
  N <- length(x)

  xmin <- fit$xmin
  tail_values <- x[x >= xmin]; n_tail <- length(tail_values)
  body_values <- x[x < xmin];  n_body <- length(body_values)  # in file order
  alpha <- 1 + n_tail / sum(log(tail_values / xmin))

  # The observed distance is recomputed at the given cutoff, as plpva.m does.
  observed_distance <- pl_ks(sort(tail_values), xmin, alpha - 1, ks_mode)
  tail_fraction <- n_tail / N

  if (!is.null(seed)) set.seed(seed, kind = "Mersenne-Twister")

  simulated_distance <- rep(NA_real_, reps)
  status <- rep(NA_character_, reps)
  alpha_boot <- rep(NA_real_, reps)       # recorded only; draws no random numbers
  xmin_boot <- rep(NA_real_, reps)
  ntail_boot <- rep(NA_integer_, reps)

  for (b in seq_len(reps)) {
    # MATLAB: n1 = sum(rand(N,1) > pz). One decision per observation (not a
    # single binomial draw); the tail count follows by subtraction.
    n1 <- sum(runif(N) > tail_fraction)
    n2 <- N - n1

    q1 <- if (n1 > 0L && n_body > 0L) body_values[ceiling(n_body * runif(n1))] else numeric(0)
    q2 <- if (n2 > 0L) xmin * (1 - runif(n2))^(-1 / (alpha - 1)) else numeric(0)
    simulated <- c(q1, q2)

    refit <- pl_fit(simulated, ks_mode = ks_mode)

    # Keep undefined fits as NA.
    alpha_boot[b] <- refit$alpha
    xmin_boot[b] <- refit$xmin
    ntail_boot[b] <- refit$n_tail

    if (identical(refit$status, "ok")) {
      simulated_distance[b] <- refit$D
      status[b] <- "ok"
    } else {
      # Count undefined fits in the total.
      status[b] <- refit$status
    }

    if (progress_every > 0L && b %% progress_every == 0L) {
      cat(sprintf("      rep %d/%d\n", b, reps)); flush.console()
    }
  }

  ok <- !is.na(simulated_distance)
  n_undefined <- sum(!ok)
  n_exceed <- sum(simulated_distance[ok] >= observed_distance)

  # The denominator is the number of repetitions requested. If any refit was
  # undefined, the ratio is only a lower bound (undefined outcomes counted as
  # non-exceedances); the upper bound counts them as exceedances, and the
  # p-value itself is left undefined. These bounds are not a confidence
  # interval.
  p_lower <- n_exceed / reps
  p_upper <- (n_exceed + n_undefined) / reps
  complete <- (n_undefined == 0L)
  p <- if (complete) p_lower else NA_real_

  # Monte Carlo standard error of a completed p-value.
  mc_se <- if (complete) sqrt(p_lower * (1 - p_lower) / reps) else NA_real_

  list(p = p, p_lower = p_lower, p_upper = p_upper, complete = complete,
       n_exceed = n_exceed, gof = observed_distance, reps = reps,
       n_failed = n_undefined, n_undefined = n_undefined,
       n_valid = sum(ok), mc_se = mc_se, nof = simulated_distance,
       status = status, alpha_boot = alpha_boot, xmin_boot = xmin_boot,
       ntail_boot = ntail_boot, alpha_at_xmin = alpha, ks_mode = ks_mode)
}

# writing numbers
# 17 significant digits identify a double exactly; missing values become "NA".
format_exact <- function(values) {
  ifelse(is.na(values), "NA", sprintf("%.17g", values))
}

# Per-replicate table of a bootstrap run: one row per simulated sample.
write_bootstrap_replicates <- function(boot, path) {
  lines <- c("replicate,alpha,xmin,n_tail,distance,status,exceeds_observed",
             sprintf("%d,%s,%s,%s,%s,%s,%s",
                     seq_len(boot$reps),
                     format_exact(boot$alpha_boot),
                     format_exact(boot$xmin_boot),
                     ifelse(is.na(boot$ntail_boot), "NA",
                            as.character(boot$ntail_boot)),
                     format_exact(boot$nof),
                     boot$status,
                     ifelse(is.na(boot$nof), "NA",
                            ifelse(boot$nof >= boot$gof, "1", "0"))))
  con <- file(path, "wb")
  writeLines(lines, con, sep = "\n")
  close(con)
}
