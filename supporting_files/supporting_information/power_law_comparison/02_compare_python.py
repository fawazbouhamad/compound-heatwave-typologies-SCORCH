"""Run the Python powerlaw comparison methods.

Reads: Prepared maximum-area tables; study settings are recorded below.
Writes: supporting_files/supporting_information/power_law_comparison/results/
Run: python supporting_files/supporting_information/power_law_comparison/02_compare_python.py
"""
from __future__ import annotations

import csv
import sys
import time
import warnings
from pathlib import Path

import numpy as np

REPO_DIR = Path(__file__).resolve().parents[3]
INPUT_DIR = REPO_DIR / "supporting_files/additional_results/modeling/power_law"
OUTPUT_DIR = REPO_DIR / "supporting_files/supporting_information/power_law_comparison/results"

DATASETS = {"daily_maxima": dict(file="daily_max_area.csv", index=1),
            "event_maxima": dict(file="event_max_area.csv", index=6)}
CONFIGURATIONS = {
    "M4_default": dict(
        id="strict-default", number=1, kwargs={},
        label="M4 powerlaw 2.0.0 (strict defaults; alpha range [0-3])"),
    "M4_adjusted": dict(
        id="adjusted-alpha", number=2,
        kwargs={"parameter_ranges": {"alpha": [1, None]}},
        label="M4 powerlaw 2.0.0 (ADJUSTED: alpha upper bound removed)"),
}

HEADER = ("dataset,method,config,n,n_tail,xmin,alpha,distance_definition,"
          "distance,p_value,reps,n_failed,n_valid,n_flagged,mc_se,seed,seconds,"
          "obs_fit_noise_flag,obs_pl_in_range,obs_pl_noise_flag,"
          "package_version,status,note\n")


def read_supplementary_settings():
    """Return the recorded base seed and bootstrap count, in that order."""
    return 20260918, 5000


def read_areas(path):
    """area_km2 column in file order, with Python's exact number parser."""
    with open(path, newline="", encoding="utf-8") as fh:
        return np.array([float(row["area_km2"]) for row in csv.DictReader(fh)])


def fmt(v):
    if v is None:
        return "NA"
    try:
        x = float(v)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if not np.isfinite(x) else repr(x)


def fit_once(data, kwargs):
    """One powerlaw.Fit under the given configuration, warnings suppressed."""
    import powerlaw
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        old = np.seterr(all="ignore")
        try:
            return powerlaw.Fit(data, verbose=False, **kwargs)
        finally:
            np.seterr(**old)


def fit_flags(fit):
    """Validity state of one fit.

    fit_noise_flag: the package found no candidate cutoff that passed its own
    validity test and fell back to the smallest distance overall.
    pl_in_range / pl_noise_flag: whether the selected fit itself lies inside the
    parameter range / was flagged by the package.
    """
    out = dict(fit_noise_flag=None, pl_in_range=None, pl_noise_flag=None)
    try:
        out["fit_noise_flag"] = bool(fit.noise_flag)
    except Exception:  # noqa: BLE001
        pass
    try:
        out["pl_in_range"] = bool(fit.power_law.in_range())
    except Exception:  # noqa: BLE001
        pass
    try:
        out["pl_noise_flag"] = bool(fit.power_law.noise_flag)
    except Exception:  # noqa: BLE001
        pass
    return out


def is_flagged(flags):
    return bool(flags.get("fit_noise_flag") or flags.get("pl_noise_flag")
                or (flags.get("pl_in_range") is False))


def run_configuration(dataset_name, configuration_name, base_seed, reps):
    import powerlaw

    dataset = DATASETS[dataset_name]
    configuration = CONFIGURATIONS[configuration_name]
    x = read_areas(INPUT_DIR / dataset["file"])
    N = x.size
    seed = base_seed + dataset["index"] * 100 + 70 + configuration["number"]

    t0 = time.time()
    status, note, lr_note = "ok", "", ""
    xmin = alpha = dist = pval = mc_se = None
    n_tail = n_failed = n_valid = n_flagged = None
    obs_flags = dict(fit_noise_flag=None, pl_in_range=None, pl_noise_flag=None)

    try:
        fit = fit_once(x, configuration["kwargs"])
        xmin = float(fit.xmin)
        alpha = float(fit.alpha)
        dist = float(fit.D)                      # the package's own distance
        n_tail = int(fit.power_law.n)
        obs_flags = fit_flags(fit)

        # Likelihood-ratio comparisons: a different test, recorded for context.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r_ln, p_ln = fit.distribution_compare("power_law", "lognormal",
                                                      normalized_ratio=True)
                r_ex, p_ex = fit.distribution_compare("power_law", "exponential",
                                                      normalized_ratio=True)
            lr_note = ("LR vs lognormal R=%.4f p=%.4f; vs exponential "
                       "R=%.4f p=%.4f" % (r_ln, p_ln, r_ex, p_ex))
        except Exception as exc:  # noqa: BLE001
            lr_note = "LR comparison unavailable: %s" % exc

        # Goodness-of-fit bootstrap for the fitted powerlaw model.
        below = x[x < xmin]                      # file order
        n_below = below.size
        p_tail = float(fit.power_law.n) / N
        rng = np.random.default_rng(seed)

        sim_d = np.full(reps, np.nan)
        sim_xmin = np.full(reps, np.nan)
        sim_alpha = np.full(reps, np.nan)
        statuses = []
        sim_flags = []

        for b in range(reps):
            n1 = int(np.sum(rng.random(N) > p_tail))   # body count
            n2 = N - n1                                # tail count by subtraction
            q1 = (rng.choice(below, size=n1, replace=True)
                  if (n1 and n_below) else np.empty(0))
            q2 = (xmin * (1.0 - rng.random(n2)) ** (-1.0 / (alpha - 1.0))
                  if n2 else np.empty(0))
            q = np.concatenate([q1, q2])
            try:
                refit = fit_once(q, configuration["kwargs"])
                d = float(refit.D)
                if not np.isfinite(d):
                    raise ValueError("non-finite distance")
                sim_d[b] = d
                sim_xmin[b] = float(refit.xmin)
                sim_alpha[b] = float(refit.alpha)
                flags = fit_flags(refit)
                sim_flags.append(flags)
                # Fallback fits are kept, but labelled.
                statuses.append("flagged" if is_flagged(flags) else "ok")
            except Exception as exc:  # noqa: BLE001
                # Count failed fits in the denominator.
                sim_flags.append(dict(fit_noise_flag=None, pl_in_range=None,
                                      pl_noise_flag=None))
                statuses.append(type(exc).__name__)

            if reps >= 500 and (b + 1) % 500 == 0:
                print("      %s/%s rep %d/%d (%.1f min)"
                      % (dataset_name, configuration_name, b + 1, reps,
                         (time.time() - t0) / 60), flush=True)

        ok = ~np.isnan(sim_d)
        n_valid = int(ok.sum())
        n_failed = int(reps - n_valid)
        n_flagged = int(sum(1 for s in statuses if s == "flagged"))
        # Denominator: the number of repetitions requested.
        pval = float(np.sum(sim_d[ok] >= dist) / reps)
        mc_se = float(np.sqrt(pval * (1 - pval) / reps))

        replicate_path = OUTPUT_DIR / "replicates" / (
            "%s_%s.csv" % (dataset_name, configuration_name))
        with open(replicate_path, "w", newline="\n", encoding="utf-8") as fh:
            fh.write("replicate,xmin,alpha,distance,status,"
                     "fit_noise_flag,pl_in_range,pl_noise_flag,"
                     "exceeds_observed\n")
            for i in range(reps):
                exceeds = ("NA" if np.isnan(sim_d[i])
                           else ("1" if sim_d[i] >= dist else "0"))
                flags = sim_flags[i] if i < len(sim_flags) else {}
                fh.write("%d,%s,%s,%s,%s,%s,%s,%s,%s\n"
                         % (i + 1, fmt(sim_xmin[i]), fmt(sim_alpha[i]),
                            fmt(sim_d[i]), statuses[i],
                            flags.get("fit_noise_flag"), flags.get("pl_in_range"),
                            flags.get("pl_noise_flag"), exceeds))
    except Exception as exc:  # noqa: BLE001
        status = "error"
        note = str(exc).replace(",", " ").replace("\n", " ")

    elapsed = time.time() - t0
    note_full = (note + (" | " if note and lr_note else "") + lr_note).strip()
    result_path = OUTPUT_DIR / "fits" / (
        "%s_%s.csv" % (dataset_name, configuration_name))
    with open(result_path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(HEADER)
        fh.write("%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%d,%.3f,"
                 "%s,%s,%s,%s,%s,%s\n"
                 % (dataset_name, configuration["label"], configuration["id"], N,
                    "NA" if n_tail is None else n_tail,
                    fmt(xmin), fmt(alpha),
                    "powerlaw native KS (D)", fmt(dist), fmt(pval), reps,
                    "NA" if n_failed is None else n_failed,
                    "NA" if n_valid is None else n_valid,
                    "NA" if n_flagged is None else n_flagged,
                    fmt(mc_se), seed, elapsed,
                    obs_flags.get("fit_noise_flag"),
                    obs_flags.get("pl_in_range"),
                    obs_flags.get("pl_noise_flag"),
                    powerlaw.__version__, status, note_full))

    if status == "ok":
        print("  %-13s %-12s xmin=%.4g alpha=%.4f ntail=%d D=%.6f p=%.4f "
              "failed=%d flagged=%d seed=%d (%.1f s)"
              % (dataset_name, configuration_name, xmin, alpha, n_tail, dist,
                 pval, n_failed, n_flagged, seed, elapsed), flush=True)
    else:
        print("  %-13s %-12s %s: %s" % (dataset_name, configuration_name,
                                        status, note), flush=True)


def main(arguments):
    import powerlaw
    base_seed, reps = read_supplementary_settings()
    datasets = [a for a in arguments if a in DATASETS] or list(DATASETS)
    configurations = ([a for a in arguments if a in CONFIGURATIONS]
                      or list(CONFIGURATIONS))
    (OUTPUT_DIR / "fits").mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "replicates").mkdir(parents=True, exist_ok=True)
    print("powerlaw %s, numpy %s, Python %s; %d bootstrap repetitions"
          % (powerlaw.__version__, np.__version__, sys.version.split()[0], reps),
          flush=True)
    for dataset_name in datasets:
        for configuration_name in configurations:
            run_configuration(dataset_name, configuration_name, base_seed, reps)


if __name__ == "__main__":
    main(sys.argv[1:])
