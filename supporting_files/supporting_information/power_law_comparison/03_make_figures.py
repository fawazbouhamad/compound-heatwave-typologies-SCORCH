"""Draw Figures S2 and S3 from the power-law comparisons.

Reads: Prepared maximum-area tables and supplementary comparison results.
Writes: supporting_files/supporting_information/power_law_comparison/results/ and supporting_files/supporting_information/figures/
Run: python supporting_files/supporting_information/power_law_comparison/03_make_figures.py
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
from figure_exports import save_figure, _grid  # noqa: E402

INPUT_DIR = REPO_DIR / "supporting_files/additional_results/modeling/power_law"
OUTPUT_DIR = REPO_DIR / "supporting_files/supporting_information/power_law_comparison/results"
FIGURES_DIR = REPO_DIR / "supporting_files" / "supporting_information" / "figures"
RESULTS_DIR = OUTPUT_DIR / "fits"

# (dataset, label in the table, N, area file, figure name)
DATASETS = [
    ("daily_maxima", "Largest ellipse per event-day", 395, "daily_max_area.csv",
     "figure_S2"),
    ("event_maxima", "Largest ellipse per event", 51, "event_max_area.csv",
     "figure_S3"),
]

# (configuration file name, id in the table, description, settings status,
#  panel letter, panel title)
CONFIGURATIONS = [
    ("M1", "M1", "M1 Clauset reference algorithm (R translation; lower-rank KS)",
     "reference algorithm", "(a)", "Lower-rank KS method"),
    ("M2", "M2", "M2 Full-KS variant (same R code; full continuous KS)",
     "reference algorithm", "(b)", "Full KS method"),
    ("M3_default", "M3 default", "M3 poweRlaw 1.0.0 (default fitting settings)",
     "default fitting settings", "(c)", "R poweRlaw: default"),
    ("M3_adjusted", "M3 adjusted", "M3 poweRlaw 1.0.0 (adjusted: xmax = Inf)",
     "ADJUSTED fitting settings", "(d)", "R poweRlaw: adjusted"),
    ("M4_default", "M4 default", "M4 powerlaw 2.0.0 (default fitting settings)",
     "default fitting settings", "(e)", "Python powerlaw: default"),
    ("M4_adjusted", "M4 adjusted",
     "M4 powerlaw 2.0.0 (adjusted: alpha bounds [1, None])",
     "ADJUSTED fitting settings", "(f)", "Python powerlaw: adjusted"),
]

TABLE_COLUMNS = [
    "dataset", "dataset_label", "n", "config_id", "configuration",
    "settings_status", "n_tail", "xmin_km2", "alpha", "distance_definition",
    "distance", "p_value", "p_lower", "p_upper", "bootstrap_complete", "reps",
    "n_undefined", "n_flagged", "mc_se", "p_ci95_low", "p_ci95_high", "seed",
    "seconds", "software", "status", "flags", "note"]

# style of Figure 11
C_EMP = "black"        # empirical points
C_FIT = "#b2182b"      # fitted tail
C_AMIN = "black"       # dashed A_min line
GRID = "0.90"          # major grid


def number_or_none(value):
    text = "" if value is None else str(value).strip()
    if text in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        x = float(text)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def wilson_interval(p, n):
    """Wilson score 95% interval for a proportion p from n repetitions."""
    if p is None or not n:
        return (None, None)
    z = 1.959963984540054
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def read_result(dataset, configuration):
    path = RESULTS_DIR / ("%s_%s.csv" % (dataset, configuration))
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return rows[0] if rows else None


def read_areas(file_name):
    """area_km2 column, with Python's exact number parser."""
    with open(INPUT_DIR / file_name, newline="", encoding="utf-8") as fh:
        return np.array([float(row["area_km2"]) for row in csv.DictReader(fh)])


# combined table
def build_table():
    rows = []
    for dataset, label, n, _file, _figure in DATASETS:
        for name, config_id, description, settings, _letter, _title in CONFIGURATIONS:
            r = read_result(dataset, name)
            row = {c: "" for c in TABLE_COLUMNS}
            row.update(dataset=dataset, dataset_label=label, n=n,
                       config_id=config_id, configuration=description,
                       settings_status=settings)
            if r is None:
                row["status"] = "not_run"
                rows.append(row)
                continue

            reps = int(number_or_none(r.get("reps")) or 0)
            complete = str(r.get("bootstrap_complete", "TRUE")).upper() != "FALSE"
            p = number_or_none(r.get("p_value"))
            low, high = wilson_interval(p, reps) if complete else (None, None)

            flags = []
            if str(r.get("obs_fit_noise_flag", "")).strip().lower() == "true":
                flags.append("package fell back (no valid cutoff candidate)")
            if str(r.get("obs_pl_noise_flag", "")).strip().lower() == "true":
                flags.append("selected fit noise-flagged")
            if str(r.get("obs_pl_in_range", "")).strip().lower() == "false":
                flags.append("selected fit outside parameter range")

            if name.startswith("M3") or name.startswith("M4"):
                software = ("poweRlaw " if name.startswith("M3") else "powerlaw ") \
                    + (r.get("package_version") or "")
            else:
                software = r.get("r_version", "")

            row.update(
                n_tail=r.get("n_tail", "") or "",
                xmin_km2=r.get("xmin", "") or "",
                alpha=r.get("alpha", "") or "",
                distance_definition=r.get("distance_definition", ""),
                distance=r.get("distance", "") or "",
                p_value=("" if p is None else repr(p)),
                p_lower=r.get("p_lower", "") or "",
                p_upper=r.get("p_upper", "") or "",
                bootstrap_complete=("" if r.get("status") == "no_usable_fit"
                                    else str(complete).upper()),
                reps=(reps or ""),
                n_undefined=(r.get("n_undefined") or r.get("n_failed") or ""),
                n_flagged=r.get("n_flagged", "") or "",
                mc_se=r.get("mc_se", "") or "",
                p_ci95_low=("" if low is None else "%.6f" % low),
                p_ci95_high=("" if high is None else "%.6f" % high),
                seed=r.get("seed", ""), seconds=r.get("seconds", ""),
                software=software.strip(), status=r.get("status", ""),
                flags="; ".join(flags), note=(r.get("note", "") or "").strip())
            rows.append(row)
    return rows


# figures
def apply_style():
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11.5,
        "axes.titlesize": 11.5,
        "axes.labelsize": 12.5,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.labelcolor": "black",
        "text.color": "black",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.antialiased": True,
    })


def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")


def empirical_ccdf(areas):
    """Sorted areas and P(area >= A) at each of them."""
    s = np.sort(areas)
    return s, 1.0 - np.arange(s.size) / s.size


def fitted_tail(xmin, alpha, n, n_tail, xmax):
    """Fitted tail from A_min to xmax, scaled to the observed tail fraction."""
    x = np.logspace(math.log10(xmin), math.log10(xmax), 200)
    return x, (n_tail / n) * (x / xmin) ** (1.0 - alpha)


def area_label(v):
    """Area as "m x 10^e km^2" with three significant figures."""
    exponent = int(math.floor(math.log10(abs(v))))
    mantissa = v / (10.0 ** exponent)
    return r"$%.2f \times 10^{%d}$ km$^2$" % (mantissa, exponent)


def p_value_lines(row, reps_default=5000):
    """p-value text; never rounds a nonzero p to zero."""
    reps = int(number_or_none(row.get("reps")) or reps_default)
    complete = str(row.get("bootstrap_complete", "TRUE")).upper() != "FALSE"
    if not complete:
        low, high = number_or_none(row.get("p_lower")), number_or_none(row.get("p_upper"))
        n_undefined = int(number_or_none(row.get("n_undefined")) or 0)
        return [r"$p$ incomplete: [%.4f, %.4f]" % (low, high),
                "(%d undefined refits)" % n_undefined]
    p = number_or_none(row.get("p_value"))
    if p is None:
        return [r"$p$ unavailable"]
    if p == 0.0:
        # No exceedance in reps repetitions: shown as a count with the
        # one-sided 95% upper bound (about 3 / reps), not as an exact zero.
        return [r"$p$ = 0/%d (95%% upper %.4f)" % (reps, 3.0 / reps)]
    return [r"$p$ = %.4f" % p]


def draw_panel(ax, letter, title, areas, row, xlim, ylim):
    x, c = empirical_ccdf(areas)
    ax.loglog(x, c, ".", ms=4.0, color=C_EMP, alpha=0.85, zorder=3)

    lines = []
    usable = (row is not None
              and row.get("status") in ("ok", "incomplete_bootstrap")
              and number_or_none(row.get("xmin")) is not None)
    if usable:
        xmin = number_or_none(row["xmin"])
        alpha = number_or_none(row["alpha"])
        n_tail = int(number_or_none(row["n_tail"]))
        n = int(number_or_none(row["n"]))
        distance = number_or_none(row["distance"])

        xt, ct = fitted_tail(xmin, alpha, n, n_tail, areas.max())
        ax.loglog(xt, ct, "-", color=C_FIT, lw=2.2, zorder=4)
        ax.axvline(xmin, ls="--", color=C_AMIN, lw=1.1, zorder=2)

        lines = [r"$n_{\mathrm{tail}}$ = %d" % n_tail,
                 r"$A_{\min}$ = " + area_label(xmin),
                 r"$\alpha$ = %.3f" % alpha,
                 r"$D$ = %.4f" % distance]
        lines += p_value_lines(row)
    else:
        # Show a status message when no usable fit is available.
        message = "configuration not run" if row is None else "No usable fit"
        ax.text(0.5, 0.55, message, transform=ax.transAxes, ha="center",
                va="center", fontsize=12, color=C_FIT)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"Area (km$^2$)")
    ax.set_ylabel(r"$P(\mathrm{area} \geq A)$")
    ax.set_title(title, fontsize=11.5, pad=7)
    ax.grid(True, which="major", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    despine(ax)

    # Label decades only: the data span about one decade, and labelled minor
    # ticks would overlap.
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

    ax.text(-0.17, 1.17, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top", ha="left", color="black")

    if lines:
        ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes,
                va="bottom", ha="left", fontsize=10, color="black",
                linespacing=1.45)


def build_figure(dataset, area_file):
    areas = read_areas(area_file)
    n = areas.size
    rows = [read_result(dataset, name) for name, *_rest in CONFIGURATIONS]

    # The same axis ranges in all six panels.
    x, c = empirical_ccdf(areas)
    xlim = (10 ** (math.log10(x.min()) - 0.06), 10 ** (math.log10(x.max()) + 0.10))
    lowest = min(c.min(), 1.0 / n)
    ylim = (10 ** (math.log10(lowest) - 0.35), 1.8)

    apply_style()
    fig, axes = plt.subplots(3, 2, figsize=(10.4, 13.2))
    for ax, (_name, _cid, _desc, _settings, letter, title), row in zip(
            axes.ravel(), CONFIGURATIONS, rows):
        draw_panel(ax, letter, title, areas, row, xlim, ylim)
    fig.subplots_adjust(left=0.095, right=0.975, top=0.965, bottom=0.055,
                        wspace=0.28, hspace=0.40)
    return fig


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    table = build_table()
    with open(OUTPUT_DIR / "combined_results_table.csv", "w", newline="\n",
              encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    print("table rows: %d; usable fits: %d; no usable fit: %d; errors: %d; "
          "not run: %d" % (
              len(table),
              sum(r["status"] in ("ok", "incomplete_bootstrap") for r in table),
              sum(r["status"] == "no_usable_fit" for r in table),
              sum(r["status"] == "error" for r in table),
              sum(r["status"] == "not_run" for r in table)))

    for dataset, _label, _n, area_file, figure_name in DATASETS:
        fig = build_figure(dataset, area_file)
        paths = save_figure(fig, FIGURES_DIR / figure_name, dpi=600,
                            panel_layout=_grid([0, 0.5, 1], [0, 1 / 3, 2 / 3, 1]),
                            bbox_inches="tight", facecolor="white")
        for path in paths:
            print("wrote %s" % path)
        plt.close(fig)

    print("matplotlib %s" % mpl.__version__)


if __name__ == "__main__":
    main()
