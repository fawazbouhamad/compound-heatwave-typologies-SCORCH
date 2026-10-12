"""Draw Figure 11 from the power-law results.

Reads: supporting_files/additional_results/modeling/power_law/
Writes: output/figures/figure_11/
Run: python scripts/modeling/power_law/03_make_figure_11.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_DIR = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
from figure_exports import save_figure  # noqa: E402

INPUT_DIR = REPO_DIR / "supporting_files/additional_results/modeling/power_law"
OUTPUT_DIR = REPO_DIR / "supporting_files/additional_results/modeling/power_law"
FIGURES_DIR = REPO_DIR / "output" / "figures"

# Fixed y-axis limits of the manuscript boxplot panels.
# Keeping these limits preserves the tick and grid positions.
BOXPLOT_Y_LIMITS = {
    "daily_maxima": (4.273542835299951, 16.793592078509292),   # panel (b)
    "event_maxima": (2.6994359805698926, 10.398797494012225),  # panel (d)
}

# colours and style of the manuscript figure
C_EMP = "black"       # empirical points
C_FIT = "#b2182b"     # fitted power-law tail
C_AMIN = "black"      # dashed A_min line
C_BOX = "0.90"        # box fill
C_BOXEDGE = "0.35"    # box edge
GRID = "0.90"         # major grid


def empirical_ccdf(areas):
    """Sorted areas and P(area >= A) at each of them."""
    a = np.sort(areas)
    n = a.size
    return a, 1.0 - np.arange(n) / n


def fitted_tail(xmin, alpha, n, n_tail, xmax):
    """Fitted tail from A_min to xmax, scaled to the observed tail fraction."""
    x = np.logspace(np.log10(xmin), np.log10(xmax), 200)
    return x, (n_tail / n) * (x / xmin) ** (1.0 - alpha)


def apply_style():
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 12.5,
        "axes.titlesize": 14,
        "axes.labelsize": 13.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11.5,
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


def panel_label(ax, letter):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=17, fontweight="bold", va="top", ha="left",
            color="black", family="Arial")


def ccdf_panel(ax, sample, stats_lines, show_amin_legend):
    x, c = empirical_ccdf(sample["areas"])
    ax.loglog(x, c, ".", ms=4.5, color=C_EMP, alpha=0.85,
              label="_nolegend_", zorder=3)

    xt, ct = fitted_tail(sample["xmin"], sample["alpha"], sample["n"],
                         sample["n_tail"], x.max())
    ax.loglog(xt, ct, "-", color=C_FIT, lw=2.2, label="_nolegend_", zorder=4)

    ax.axvline(sample["xmin"], ls="--", color=C_AMIN, lw=1.1,
               label=r"$A_{\min}$", zorder=2)

    ax.set_xlabel(r"Area (km$^2$)")
    ax.set_ylabel(r"$P(\mathrm{area} \geq A)$")
    ax.grid(True, which="major", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    despine(ax)

    if show_amin_legend:
        ax.legend(loc="upper right", frameon=False, handlelength=1.6,
                  borderaxespad=0.3)

    ax.text(0.03, 0.04, "\n".join(stats_lines), transform=ax.transAxes,
            va="bottom", ha="left", fontsize=11.5, color="black",
            linespacing=1.4)
    return (x, c), (xt, ct)


def alpha_box_panel(ax, boot_alpha):
    ax.boxplot(boot_alpha, positions=[1], widths=0.45, showfliers=False,
               patch_artist=True, whis=(2.5, 97.5),
               medianprops=dict(color="black", lw=1.6),
               boxprops=dict(facecolor=C_BOX, edgecolor=C_BOXEDGE, lw=1.0),
               whiskerprops=dict(color=C_BOXEDGE, lw=1.0),
               capprops=dict(color=C_BOXEDGE, lw=1.0))
    ax.set_xlim(0.4, 1.6)
    ax.set_xticks([])
    ax.tick_params(axis="x", which="both", length=0)
    ax.set_ylabel(r"Bootstrap $\alpha$ estimate")
    ax.grid(True, axis="y", which="major", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    despine(ax)


def load_sample(name, areas_file, fits):
    """Areas, fit and bootstrap exponents of one sample, with basic checks."""
    fit = fits.loc[name]
    areas = pd.read_csv(INPUT_DIR / areas_file,
                        converters={"area_km2": float})["area_km2"].to_numpy()
    replicates = pd.read_csv(OUTPUT_DIR / f"bootstrap_replicates_{name}.csv",
                             converters={"alpha": float})
    # Every simulated sample was fitted, none dropped, and the exceedance
    # count agrees with the p-value in power_law_fits.csv.
    assert (replicates["status"] == "ok").all(), f"{name}: undefined refits"
    assert len(replicates) == int(fit["reps"])
    assert int((replicates["exceeds_observed"] == 1).sum()) == int(fit["n_exceed"])
    assert areas.size == int(fit["n"])
    return dict(name=name, areas=areas, xmin=float(fit["xmin_km2"]),
                alpha=float(fit["alpha"]), n=int(fit["n"]),
                n_tail=int(fit["n_tail"]), p_value=float(fit["p_value"]),
                boot_alpha=replicates["alpha"].to_numpy())


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    apply_style()

    fits = pd.read_csv(OUTPUT_DIR / "power_law_fits.csv", dtype=str).set_index("dataset")
    daily = load_sample("daily_maxima", "daily_max_area.csv", fits)
    event = load_sample("event_maxima", "event_max_area.csv", fits)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.4))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    def stats_lines(sample):
        return ["$n$ = %d" % sample["n"],
                "$n_{\\mathrm{tail}}$ = %d" % sample["n_tail"],
                r"$\alpha \approx$ " + "%.2f" % sample["alpha"],
                r"$p \approx$ " + "%.3f" % sample["p_value"]]

    for letter, ax_ccdf, ax_box, sample, legend in (
            ("a", ax_a, ax_b, daily, True), ("c", ax_c, ax_d, event, False)):
        box_letter = "b" if letter == "a" else "d"
        ccdf_panel(ax_ccdf, sample, stats_lines(sample), legend)
        panel_label(ax_ccdf, "(%s)" % letter)
        alpha_box_panel(ax_box, sample["boot_alpha"])
        panel_label(ax_box, "(%s)" % box_letter)

    fig.subplots_adjust(left=0.085, right=0.975, top=0.955, bottom=0.075,
                        wspace=0.28, hspace=0.30)

    ax_b.set_ylim(*BOXPLOT_Y_LIMITS["daily_maxima"])
    ax_d.set_ylim(*BOXPLOT_Y_LIMITS["event_maxima"])

    # The whiskers must lie inside the fixed ranges; nothing may be cut off.
    for ax, sample in ((ax_b, daily), (ax_d, event)):
        low, high = np.percentile(sample["boot_alpha"], [2.5, 97.5])
        y_low, y_high = ax.get_ylim()
        assert y_low <= low and high <= y_high, (
            "%s: whiskers [%.4f, %.4f] outside the fixed axis range [%.4f, %.4f]"
            % (sample["name"], low, high, y_low, y_high))

    paths = save_figure(fig, FIGURES_DIR / "figure_11", dpi=600,
                        bbox_inches="tight")
    plt.close(fig)

    for path in paths:
        print("wrote %s" % path)


if __name__ == "__main__":
    main()
