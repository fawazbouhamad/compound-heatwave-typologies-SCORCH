#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared fonts, colours, panel labels and map layers for the figures.

Original fonts can be supplied with SCORCH_SLIDE_FONT_DIR; map layers are bundled.
"""
from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import warnings

from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

# fonts
# Typology headings use Abadi.
# Strict matching is optional; set SCORCH_REQUIRE_SLIDE_FONTS=1 to enable it.
_SLIDE_FONT_DIR_VAR = "SCORCH_SLIDE_FONT_DIR"
_REQUIRE_VAR = "SCORCH_REQUIRE_SLIDE_FONTS"
_FALLBACK_FAMILY = "Arial"

# Font files recorded for this run.
SLIDE_FONT_FILES: list = []


def slide_font_dir():
    """Return the explicitly supplied font directory, if any."""
    explicit = os.environ.get(_SLIDE_FONT_DIR_VAR, "").strip()
    return Path(explicit).expanduser().resolve() if explicit else None


def register_slide_fonts() -> str:
    """Use original fonts when available; require exact files only in strict mode."""
    del SLIDE_FONT_FILES[:]
    directory = slide_font_dir()
    strict = os.environ.get(_REQUIRE_VAR, "").strip().lower() in {"1", "true", "yes"}
    if strict:
        if directory is None:
            raise RuntimeError("Exact figures require SCORCH_SLIDE_FONT_DIR to point to the recorded Abadi fonts.")
        manifest = Path(__file__).resolve().parents[3] / "supporting_files/additional_inputs/figure_assets/font_files.csv"
        with manifest.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                path = directory / Path(row["file"]).name
                if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
                    raise RuntimeError(f"The recorded font is missing or different: {path}")
    for ttf in (sorted(directory.glob("*.ttf")) if directory else []):
        try:
            fm.fontManager.addfont(str(ttf))
            SLIDE_FONT_FILES.append(str(ttf))
        except Exception as error:
            if strict:
                raise RuntimeError(f"Could not load the recorded font: {ttf}") from error
    names = {f.name for f in fm.fontManager.ttflist}
    if "Abadi" in names and (SLIDE_FONT_FILES or not strict):
        return "Abadi"
    if strict:
        raise RuntimeError(f"The recorded Abadi fonts could not be registered from {directory}.")
    fallback = _FALLBACK_FAMILY if _FALLBACK_FAMILY in names else "DejaVu Sans"
    warnings.warn(f"Abadi is unavailable; figure headings use {fallback}. Scientific calculations are unchanged.",
                  stacklevel=2)
    return fallback


HEADING_FAMILY = register_slide_fonts()


# map data
def use_distributed_map_data() -> None:
    """Point cartopy to the Natural Earth layers in supporting_files/additional_inputs/maps/cartopy/."""
    import cartopy
    repo = Path(__file__).resolve().parents[3]
    import functions as f
    config = f.load_config()
    map_dir = f.supporting_input_dir(config) / "maps" / "cartopy"
    cartopy.config["pre_existing_data_dir"] = map_dir if map_dir.is_absolute() else repo / map_dir


use_distributed_map_data()

# shared constants
PANEL_LABEL = dict(fontsize=17, fontweight="bold", color="black",
                   family="Arial")
LABEL_PAD_PT = 6          # gap between panel edge and label, in points
TYPE_COLORS = {1: "#FF0000", 2: "#FF7401", 3: "#FFBC01", 4: "#7030A0"}
DATE_FS = 12
LEGEND_FS = 12


def _pt_to_fig(fig, pts: float, axis: str) -> float:
    size_in = fig.get_size_inches()[0 if axis == "x" else 1]
    return (pts / 72.0) / size_in


def label_above(fig, ax, letter: str, pad_pt: float = LABEL_PAD_PT) -> None:
    """Panel label such as "(a)" centred above a panel."""
    b = ax.get_position()
    fig.text((b.x0 + b.x1) / 2.0, b.y1 + _pt_to_fig(fig, pad_pt, "y"),
             letter, ha="center", va="bottom", **PANEL_LABEL)


def type_heading(fig, x, y, type_n: int, fontsize: float = 26,
                 text: str | None = None) -> None:
    """Typology heading in bold Abadi."""
    fig.text(x, y, text or f"Type {type_n}", ha="center", va="center",
             fontsize=fontsize, fontweight="bold", family=HEADING_FAMILY,
             color=TYPE_COLORS[type_n])


def date_label(ax, text: str, fontsize: float = DATE_FS) -> None:
    """Date inside the bottom-right of a snapshot: black, bold, fixed inset."""
    ax.text(0.975, 0.03, text, transform=ax.transAxes, ha="right",
            va="bottom", fontsize=fontsize, fontweight="bold",
            color="black", zorder=20,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none",
                      pad=1.5))


def snapshot_legend_handles(subscript_axis_labels: bool = False,
                            axis_labels=None):
    """L1 / L2 / Centroid legend entries (line styles of the snapshot maps).

    The two ellipse-axis entries are written as plain "L1" / "L2" by
    default. subscript_axis_labels=True writes them with upright subscript
    numerals (as in Figure S1); axis_labels=(text1, text2) gives the two
    texts explicitly (Figures 5-7 use the italic mathtext "$L_1$", "$L_2$").
    """
    if axis_labels is not None:
        l1, l2 = axis_labels
    elif subscript_axis_labels:
        l1, l2 = r"$\mathregular{L_1}$", r"$\mathregular{L_2}$"
    else:
        l1, l2 = "L1", "L2"
    return [
        Line2D([0], [0], color="black", lw=2.4, label=l1),
        Line2D([0], [0], color="black", lw=1.5, label=l2),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markeredgecolor="k", markersize=7, label="Centroid"),
    ]


def unified_snapshot_legend(fig, y: float = 0.012,
                            fontsize: float = LEGEND_FS,
                            subscript_axis_labels: bool = False,
                            axis_labels=None):
    """One horizontal legend centred below the whole figure (Figs 5-7, S1)."""
    return fig.legend(
        handles=snapshot_legend_handles(subscript_axis_labels, axis_labels),
        loc="lower center",
        bbox_to_anchor=(0.5, y), ncol=3, frameon=True,
        fontsize=fontsize, edgecolor="black", facecolor="white",
        handlelength=2.2, columnspacing=1.6, borderpad=0.55)
