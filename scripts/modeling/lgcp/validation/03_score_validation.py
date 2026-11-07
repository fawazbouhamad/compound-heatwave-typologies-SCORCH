"""Score the held-out centroids against each fold's intensity zones.

Reads: Prepared LGCP inputs and Appendix C fold predictions.
Writes: supporting_files/additional_results/modeling/lgcp/validation/
Run: python scripts/modeling/lgcp/validation/03_score_validation.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
sys.path.insert(0, str(REPO_DIR / "scripts" / "modeling" / "lgcp"))
import functions as f  # noqa: E402
import zone_scoring as zs  # noqa: E402

OUTPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp" / "validation"
LGCP_INPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp"


def add_centroid_columns(per_centroid, centroids):
    """Put identifying columns of the centroids in front of the scores."""
    info = centroids.set_index("ellipse_id").loc[per_centroid["ellipse_id"]]
    for position, column in enumerate(["compound_event_id", "date", "centroid_lon",
                                       "centroid_lat", "is_daily_largest", "is_event_largest"],
                                      start=1):
        per_centroid.insert(position, column, info[column].to_numpy())
    return per_centroid


def main():
    started = time.time()
    config = f.load_config()
    zone_fractions = config["appendix_c"]["zone_area_fractions"]
    n_folds = config["appendix_c"]["n_folds"]

    cells = pd.read_csv(LGCP_INPUT_DIR / "temperature_covariates.csv")
    cells["area_km2"] = zs.cell_area_km2(cells["lat"])
    # The manuscript distances use coordinates parsed by pandas' default
    # number reader, so keep that precision here.
    centroids = pd.read_csv(LGCP_INPUT_DIR / "centroids.csv")
    row_of_cell = pd.Series(np.arange(len(cells)), index=cells["cell_id"])
    centroids["cell_row"] = row_of_cell.reindex(centroids["nearest_cell_id"]).to_numpy()
    folds = pd.read_csv(OUTPUT_DIR / "fold_assignment.csv")
    if not np.array_equal(folds["ellipse_id"], centroids["ellipse_id"]):
        raise ValueError("fold_assignment.csv is not in centroid order")
    centroids["fold"] = folds["fold"].to_numpy()
    fold_surfaces = pd.read_csv(OUTPUT_DIR / "fold_mean_intensity_1deg.csv")
    full_surface = pd.read_csv(OUTPUT_DIR / "full_data_mean_intensity_1deg.csv")

    # cross-validation: each fold's withheld centroids on that fold's model
    per_fold, zone_rows = [], []
    for fold in range(1, n_folds + 1):
        surface = (fold_surfaces[fold_surfaces["fold"] == fold]
                   .set_index("cell_id").reindex(cells["cell_id"]))
        intensity = surface["mean_intensity_per_km2"].to_numpy(float)
        withheld = centroids[centroids["fold"] == fold]
        scores, zones = zs.score_centroids(intensity, cells, withheld, zone_fractions)
        scores.insert(1, "fold", fold)
        zones.insert(0, "fold", fold)
        per_fold.append(scores)
        zone_rows.append(zones)
    validation = pd.concat(per_fold).sort_values("ellipse_id").reset_index(drop=True)
    if not (len(validation) == len(centroids) and validation["ellipse_id"].is_unique):
        raise ValueError("every centroid must be withheld exactly once")
    validation = add_centroid_columns(validation, centroids)
    validation_zones = pd.concat(zone_rows).reset_index(drop=True)
    validation_groups = zs.summarize_groups(validation, zone_fractions)

    # the model fitted to all 753 centroids, scored on the same 753
    # Convert the stored intensity per million km^2 to intensity per km^2.
    full_intensity = full_surface["mean_intensity_per_million_km2"].to_numpy(float) / 1e6
    full_data, _ = zs.score_centroids(full_intensity, cells, centroids, zone_fractions)
    full_data = add_centroid_columns(full_data, centroids)
    full_groups = zs.summarize_groups(full_data, zone_fractions)

    # full-data fit against cross-validation
    rows = []
    for group in ("all_centroids", "daily_largest", "event_largest"):
        full_row = full_groups.set_index("group").loc[group]
        validation_row = validation_groups.set_index("group").loc[group]
        statistics = [f"share_in_{zs.zone_label(x)}_zone" for x in zone_fractions] + \
                     [f"mean_distance_to_{zs.zone_label(x)}_zone_km" for x in zone_fractions] + \
                     ["mean_area_weighted_rank"]
        for statistic in statistics:
            rows.append({"group": group, "statistic": statistic,
                         "full_data_fit": full_row[statistic],
                         "cross_validation": validation_row[statistic],
                         "cross_validation_minus_full_data_fit":
                             validation_row[statistic] - full_row[statistic]})
    full_vs_validation = pd.DataFrame(rows)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {"validation_per_centroid.csv": validation,
               "validation_group_summary.csv": validation_groups,
               "validation_zones.csv": validation_zones,
               "full_data_vs_validation.csv": full_vs_validation}
    for name, table in outputs.items():
        table.to_csv(OUTPUT_DIR / name, index=False, lineterminator="\n")

    numbers = appendix_c_numbers(full_groups, validation_groups)
    numbers.to_csv(OUTPUT_DIR / "appendix_c_numbers.csv", index=False, lineterminator="\n")
    print(numbers.to_string(index=False))

    numbers_ok = bool(numbers["matches"].all())
    print(f"Appendix C numbers reproduced: {numbers_ok}")
    print(f"finished in {time.time() - started:.1f} s")
    if not numbers_ok:
        sys.exit(1)


def appendix_c_numbers(full_groups, validation_groups):
    """The Appendix C numbers, rounded as printed, against the manuscript (final corrected analysis)."""
    full = full_groups.set_index("group")
    validation = validation_groups.set_index("group")
    printed = [
        # (statistic, group, model: full or validation, text in the manuscript)
        ("share in top 10% zone", "all_centroids", "full", "19.0%"),
        ("share in top 20% zone", "all_centroids", "full", "35.5%"),
        ("share in top 30% zone", "all_centroids", "full", "46.3%"),
        ("share in top 50% zone", "all_centroids", "full", "71.3%"),
        ("share in top 10% zone", "all_centroids", "validation", "19.0%"),
        ("share in top 20% zone", "all_centroids", "validation", "35.3%"),
        ("share in top 30% zone", "all_centroids", "validation", "46.5%"),
        ("share in top 50% zone", "all_centroids", "validation", "70.7%"),
        ("mean distance to top 10% zone", "all_centroids", "full", "332.3 km"),
        ("mean distance to top 20% zone", "all_centroids", "full", "214.7 km"),
        ("mean distance to top 30% zone", "all_centroids", "full", "147.1 km"),
        ("mean distance to top 50% zone", "all_centroids", "full", "62.8 km"),
        ("mean distance to top 10% zone", "all_centroids", "validation", "347.4 km"),
        ("mean distance to top 20% zone", "all_centroids", "validation", "215.9 km"),
        ("mean distance to top 30% zone", "all_centroids", "validation", "148.9 km"),
        ("mean distance to top 50% zone", "all_centroids", "validation", "64.8 km"),
        ("share in top 20% zone (conclusion)", "all_centroids", "validation", "35%"),
        ("share in top 20% zone (conclusion)", "event_largest", "validation", "55%"),
    ]
    rows = []
    for statistic, group, source, manuscript in printed:
        table = full if source == "full" else validation
        zone = statistic.split("top ")[1].split("%")[0]
        if statistic.startswith("share"):
            value = 100 * table.loc[group, f"share_in_top{zone}_zone"]
            text = f"{value:.0f}%" if "conclusion" in statistic else f"{value:.1f}%"
        else:
            value = table.loc[group, f"mean_distance_to_top{zone}_zone_km"]
            text = f"{value:.1f} km"
        rows.append({"statistic": statistic, "group": group,
                     "model": "fitted to all 753 centroids" if source == "full"
                     else "cross-validation (withheld centroids)",
                     "reproduced_value": value, "reproduced_as_printed": text,
                     "manuscript": manuscript, "matches": text == manuscript})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
