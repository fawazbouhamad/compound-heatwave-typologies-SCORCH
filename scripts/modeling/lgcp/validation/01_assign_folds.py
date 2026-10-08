"""Assign the centroid observations to five validation folds.

Reads: supporting_files/additional_results/modeling/lgcp/centroids.csv and shared Python study settings.
Writes: supporting_files/additional_results/modeling/lgcp/validation/
Run: python scripts/modeling/lgcp/validation/01_assign_folds.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_DIR / "scripts" / "main_workflow" / "helpers"))
import functions as f  # noqa: E402

OUTPUT_DIR = REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp" / "validation"


def main():
    started = time.time()
    config = f.load_config()
    n_folds = config["appendix_c"]["n_folds"]
    seed = config["appendix_c"]["fold_seed"]

    # Use pandas' default number reader for the centroid coordinates.
    # It can change the last binary digit (here by up to 5e-13 km).
    # The manuscript fold fits use these parsed positions, so keep this
    # reading step before writing the fold assignments.
    centroids = pd.read_csv(REPO_DIR / "supporting_files" / "additional_results" / "modeling" / "lgcp" / "centroids.csv")
    n_centroids = len(centroids)

    shuffled_rows = np.random.default_rng(seed).permutation(n_centroids)
    fold_of_row = np.empty(n_centroids, dtype=int)
    for fold, rows in enumerate(np.array_split(shuffled_rows, n_folds), start=1):
        fold_of_row[rows] = fold

    folds = pd.DataFrame({
        "ellipse_id": centroids["ellipse_id"],
        "row_in_catalog": np.arange(n_centroids),
        "fold": fold_of_row,
        "compound_event_id": centroids["compound_event_id"],
        "date": centroids["date"],
        "centroid_lon": centroids["centroid_lon"],
        "centroid_lat": centroids["centroid_lat"],
        "x_km": centroids["x_km"],
        "y_km": centroids["y_km"],
        "nearest_cell_id": centroids["nearest_cell_id"],
    })
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    folds.to_csv(OUTPUT_DIR / "fold_assignment.csv", index=False, lineterminator="\n")
    sizes = folds["fold"].value_counts().sort_index()
    print("centroids per fold:", {int(k): int(v) for k, v in sizes.items()})

    # np.array_split gives fold sizes that differ by at most one (753 centroids: 151, 151, 151, 150, 150).
    if not (len(sizes) == n_folds and sizes.max() - sizes.min() <= 1):
        raise SystemExit("the fold sizes differ by more than one")
    print(f"finished in {time.time() - started:.1f} s")


if __name__ == "__main__":
    main()
