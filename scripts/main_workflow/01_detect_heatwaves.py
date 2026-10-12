# SPDX-License-Identifier: GPL-3.0-only
"""Find local heatwaves using each cell's warm-season 95th percentile.

Reads input/tmax_processed.nc and supporting_files/additional_inputs/maps/grid_cells.csv.
Writes the local-event catalog; saves flags and thresholds as supporting results.
The calculations are in src/scorch/stages/detect_heatwaves.py (stage 1).
Run: python scripts/main_workflow/01_detect_heatwaves.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions  # noqa: E402,F401  (makes the scorch package importable)
from scorch.stages.detect_heatwaves import main  # noqa: E402

if __name__ == "__main__":
    main()
