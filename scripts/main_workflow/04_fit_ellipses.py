# SPDX-License-Identifier: GPL-3.0-only
"""Fit a PCA ellipse and temperature-weighted centroid to each cluster.

Reads Catalog 3 and input temperatures; writes Catalog 4.
The calculations are in src/scorch/stages/fit_ellipses.py (stage 4).
Run: python scripts/main_workflow/04_fit_ellipses.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions  # noqa: E402,F401  (makes the scorch package importable)
from scorch.stages.fit_ellipses import main  # noqa: E402

if __name__ == "__main__":
    main()
