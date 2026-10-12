# SPDX-License-Identifier: GPL-3.0-only
"""Choose DBSCAN settings and group the heatwave cells on selected days.

Reads Catalog 2; writes clusters and cell membership in Catalog 3.
Parameter searches and consecutive-date sequence settings are supporting results.
The calculations are in src/scorch/stages/cluster_heatwaves.py (stage 3).
Run: python scripts/main_workflow/03_cluster_heatwaves.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions  # noqa: E402,F401  (makes the scorch package importable)
from scorch.stages.cluster_heatwaves import main  # noqa: E402

if __name__ == "__main__":
    main()
