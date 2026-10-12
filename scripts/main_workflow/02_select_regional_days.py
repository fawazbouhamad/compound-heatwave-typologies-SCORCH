# SPDX-License-Identifier: GPL-3.0-only
"""Select days with regionally extensive heatwaves and their heatwave cells.

Reads step-1 heatwave flags; writes selected days and cells in Catalog 2.
The full daily summary is saved with the supporting results.
The calculations are in src/scorch/stages/select_regional_days.py (stage 2).
Run: python scripts/main_workflow/02_select_regional_days.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions  # noqa: E402,F401  (makes the scorch package importable)
from scorch.stages.select_regional_days import main  # noqa: E402

if __name__ == "__main__":
    main()
