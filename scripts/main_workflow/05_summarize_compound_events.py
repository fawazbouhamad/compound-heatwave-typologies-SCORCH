# SPDX-License-Identifier: GPL-3.0-only
"""Classify compound events and calculate the manuscript summaries.

Reads the preceding catalogs; writes Catalog 5, output/tables/table_1.csv,
and supporting calculations in supporting_files/additional_results/main_workflow/.
The calculations are in src/scorch/stages/summarize_compound_events.py (stage 5).
Run: python scripts/main_workflow/05_summarize_compound_events.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
import functions  # noqa: E402,F401  (makes the scorch package importable)
from scorch.stages.summarize_compound_events import main  # noqa: E402

if __name__ == "__main__":
    main()
