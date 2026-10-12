# SPDX-License-Identifier: GPL-3.0-only
"""The five SCORCH stages, in order.

Each stage module has a ``main(config=None, log=print)`` function that reads
the results of the earlier stages from the configured folders and writes its
own. ``run_stages`` runs a selection of them in order with the same settings.
"""
from __future__ import annotations

import importlib
import time


class StageError(RuntimeError):
    """A stage found inconsistent data and stopped."""


# (number, record name, module, description)
STAGES = (
    (1, "01_detect_heatwaves", "scorch.stages.detect_heatwaves",
     "Detect local heatwaves at every grid cell (Catalog 1)"),
    (2, "02_select_regional_days", "scorch.stages.select_regional_days",
     "Select regionally extensive heatwave days (Catalog 2)"),
    (3, "03_cluster_heatwaves", "scorch.stages.cluster_heatwaves",
     "Cluster the heatwave cells of each selected day with DBSCAN (Catalog 3)"),
    (4, "04_fit_ellipses", "scorch.stages.fit_ellipses",
     "Fit PCA ellipses and temperature-weighted centroids (Catalog 4)"),
    (5, "05_summarize_compound_events", "scorch.stages.summarize_compound_events",
     "Classify compound events and write Table 1 (Catalog 5)"),
)


def stage_module(number):
    for stage_number, _, module_name, _ in STAGES:
        if stage_number == int(number):
            return importlib.import_module(module_name)
    raise ValueError(f"there is no stage {number}; the stages are 1 to 5")


def parse_stage_selection(text):
    """Turn '1-3', '2,4', 'all' or '5' into a sorted list of stage numbers."""
    text = (text or "all").strip().lower()
    if text in ("all", ""):
        return [number for number, *_ in STAGES]
    chosen = set()
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            first, last = part.split("-", 1)
            chosen.update(range(int(first), int(last) + 1))
        elif part:
            chosen.add(int(part))
    for number in chosen:
        if number not in range(1, len(STAGES) + 1):
            raise ValueError(f"there is no stage {number}; the stages are 1 to 5")
    return sorted(chosen)


def run_stages(config, stage_numbers=None, log=print):
    """Run the chosen stages in order; return a list of (number, summary)."""
    stage_numbers = stage_numbers or [number for number, *_ in STAGES]
    results = []
    for number, record_name, module_name, description in STAGES:
        if number not in stage_numbers:
            continue
        log(f"\nStage {number}/5: {description}")
        started = time.time()
        module = importlib.import_module(module_name)
        summary = module.main(config, log=log)
        log(f"stage {number} finished in {time.time() - started:.1f} s")
        results.append((number, summary))
    return results
