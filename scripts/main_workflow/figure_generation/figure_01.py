"""Figure 1: workflow schematic.

Exports the editable PowerPoint to output/figures/figure_1/figure_1.png.
Requires desktop PowerPoint on Windows. Figure 1 has no separate panels.
Run: python scripts/main_workflow/figure_generation/figure_01.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("SCORCH_REQUIRE_SLIDE_FONTS", "0")
WORKFLOW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW))
sys.path.insert(0, str(WORKFLOW / "helpers"))

import matplotlib as mpl
import functions as f

from figure_exports import export_powerpoint


def main():
    config = f.load_config()
    inputs = {
        "artwork": f.supporting_input_dir(config) / "figure_assets/figure_01_source.pptx",
    }

    out_stem = f.output_subdir(config, "figures") / "figure_1"
    with mpl.rc_context():
        mpl.rcdefaults()
        written = export_powerpoint(inputs["artwork"], out_stem)
    f.write_run_record(config, __file__, list(inputs.values()), written)
    for path in written:
        print(f.relative_to_repo(path))
    return written


if __name__ == "__main__":
    main()
