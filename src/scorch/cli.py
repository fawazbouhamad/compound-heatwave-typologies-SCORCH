# SPDX-License-Identifier: GPL-3.0-only
"""The ``scorch`` command.

    scorch run          run the five stages on a prepared NetCDF file
    scorch demo         build a small synthetic example and run it end to end
    scorch check-input  describe a NetCDF file and check that SCORCH can read it
    scorch settings     print the effective settings (or write a template)
    scorch info         versions of SCORCH and its main packages
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from . import __version__
from . import config as cfg
from . import io
from .stages import StageError, parse_stage_selection, run_stages

PROGRAM = "scorch"


def _add_folder_options(parser):
    parser.add_argument("--settings", metavar="FILE",
                        help="YAML, TOML or JSON file overriding the manuscript settings "
                             "(see config/manuscript_settings.yml)")
    parser.add_argument("--input-dir", metavar="DIR",
                        help="folder holding the processed temperature file (default: input)")
    parser.add_argument("--input-file", metavar="NAME",
                        help="name of the temperature file inside --input-dir "
                             "(default: tmax_processed.nc)")
    parser.add_argument("--output-dir", metavar="DIR",
                        help="folder for catalogs and tables (default: output). When given "
                             "without --supporting-results-dir, supporting results and run "
                             "records are written under it too")
    parser.add_argument("--supporting-input-dir", metavar="DIR",
                        help="folder with maps/grid_cells.csv and the country boundaries "
                             "(default: supporting_files/additional_inputs)")
    parser.add_argument("--supporting-results-dir", metavar="DIR",
                        help="folder for intermediate results "
                             "(default: supporting_files/additional_results)")
    parser.add_argument("--run-records-dir", metavar="DIR",
                        help="folder for the JSON run records (default: .scorch/logs/run-records)")


def config_from_arguments(args):
    """Build the settings from command-line options (options win over the file)."""
    paths = {}
    if args.input_dir:
        paths["input_dir"] = args.input_dir
    if args.input_file:
        paths["input_file"] = args.input_file
    if args.supporting_input_dir:
        paths["supporting_input_dir"] = args.supporting_input_dir
    if args.output_dir:
        paths["output_dir"] = args.output_dir
        if not args.supporting_results_dir:
            paths["supporting_results_dir"] = str(Path(args.output_dir) / "supporting_results")
        if not args.run_records_dir:
            paths["run_records_dir"] = str(Path(args.output_dir) / "run_records")
    if args.supporting_results_dir:
        paths["supporting_results_dir"] = args.supporting_results_dir
    if args.run_records_dir:
        paths["run_records_dir"] = args.run_records_dir
    overrides = {"paths": paths} if paths else None
    return cfg.load_config(settings_file=args.settings, overrides=overrides)


def command_run(args):
    config = config_from_arguments(args)
    stages = parse_stage_selection(args.stages)
    print(f"SCORCH {__version__}")
    print(f"input file:         {cfg.input_file(config)}")
    print(f"output folder:      {cfg.output_dir(config)}")
    print(f"supporting results: {cfg.supporting_results_dir(config)}")
    print(f"run records:        {cfg.run_records_dir(config)}")
    print(f"stages:             {', '.join(map(str, stages))}")
    if stages[0] == 1 and not cfg.input_file(config).is_file():
        raise io.InputError(
            f"input file not found: {cfg.input_file(config)}\n"
            "For the study, download SCORCH_data.zip from https://doi.org/10.5281/zenodo.21717874, "
            "extract it and copy scorch_data/tmax_processed.nc into the input folder. "
            "To try SCORCH without the dataset, run: scorch demo")
    results = run_stages(config, stages)
    print("\nDone. Catalogs are in", cfg.output_dir(config) / "catalogs")
    if 5 in stages:
        print("Table 1 is in", cfg.output_dir(config) / "tables" / "table_1.csv")
        summary = dict(results)[5]
        print("compound events by type:", summary.get("events_by_type"))
    return 0


def command_demo(args):
    from .demo import run_demo
    run_demo(Path(args.dir), seed=args.seed, make_figure=not args.no_figure,
             with_r=args.with_r)
    return 0


def command_check_input(args):
    description = io.describe_netcdf(args.file)
    width = max(len(k) for k in description)
    for key, value in description.items():
        print(f"{key:<{width}}  {value}")
    if args.settings or args.strict:
        config = cfg.load_config(settings_file=args.settings)
        data = io.read_netcdf_tmax(args.file)
        io.check_study_calendar(data["dates"], config)
        expected_cells = config["study_area"].get("expected_cells")
        if expected_cells is not None and description["cells"] != int(expected_cells):
            raise io.InputError(f"the file has {description['cells']} cells; the settings "
                                f"expect {expected_cells}")
        print("\nThe file matches the study calendar and grid of the settings.")
    else:
        print("\nSCORCH can read this file. Add --strict (or --settings FILE) to check it "
              "against the study calendar and grid.")
    return 0


def command_settings(args):
    config = cfg.load_config(settings_file=args.settings)
    text = cfg.settings_as_yaml_text(config)
    if args.write:
        Path(args.write).write_text(text, encoding="utf-8", newline="\n")
        print(f"wrote {args.write}")
    else:
        print(text, end="")
    return 0


def command_info(args):
    versions = io.software_versions()
    width = max(len(k) for k in versions)
    for key, value in versions.items():
        print(f"{key:<{width}}  {value}")
    root = cfg.repository_root()
    print(f"{'repository':<{width}}  {root if root else 'not running from a checkout'}")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description="SCORCH: Spatiotemporal Classification of Regional Compound Heatwaves.",
        epilog="Documentation: README.md and the docs/ folder of the repository.")
    parser.add_argument("--version", action="version", version=f"{PROGRAM} {__version__}")
    parser.add_argument("--debug", action="store_true",
                        help="show the full Python traceback when a command fails")
    commands = parser.add_subparsers(dest="command", metavar="COMMAND")
    commands.required = True

    run = commands.add_parser("run", help="run the five stages on a prepared temperature file",
                              description="Run the SCORCH stages in order. Without options this "
                                          "reproduces the manuscript catalogs from input/tmax_processed.nc.")
    run.add_argument("--stages", default="all", metavar="SELECTION",
                     help="which stages to run: 'all' (default), '1-3', '2,4' or '5'")
    _add_folder_options(run)
    run.set_defaults(handler=command_run)

    demo = commands.add_parser("demo", help="build and run the small synthetic example",
                               description="Create a small synthetic temperature dataset with "
                                           "known heatwave events, run the five stages on it and "
                                           "draw an overview figure. Needs no download.")
    demo.add_argument("--dir", default="demo/run", metavar="DIR",
                      help="folder for the synthetic input and the results (default: demo/run)")
    demo.add_argument("--seed", type=int, default=None,
                      help="random seed of the synthetic weather (default: the recorded demo seed)")
    demo.add_argument("--no-figure", action="store_true", help="skip the overview figure")
    demo.add_argument("--with-r", action="store_true",
                      help="also fit the power law to the demo ellipse areas with R (needs Rscript)")
    demo.set_defaults(handler=command_demo)

    check = commands.add_parser("check-input", help="describe a NetCDF input and check that SCORCH can read it")
    check.add_argument("file", help="path to a processed temperature NetCDF file")
    check.add_argument("--settings", metavar="FILE", help="settings file to check the calendar against")
    check.add_argument("--strict", action="store_true",
                       help="also check the calendar and grid against the (manuscript) settings")
    check.set_defaults(handler=command_check_input)

    settings = commands.add_parser("settings", help="print the effective settings as YAML")
    settings.add_argument("--settings", metavar="FILE", help="settings file to apply first")
    settings.add_argument("--write", metavar="FILE", help="write the settings to this file instead")
    settings.set_defaults(handler=command_settings)

    info = commands.add_parser("info", help="show versions of SCORCH and its main packages")
    info.set_defaults(handler=command_info)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (io.InputError, cfg.SettingsError, StageError) as error:
        print(f"\n{PROGRAM}: error: {error}", file=sys.stderr)
        if args.debug:
            traceback.print_exc()
        return 1
    except (ValueError, FileNotFoundError, ImportError) as error:
        print(f"\n{PROGRAM}: error: {error}", file=sys.stderr)
        if args.debug:
            traceback.print_exc()
        else:
            print("(run again with --debug for the full traceback)", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
