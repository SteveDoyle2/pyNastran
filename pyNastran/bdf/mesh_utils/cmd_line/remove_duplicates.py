"""
Removes cards listed in a JSON/CSV report from a BDF/punch file.

bdf remove_duplicates fem.blk duplicates.json --punch
bdf remove_duplicates fem.bdf duplicates.json -o fem_v2.bdf
"""

from __future__ import annotations
import csv
import json
import os
import sys
from typing import TYPE_CHECKING

from cpylog import SimpleLogger

import pyNastran
from .utils import add_argparse_arguments

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_remove_duplicates(argv=None, quiet: bool = False) -> None:
    """command line interface to ``bdf remove_duplicates``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse

    parent_parser = argparse.ArgumentParser(
        description="Remove cards listed in a JSON/CSV report from a BDF file.",
    )
    parent_parser.add_argument("remove_duplicates", type=str)
    parent_parser.add_argument("INPUT", help="path to input BDF/DAT/NAS/BLK file", type=str)
    parent_parser.add_argument("REPORT", help="path to duplicates JSON/CSV file", type=str)
    parent_parser.add_argument(
        "-o", "--output", help="path to output file (default: INPUT_v2.ext)", type=str
    )
    add_argparse_arguments(parent_parser, ["--punch", "--lax", "--allow_dup"])
    parent_parser.add_argument("-v", "--version", action="version", version=pyNastran.__version__)
    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename = args.INPUT
    report_filename = args.REPORT
    bdf_filename_out = args.output
    punch = args.punch
    duplicate_cards = args.allow_dup.split(",") if args.allow_dup else []

    if bdf_filename_out is None:
        base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f"{base}_v2{ext}"

    level = "debug" if not quiet else "warning"
    log = SimpleLogger(level=level, encoding="utf-8")

    remove_duplicates(
        bdf_filename,
        report_filename,
        bdf_filename_out=bdf_filename_out,
        punch=punch,
        duplicate_cards=duplicate_cards,
        log=log,
    )


def remove_duplicates(
    bdf_filename: str,
    report_filename: str,
    bdf_filename_out: str = "",
    punch: bool = False,
    duplicate_cards: list[str] | None = None,
    log: SimpleLogger | None = None,
) -> None:
    """
    Reads a BDF, removes cards listed in a JSON/CSV report,
    and writes a new output file.

    Parameters
    ----------
    bdf_filename : str
        path to the input BDF/punch file
    report_filename : str
        path to the duplicates JSON/CSV report
    bdf_filename_out : str; default=''
        path to the output file; empty appends '_v2'
    punch : bool; default=False
        is this a punch file
    duplicate_cards : list[str]; default=None
        cards to allow duplicates for
    log : SimpleLogger; default=None
        logging object

    """
    from pyNastran.bdf.bdf import BDF

    if log is None:
        log = SimpleLogger(level="warning", encoding="utf-8")

    if not bdf_filename_out:
        base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f"{base}_v2{ext}"

    remove_set = _load_remove_set(report_filename, log)
    if not remove_set:
        log.warning("no cards to remove")
        return

    model = BDF(log=log)
    if duplicate_cards:
        model.set_allow_duplicates(set(duplicate_cards))
    model._stop_on_duplicate_error = False
    model.read_bdf(
        bdf_filename,
        punch=punch,
        xref=False,
        save_file_structure=False,
    )

    n_removed = _remove_cards(model, remove_set, log)

    model.write_bdf(bdf_filename_out)
    log.info(f"removed {n_removed} cards")
    log.info(f"wrote {bdf_filename_out}")


OBJ_NAME_TO_DICT = {
    "node": "nodes",
    "element": "elements",
    "rigid_element": "rigid_elements",
    "property": "properties",
    "material": "materials",
    "thermal_material": "thermal_materials",
    "mass": "masses",
    "coord": "coords",
}


def _load_remove_set(report_filename: str, log: SimpleLogger) -> set[tuple[str, int]]:
    """
    Loads the report and returns a set of (obj_name, id) pairs to remove.

    Only entries with action != 'skip' are included.
    """
    _base, ext = os.path.splitext(report_filename)
    if ext.lower() == ".csv":
        return _load_csv_remove_set(report_filename, log)
    return _load_json_remove_set(report_filename, log)


def _load_json_remove_set(json_filename: str, log: SimpleLogger) -> set[tuple[str, int]]:
    """load card IDs to remove from a JSON report"""
    with open(json_filename, "r", encoding="utf-8") as f:
        report = json.load(f)

    remove_set: set[tuple[str, int]] = set()
    for dup in report.get("duplicates", []):
        action = dup.get("action", "skip")
        if action == "skip":
            continue
        obj_name = dup["obj_name"]
        card_id = dup["id"]
        remove_set.add((obj_name, card_id))
    log.info(f"loaded {len(remove_set)} cards to remove from {json_filename}")
    return remove_set


def _load_csv_remove_set(csv_filename: str, log: SimpleLogger) -> set[tuple[str, int]]:
    """load card IDs to remove from a CSV report"""
    remove_set: set[tuple[str, int]] = set()
    with open(csv_filename, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            action = row.get("action", "skip")
            if action == "skip":
                continue
            obj_name = row["obj_name"]
            card_id = int(row["id"])
            remove_set.add((obj_name, card_id))
    log.info(f"loaded {len(remove_set)} cards to remove from {csv_filename}")
    return remove_set


def _remove_cards(
    model: BDF,
    remove_set: set[tuple[str, int]],
    log: SimpleLogger,
) -> int:
    """delete the listed cards from the model dicts"""
    n_removed = 0
    n_missing = 0
    for obj_name, card_id in sorted(remove_set):
        dict_name = OBJ_NAME_TO_DICT.get(obj_name, "")
        if not dict_name:
            log.warning(f"unknown obj_name={obj_name!r} for id={card_id}, skipping")
            continue

        obj_dict = getattr(model, dict_name, None)
        if obj_dict is None or card_id not in obj_dict:
            n_missing += 1
            continue

        del obj_dict[card_id]
        n_removed += 1

    if n_missing:
        log.info(f"{n_missing} cards not found in this file (expected)")
    return n_removed
