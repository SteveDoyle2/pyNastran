"""
Finds duplicate nodes, elements, properties, etc. across include files
and writes a JSON report.

bdf find_duplicates fem.bdf
bdf find_duplicates fem.bdf -o duplicates.json
"""

from __future__ import annotations
import csv
import json
import os
import sys
from collections import defaultdict
from typing import TYPE_CHECKING

from cpylog import SimpleLogger

import pyNastran
from .utils import add_argparse_arguments

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_find_duplicates(argv=None, quiet: bool = False) -> str:
    """command line interface to ``bdf find_duplicates``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse

    parent_parser = argparse.ArgumentParser(
        description="Find duplicate cards across include files and write a JSON report.",
    )
    parent_parser.add_argument("find_duplicates", type=str)
    parent_parser.add_argument("INPUT", help="path to input BDF/DAT/NAS file", type=str)
    parent_parser.add_argument("-o", "--output", help="path to output JSON file", type=str)
    parent_parser.add_argument("--csv", help="also write a CSV report", action="store_true")
    add_argparse_arguments(parent_parser, ["--punch", "--lax", "--allow_dup"])
    parent_parser.add_argument("-v", "--version", action="version", version=pyNastran.__version__)
    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename = args.INPUT
    json_filename = args.output
    punch = args.punch
    is_strict_card_parser = not args.lax
    duplicate_cards = args.allow_dup.split(",") if args.allow_dup else []

    if json_filename is None:
        base, _ext = os.path.splitext(bdf_filename)
        json_filename = f"{base}.duplicates.json"

    level = "debug" if not quiet else "warning"
    log = SimpleLogger(level=level, encoding="utf-8")

    write_csv = args.csv
    json_filename = find_duplicates(
        bdf_filename,
        json_filename,
        punch=punch,
        is_strict_card_parser=is_strict_card_parser,
        duplicate_cards=duplicate_cards,
        write_csv=write_csv,
        log=log,
    )
    return json_filename


def find_duplicates(
    bdf_filename: str,
    json_filename: str,
    punch: bool = False,
    is_strict_card_parser: bool = True,
    duplicate_cards: list[str] | None = None,
    write_csv: bool = False,
    log: SimpleLogger | None = None,
) -> str:
    """
    Reads a BDF with save_file_structure=True, collects all
    equivalent duplicate cards, and writes a JSON report.

    Parameters
    ----------
    bdf_filename : str
        path to the input BDF
    json_filename : str
        path to the output JSON report
    punch : bool; default=False
        is this a punch file
    is_strict_card_parser : bool; default=True
        strict card parsing
    duplicate_cards : list[str]; default=None
        cards to allow duplicates for
    log : SimpleLogger; default=None
        logging object

    Returns
    -------
    json_filename : str
        path to the written JSON report

    """
    from pyNastran.bdf.bdf import BDF

    if log is None:
        log = SimpleLogger(level="warning", encoding="utf-8")

    model = BDF(log=log)
    if duplicate_cards:
        model.set_allow_duplicates(set(duplicate_cards))
    model._stop_on_duplicate_error = False
    model.read_bdf(
        bdf_filename,
        punch=punch,
        xref=False,
        save_file_structure=True,
    )

    duplicates = _collect_duplicates(model)

    report = {
        "bdf_filename": os.path.abspath(bdf_filename),
        "files": {str(i): fname for i, fname in enumerate(model.active_filenames)},
        "duplicates": duplicates,
    }

    with open(json_filename, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    n = len(duplicates)
    log.info(f"found {n} duplicate cards")
    log.info(f"wrote {json_filename}")

    if write_csv:
        csv_filename = os.path.splitext(json_filename)[0] + ".csv"
        _write_csv_report(duplicates, csv_filename, log)

    _print_summary(duplicates, model.active_filenames, log)
    return json_filename


def _collect_duplicates(model: BDF) -> list[dict]:
    """
    Collects duplicate info from the model into a list of dicts.

    Each entry has:
        obj_name : str  ('node', 'element', 'property', ...)
        id : int
        card_type : str  ('GRID', 'CQUAD4', 'PSHELL', ...)
        equivalent : bool
        ifile_kept : int
        ifile_discarded : int
        file_kept : str
        file_discarded : str
        action : str  ('remove_from_file2')

    """
    duplicates = []
    filenames = model.active_filenames

    nfiles = len(filenames)
    for obj_name, card_id, card_type, ifile_kept, ifile_discarded in model._duplicate_equivalents:
        file_kept = filenames[ifile_kept] if ifile_kept < nfiles else f"ifile={ifile_kept}"
        file_discarded = (
            filenames[ifile_discarded] if ifile_discarded < nfiles else f"ifile={ifile_discarded}"
        )
        duplicates.append(
            {
                "obj_name": obj_name,
                "id": card_id,
                "card_type": card_type,
                "equivalent": True,
                "ifile_kept": ifile_kept,
                "ifile_discarded": ifile_discarded,
                "file_kept": file_kept,
                "file_discarded": file_discarded,
                "action": "remove_from_file2",
            }
        )

    for category, dup_list in model._duplicate.items():
        for obj in dup_list:
            card_id = _get_card_id(obj)
            ifile_discarded = obj.ifile if hasattr(obj, "ifile") else -1
            file_discarded = ""
            if ifile_discarded >= 0 and ifile_discarded < len(filenames):
                file_discarded = filenames[ifile_discarded]

            obj_name = category.rstrip("s")
            if category == "properties":
                obj_name = "property"
            elif category == "masses":
                obj_name = "mass"

            duplicates.append(
                {
                    "obj_name": obj_name,
                    "id": card_id,
                    "card_type": obj.type,
                    "equivalent": False,
                    "ifile_kept": -1,
                    "ifile_discarded": ifile_discarded,
                    "file_kept": "",
                    "file_discarded": file_discarded,
                    "action": "skip",
                }
            )
    return duplicates


def _get_card_id(obj) -> int:
    """get the ID from a card object"""
    for attr in ("nid", "eid", "pid", "mid", "cid"):
        if hasattr(obj, attr):
            return getattr(obj, attr)
    return -1


def _print_summary(duplicates: list[dict], filenames: list[str], log: SimpleLogger) -> None:
    """print a summary of duplicates by file and type"""
    by_file: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    for dup in duplicates:
        if dup["action"] == "skip":
            continue
        fname = dup["file_discarded"]
        by_file[fname][dup["card_type"]].append(dup["id"])

    if not by_file:
        log.info("no actionable duplicates found")
        return

    log.info("--- duplicate summary (to remove) ---")
    for fname, type_ids in sorted(by_file.items()):
        log.info(f"  {fname}:")
        for card_type, ids in sorted(type_ids.items()):
            if len(ids) <= 10:
                id_str = ", ".join(str(i) for i in sorted(ids))
            else:
                id_str = f"{len(ids)} cards"
            log.info(f"    {card_type}: {id_str}")


def _write_csv_report(duplicates: list[dict], csv_filename: str, log: SimpleLogger) -> None:
    """write the duplicates report as CSV"""
    fieldnames = [
        "obj_name",
        "id",
        "card_type",
        "equivalent",
        "ifile_kept",
        "ifile_discarded",
        "file_kept",
        "file_discarded",
        "action",
    ]
    with open(csv_filename, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for dup in duplicates:
            writer.writerow(dup)
    log.info(f"wrote {csv_filename}")
