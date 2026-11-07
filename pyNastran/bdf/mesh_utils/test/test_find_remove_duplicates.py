"""tests for bdf find_duplicates / bdf remove_duplicates"""

import json
import os
import unittest
from pathlib import Path

from cpylog import SimpleLogger

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.cmd_line.find_duplicates import (
    cmd_line_find_duplicates,
    find_duplicates,
)
from pyNastran.bdf.mesh_utils.cmd_line.remove_duplicates import (
    cmd_line_remove_duplicates,
    remove_duplicates,
)

TEST_DIR = Path(__file__).resolve().parent
LOG = SimpleLogger(level="warning", encoding="utf-8")


def _write_test_model(test_dir: Path) -> tuple[str, str, str]:
    """
    Creates a small test model with two include files that share
    duplicate GRID and PSHELL cards.

    Returns (main_bdf, include1, include2) paths.

    Duplicated cards (equivalent):
        GRID 100, 101, 102  (in both includes)
        PSHELL 10           (in both includes)

    Unique to include1:
        GRID 200, 201
        CQUAD4 1
        PSHELL 20

    Unique to include2:
        GRID 300, 301
        CQUAD4 2
        PSHELL 30
    """
    main_bdf = str(test_dir / "tmp_dup_main.bdf")
    inc1 = str(test_dir / "tmp_dup_inc1.blk")
    inc2 = str(test_dir / "tmp_dup_inc2.blk")

    with open(main_bdf, "w") as f:
        f.write("SOL 101\n")
        f.write("CEND\n")
        f.write("BEGIN BULK\n")
        f.write(f"INCLUDE '{os.path.basename(inc1)}'\n")
        f.write(f"INCLUDE '{os.path.basename(inc2)}'\n")
        f.write("ENDDATA\n")

    # Use BDF API to write properly formatted cards
    m1 = BDF(log=SimpleLogger(level='error'))
    m1.add_grid(100, [0.0, 0.0, 0.0])
    m1.add_grid(101, [1.0, 0.0, 0.0])
    m1.add_grid(102, [1.0, 1.0, 0.0])
    m1.add_grid(200, [2.0, 0.0, 0.0])
    m1.add_grid(201, [2.0, 1.0, 0.0])
    m1.add_pshell(10, mid1=1, t=0.1)
    m1.add_pshell(20, mid1=1, t=0.2)
    m1.add_cquad4(1, 10, [100, 101, 102, 200])
    m1.add_mat1(1, 1.0e7, None, 0.3)
    m1.write_bdf(inc1, enddata=False)

    m2 = BDF(log=SimpleLogger(level='error'))
    m2.add_grid(100, [0.0, 0.0, 0.0])
    m2.add_grid(101, [1.0, 0.0, 0.0])
    m2.add_grid(102, [1.0, 1.0, 0.0])
    m2.add_grid(300, [3.0, 0.0, 0.0])
    m2.add_grid(301, [3.0, 1.0, 0.0])
    m2.add_pshell(10, mid1=1, t=0.1)
    m2.add_pshell(30, mid1=1, t=0.3)
    m2.add_cquad4(2, 10, [100, 101, 102, 300])
    m2.write_bdf(inc2, enddata=False)

    return main_bdf, inc1, inc2


def _cleanup(paths: list[str]) -> None:
    for p in paths:
        if os.path.exists(p):
            os.remove(p)


class TestFindDuplicates(unittest.TestCase):
    """
    Tests for bdf find_duplicates.

    Tolerances: N/A (integer ID comparisons).
    Tests: equivalent duplicate detection, JSON/CSV report content,
           non-equivalent duplicates marked as skip.
    """

    def test_find_duplicates_api(self) -> None:
        """find_duplicates() detects equivalent GRID and PSHELL cards"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_dup_report.json")
        csv_file = str(TEST_DIR / "tmp_dup_report.csv")
        try:
            find_duplicates(
                main_bdf,
                json_file,
                write_csv=True,
                log=LOG,
            )
            assert os.path.exists(json_file), json_file
            assert os.path.exists(csv_file), csv_file

            with open(json_file, "r") as f:
                report = json.load(f)

            dups = report["duplicates"]
            # 3 GRID + 1 PSHELL = 4 equivalent duplicates
            assert len(dups) == 4, f"expected 4, got {len(dups)}: {dups}"

            for d in dups:
                assert d["equivalent"] is True, d
                assert d["action"] == "remove_from_file2", d

            # check IDs
            dup_ids = {(d["obj_name"], d["id"]) for d in dups}
            expected = {
                ("node", 100),
                ("node", 101),
                ("node", 102),
                ("property", 10),
            }
            assert dup_ids == expected, f"{dup_ids} != {expected}"

            # file_kept should be inc1, file_discarded should be inc2
            for d in dups:
                assert os.path.basename(d["file_kept"]) == "tmp_dup_inc1.blk", d
                assert os.path.basename(d["file_discarded"]) == "tmp_dup_inc2.blk", d

            # files dict should list all 3 files
            assert len(report["files"]) == 3, report["files"]

        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, csv_file])

    def test_find_duplicates_cmdline(self) -> None:
        """cmd_line_find_duplicates() works via argv"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_dup_cmdline.json")
        try:
            argv = ["bdf", "find_duplicates", main_bdf, "-o", json_file]
            cmd_line_find_duplicates(argv=argv, quiet=True)
            assert os.path.exists(json_file), json_file

            with open(json_file, "r") as f:
                report = json.load(f)
            assert len(report["duplicates"]) == 4
        finally:
            _cleanup([main_bdf, inc1, inc2, json_file])

    def test_find_duplicates_no_duplicates(self) -> None:
        """find_duplicates on a model with no duplicates returns empty list"""
        bdf_file = str(TEST_DIR / "tmp_no_dup.bdf")
        json_file = str(TEST_DIR / "tmp_no_dup_report.json")
        try:
            m = BDF(log=SimpleLogger(level="error"))
            m.add_grid(1, [0.0, 0.0, 0.0])
            m.add_grid(2, [1.0, 0.0, 0.0])
            m.write_bdf(bdf_file, enddata=False)
            find_duplicates(bdf_file, json_file, punch=True, log=LOG)

            with open(json_file, "r") as f:
                report = json.load(f)
            assert len(report["duplicates"]) == 0
        finally:
            _cleanup([bdf_file, json_file])

    def test_find_duplicates_punch(self) -> None:
        """find_duplicates works on a single punch file with no includes"""
        bdf_file = str(TEST_DIR / "tmp_punch_dup.bdf")
        json_file = str(TEST_DIR / "tmp_punch_dup_report.json")
        try:
            m = BDF(log=SimpleLogger(level="error"))
            m.add_grid(1, [0.0, 0.0, 0.0])
            m.add_grid(2, [1.0, 0.0, 0.0])
            m.write_bdf(bdf_file, enddata=False)
            # append a duplicate GRID 1
            with open(bdf_file, "a") as f:
                f.write("GRID           1              0.      0.      0.\n")
            find_duplicates(bdf_file, json_file, punch=True, log=LOG)

            with open(json_file, "r") as f:
                report = json.load(f)
            dups = report["duplicates"]
            assert len(dups) == 1, dups
            assert dups[0]["id"] == 1
            assert dups[0]["card_type"] == "GRID"
            assert dups[0]["equivalent"] is True
        finally:
            _cleanup([bdf_file, json_file])


class TestRemoveDuplicates(unittest.TestCase):
    """
    Tests for bdf remove_duplicates.

    Tolerances: N/A (integer ID comparisons).
    Tests: card deletion from include file, skip of missing cards,
           CSV report input, cmdline interface.
    """

    def test_remove_duplicates_api(self) -> None:
        """remove_duplicates() removes listed cards from a file"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_rm_report.json")
        inc2_v2 = str(TEST_DIR / "tmp_dup_inc2_v2.blk")
        try:
            find_duplicates(main_bdf, json_file, log=LOG)

            remove_duplicates(
                inc2,
                json_file,
                bdf_filename_out=inc2_v2,
                punch=True,
                log=LOG,
            )
            assert os.path.exists(inc2_v2), inc2_v2

            # Read the v2 include file
            model_v2 = BDF(log=LOG)
            model_v2.read_bdf(inc2_v2, punch=True, xref=False)

            # 3 duplicate GRIDs removed, 2 unique remain
            assert len(model_v2.nodes) == 2, f"nodes: {list(model_v2.nodes.keys())}"
            assert 300 in model_v2.nodes
            assert 301 in model_v2.nodes
            assert 100 not in model_v2.nodes

            # 1 duplicate PSHELL removed, 1 unique remains
            assert len(model_v2.properties) == 1, f"props: {list(model_v2.properties.keys())}"
            assert 30 in model_v2.properties
            assert 10 not in model_v2.properties

            # Element should still be there
            assert 2 in model_v2.elements

        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, inc2_v2])

    def test_remove_duplicates_cmdline(self) -> None:
        """cmd_line_remove_duplicates() works via argv"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_rm_cmd_report.json")
        inc2_v2 = str(TEST_DIR / "tmp_dup_inc2_v2.blk")
        try:
            find_duplicates(main_bdf, json_file, log=LOG)

            argv = [
                "bdf",
                "remove_duplicates",
                inc2,
                json_file,
                "--punch",
                "-o",
                inc2_v2,
            ]
            cmd_line_remove_duplicates(argv=argv, quiet=True)

            model_v2 = BDF(log=LOG)
            model_v2.read_bdf(inc2_v2, punch=True, xref=False)
            assert len(model_v2.nodes) == 2
            assert len(model_v2.properties) == 1
        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, inc2_v2])

    def test_remove_duplicates_csv(self) -> None:
        """remove_duplicates() accepts a CSV report"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_rm_csv_report.json")
        csv_file = str(TEST_DIR / "tmp_rm_csv_report.csv")
        inc2_v2 = str(TEST_DIR / "tmp_dup_inc2_v2.blk")
        try:
            find_duplicates(
                main_bdf,
                json_file,
                write_csv=True,
                log=LOG,
            )

            remove_duplicates(
                inc2,
                csv_file,
                bdf_filename_out=inc2_v2,
                punch=True,
                log=LOG,
            )

            model_v2 = BDF(log=LOG)
            model_v2.read_bdf(inc2_v2, punch=True, xref=False)
            assert len(model_v2.nodes) == 2
            assert len(model_v2.properties) == 1
        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, csv_file, inc2_v2])

    def test_remove_skips_missing_cards(self) -> None:
        """remove_duplicates on inc1 skips cards not in that file"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_rm_skip_report.json")
        inc1_v2 = str(TEST_DIR / "tmp_dup_inc1_v2.blk")
        try:
            find_duplicates(main_bdf, json_file, log=LOG)

            # Run on inc1 — the duplicates are also in inc1, so they
            # will be removed. But unique inc2 cards are not in inc1,
            # so those are skipped.
            remove_duplicates(
                inc1,
                json_file,
                bdf_filename_out=inc1_v2,
                punch=True,
                log=LOG,
            )

            model_orig = BDF(log=LOG)
            model_orig.read_bdf(inc1, punch=True, xref=False)

            model_v2 = BDF(log=LOG)
            model_v2.read_bdf(inc1_v2, punch=True, xref=False)

            # duplicate GRIDs 100,101,102 removed; unique 200,201 remain
            assert len(model_v2.nodes) == 2, f"nodes: {list(model_v2.nodes.keys())}"
            assert 200 in model_v2.nodes
            assert 201 in model_v2.nodes

            # duplicate PSHELL 10 removed; unique PSHELL 20 remains
            assert len(model_v2.properties) == 1
            assert 20 in model_v2.properties

        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, inc1_v2])

    def test_remove_no_actions(self) -> None:
        """remove_duplicates with all-skip report does nothing"""
        bdf_file = str(TEST_DIR / "tmp_no_action.bdf")
        json_file = str(TEST_DIR / "tmp_no_action_report.json")
        try:
            m = BDF(log=SimpleLogger(level="error"))
            m.add_grid(1, [0.0, 0.0, 0.0])
            m.write_bdf(bdf_file, enddata=False)

            report = {
                "duplicates": [
                    {"obj_name": "node", "id": 999, "action": "skip"},
                ],
            }
            with open(json_file, "w") as f:
                json.dump(report, f)

            # should not crash, should not write output
            remove_duplicates(
                bdf_file,
                json_file,
                punch=True,
                log=LOG,
            )
            # default v2 file should not be created since no actions
            v2_file = str(TEST_DIR / "tmp_no_action_v2.bdf")
            assert not os.path.exists(v2_file)
        finally:
            _cleanup([bdf_file, json_file])


class TestRoundTrip(unittest.TestCase):
    """
    End-to-end: find_duplicates -> remove_duplicates -> re-read shows
    zero duplicates.

    Tolerances: N/A.
    Tests: full pipeline roundtrip.
    """

    def test_roundtrip(self) -> None:
        """find → remove → re-read has zero duplicates"""
        main_bdf, inc1, inc2 = _write_test_model(TEST_DIR)
        json_file = str(TEST_DIR / "tmp_rt_report.json")
        inc2_v2 = str(TEST_DIR / "tmp_dup_inc2_v2.blk")
        main_v2 = str(TEST_DIR / "tmp_dup_main_v2.bdf")
        try:
            find_duplicates(main_bdf, json_file, log=LOG)
            remove_duplicates(
                inc2,
                json_file,
                bdf_filename_out=inc2_v2,
                punch=True,
                log=LOG,
            )

            # Write a new main file pointing to inc1 and inc2_v2
            with open(main_v2, "w") as f:
                f.write("SOL 101\n")
                f.write("CEND\n")
                f.write("BEGIN BULK\n")
                f.write(f"INCLUDE '{os.path.basename(inc1)}'\n")
                f.write(f"INCLUDE '{os.path.basename(inc2_v2)}'\n")
                f.write("ENDDATA\n")

            # Re-read and check for zero duplicates
            model = BDF(log=LOG)
            model._stop_on_duplicate_error = False
            model.read_bdf(main_v2, xref=False, save_file_structure=True)
            assert len(model._duplicate_equivalents) == 0, (
                f"still have {len(model._duplicate_equivalents)} duplicates"
            )

            # Total cards: inc1 has 5 grids, inc2_v2 has 2 grids = 7
            assert len(model.nodes) == 7, f"nodes: {list(model.nodes.keys())}"
            # inc1: PSHELL 10, 20; inc2_v2: PSHELL 30 = 3
            assert len(model.properties) == 3
            # 2 elements total
            assert len(model.elements) == 2

        finally:
            _cleanup([main_bdf, inc1, inc2, json_file, inc2_v2, main_v2])


if __name__ == "__main__":
    unittest.main()
