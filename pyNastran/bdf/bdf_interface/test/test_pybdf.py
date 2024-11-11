# encoding: utf8
"""tests PyBDF"""
# pylint: disable=W0212
import os
import sys
import unittest
from io import StringIO
from cpylog import get_logger

import numpy as np
from pyNastran.bdf.bdf_interface.pybdf import (
    BDFInputPy, _show_bad_file, _lines_to_decks,
    MissingDeckSections,
    lines_to_decks2,
)


class TestPyBDF(unittest.TestCase):
    """tests PyBDF"""

    def test_pybdf_open_file_checks(self):
        """tests _open_file_checks"""
        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           log=None, debug=False)

        op2_filename = 'spike.op2'
        if os.path.exists(op2_filename):  # pragma: no cover
            os.remove(op2_filename)

        bdf_filename = 'spike.bdf'
        if os.path.exists(bdf_filename):  # pragma: no cover
            os.remove(bdf_filename)

        #------------------------------------------------------
        with self.assertRaises(IOError):
            pybdf._open_file_checks(op2_filename, basename=False)

        with open(op2_filename, 'w'):
            pass

        with self.assertRaises(IOError):
            pybdf._open_file_checks(op2_filename, basename=False)

        with self.assertRaises(IOError):
            pybdf._validate_open_file(op2_filename, check=True)

        with self.assertRaises(RuntimeError):
            _show_bad_file(pybdf, op2_filename, encoding, nlines_previous=10)

        #------------------------------------------------------
        with self.assertRaises(IOError):
            pybdf._validate_open_file(op2_filename, check=True)
        with open(bdf_filename, 'w', encoding='utf8') as bdf_file:
            bdf_file.write('helló wörld from two\n')

        ## TODO: why doesn't this fail with a unicode error?
        with self.assertRaises(RuntimeError):
            _show_bad_file(pybdf, bdf_filename, encoding='latin1', nlines_previous=10)

        #------------------------------------------------------
        bdf_dir = 'pybdf_dir'
        if not os.path.exists(bdf_dir):
            os.makedirs(bdf_dir)
        with self.assertRaises(IOError):
            pybdf._open_file_checks(bdf_dir, basename=False)

        with self.assertRaises(IOError):
            pybdf._validate_open_file(bdf_dir, check=True)

        os.remove('spike.op2')
        os.remove('spike.bdf')
        os.rmdir(bdf_dir)

    def test_get_lines_1(self):
        """tests the basic deck sections, with a GRID and a skipped POST"""
        punch = False
        lines = [
            'CEND',
            'BEGIN BULK',
            'GRID,1',
            'ENDDATA',
            'POST',
        ]
        ilines = np.array(
            [[0, 1],
             [0, 2],
             [0, 3],
             [0, 4],
             [0, 5],
             ])
        log = get_logger(log=None, level='debug', encoding='utf-8')
        out = _lines_to_decks(lines, ilines, punch, log,
                              keep_enddata=False, consider_superelements=False)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines, superelement_lines, superelement_ilines) = out
        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 0, case_control_lines
        assert len(bulk_data_lines) == 1, bulk_data_lines  #  TODO: make this 2
        #assert len(bulk_data_ilines) == 0, bulk_data_ilines
        assert len(superelement_lines) == 0, superelement_lines
        #assert len(superelement_ilines) == 0, superelement_ilines
        #for line in bulk_data_ilines:
            #print(line)

    def test_get_lines_2(self):
        """tests system control lines"""
        with open('junk.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1\n')
            bdf_file.write('ENDDATA\n')

        bdf_filename = StringIO()
        bdf_filename.write(
            'ASSIGN FEM="junk.bdf"\n'
            'ASSIGN MATRIX=junk.mat\n'
            'ASSIGN OUTPUT4=junk.op4\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,2'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           consider_superelements=False, log=None, debug=False)
        bulk_data_lines = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)[3]
        #(unused_system_lines,
         #unused_executive_control_lines,
         #unused_case_control_lines,
         #bulk_data_lines,
         #unused_bulk_data_ilines) = out
        #print('bulk_data_linesA =', bulk_data_lines)
        assert bulk_data_lines == ['GRID,1', 'GRID,2'], bulk_data_lines

        #-----------------------------------------------------------------------
        bdf_filename = StringIO()
        bdf_filename.write(
            'ASSIGN FEM=junk.f06\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,2\n'
            #'ENDDATA'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           consider_superelements=False, log=None, debug=False)
        bulk_data_lines = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)[3]
        #print('bulk_data_linesB =', bulk_data_lines)

        #-----------------------------------------------------------------------
        bdf_filename = StringIO()
        bdf_filename.write(
            "ASSIGN FEM='junk.f06'\n"
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,3\n'
            #'ENDDATA'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           consider_superelements=False, log=None, debug=False)
        bulk_data_lines = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)[3]
        #print('bulk_data_linesC =', bulk_data_lines)
        os.remove('junk.bdf')

    def test_no_punch2(self):
        """tests not definng punch"""
        bdf_filename = StringIO()
        bdf_filename.write(
            'GRID,1,,0.,0.,0.\n'
            'GRID.2,,1.,0.,0.\n'
            'GRID,3,,1.,1.,0.\n'
            'GRID,4,,0.,1.,0.\n'
            'CQUAD4,1,2,3,4,5',
            #'ENDDATA'
        )
        bdf_filename.seek(0)
        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        pybdf.use_new_parser = True
        decks = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)

        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 0, executive_control_lines
        assert len(case_control_lines) == 0, case_control_lines
        assert len(bulk_data_lines) == 5, bulk_data_lines  #  TODO: make this 2
        #assert len(bulk_data_ilines) == 0, bulk_data_ilines
        assert len(additional_deck_lines) == 0, additional_deck_lines

    def test_no_punch1(self):
        """tests not definng punch"""
        bdf_filename = StringIO()
        bdf_filename.write(
            'GRID,1,,0.,0.,0.\n'
            'GRID.2,,1.,0.,0.\n'
            'GRID,3,,1.,1.,0.\n'
            'GRID,4,,0.,1.,0.\n'
            'CQUAD4,1,2,3,4,5',
            #'ENDDATA'
        )
        bdf_filename.seek(0)
        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        with self.assertRaises(MissingDeckSections):
            bulk_data_lines = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)[3]

        #pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           #consider_superelements=False, log=None, debug=False)
        bdf_filename.seek(0)
        decks = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)

        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 0, executive_control_lines
        assert len(case_control_lines) == 0, case_control_lines
        assert len(bulk_data_lines) == 5, bulk_data_lines  #  TODO: make this 2
        #assert len(bulk_data_ilines) == 0, bulk_data_ilines
        assert len(additional_deck_lines) == 0, additional_deck_lines
        #print(bulk_data_lines)
        # -----------------------------------
        bdf_filename = StringIO()
        bdf_filename.write(
            'CEND\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID.2,,1.,0.,0.\n'
            'GRID,3,,1.,1.,0.\n'
            'GRID,4,,0.,1.,0.\n'
            'CQUAD4,1,2,3,4,5',
            #'ENDDATA'
        )
        bdf_filename.seek(0)
        decks = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)

        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 0, case_control_lines
        assert len(bulk_data_lines) == 5, bulk_data_lines  #  TODO: make this 2
        #assert len(bulk_data_ilines) == 0, bulk_data_ilines
        assert len(additional_deck_lines) == 0, additional_deck_lines

    def test_unicode_errors1(self):
        """tests some error handling"""
        bdf_filename = 'unicode.bdf'
        _write_unicode_deck(bdf_filename)

        read_includes = True
        dumplines = True
        encoding = 'utf8'
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        #if sys.platform == 'win32':
            #with self.assertRaises(RuntimeError):
                #pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)
        #else:
        pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)

        #with self.assertRaises(RuntimeError):
        #with self.assertRaises(UnicodeDecodeError):
        #unused_bulk_data_lines1 = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)[3]

    def test_unicode_errors2(self):
        """tests some error handling"""
        bdf_filename = 'unicode.bdf'
        _write_unicode_deck(bdf_filename)

        read_includes = True
        dumplines = True
        #unused_bulk_data_lines1 = ...

        encoding = 'latin1'
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        unused_bulk_data_lines2 = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)[3]
        os.remove(bdf_filename)

        #--------------
        bdf_filename = 'inc.inc'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(
                'GRID,3,,1.,1.,0.\n'
                'GRID,4,,0.,1.,0.\n'
                "INCLUDE 'main.bdf'\n"
            )

        bdf_filename = 'main.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(
                'GRID,1,,0.,0.,0.\n'
                "INCLUDE 'inc.inc'\n"
            )

        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        with self.assertRaises(RuntimeError):
            unused_bulk_data_lines3 = pybdf.get_lines(bdf_filename, punch=True, make_ilines=True)
        os.remove('main.bdf')
        os.remove('inc.inc')

    def test_line_to_decks_nominal(self):
        read_includes = True
        dumplines = False
        encoding = None
        #consider_superelements = True
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=True, log=None, debug=False)
        main_lines = [
            'SOL 101',                              # executive
            'CEND',                                 # executive
            'SUBCASE 1',                            # case
            '  LOAD = 1',                           # case
            'BEGIN BULK',                           # tag
            'GRID,1',                               # bulk
            'ENDDATA',                              # bulk
        ]
        bdf_filename = StringIO()
        bdf_filename.write('\n'.join(main_lines))
        bdf_filename.seek(0)

        log = get_logger(log=None, level='debug', encoding='utf-8')
        iline = None
        punch = None
        decks = lines_to_decks2(main_lines, iline, punch, log)
        out = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         additional_deck_lines) = out

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 2, executive_control_lines
        assert len(case_control_lines) == 3, case_control_lines
        assert len(bulk_data_lines) == 2, bulk_data_lines
        assert len(additional_deck_lines) == 0, additional_deck_lines

        decks = lines_to_decks2(main_lines, iline, punch, log)
        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 2, case_control_lines
        assert len(bulk_data_lines) == 2, bulk_data_lines
        assert len(additional_deck_lines) == 0, additional_deck_lines

    def test_line_to_decks_massid(self):
        #read_includes = True
        #dumplines = False
        #encoding = None
        #consider_superelements = True
        #pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           #consider_superelements=True, log=None, debug=False)
        main_lines = [
            'SOL 101',                         # executive
            'CEND',                            # tag
            'SUBCASE 1',                       # case
            '  LOAD = 1',                      # case
            "begin massid=1 label='cat dog' ", # tag
            'CONM2,1',                         # MASSID=1
            "begin massid=2 label='cat dog' ", # tag
            'CONM2,2',                         # MASSID=2
            'GRID,1',
            'ENDDATA',
        ]
        log = get_logger(log=None, level='debug', encoding='utf-8')
        iline = None
        punch = None
        decks = lines_to_decks2(main_lines, iline, punch, log)
        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks
        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 2, case_control_lines
        assert len(bulk_data_lines) == 0, bulk_data_lines
        assert len(additional_deck_lines) == 2, additional_deck_lines
        assert len(additional_deck_lines[('MASSID', 1, 'CAT DOG')]) == 1
        assert len(additional_deck_lines[('MASSID', 2, 'CAT DOG')]) == 3

    def test_line_to_decks_super_1(self):
        log = get_logger(log=None, level='debug', encoding='utf-8')
        read_includes = True
        dumplines = False
        encoding = None
        #consider_superelements = True
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=True, log=None, debug=False)
        main_lines = [
            # system
            'SOL 101',       # executive
            'CEND',          # tag
            'SUBCASE 1',     # case
            '  LOAD = 1',    # case
            "begin super=1", # tag
            'CONM2,1',       # SUPER=1
            "begin super=2", # tag
            'CONM2,2',       # SUPER=2
            'GRID,1',        # SUPER=2
            'GRID,2',        # SUPER=2
            'begin bulk',    # tag
            'GRID,1',        # bulk
            'ENDDATA',       # bulk
        ]
        bdf_filename = StringIO()
        bdf_filename.write('\n'.join(main_lines))
        bdf_filename.seek(0)

        out = pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         additional_deck_lines) = out

        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 2, executive_control_lines
        assert len(case_control_lines) == 3, case_control_lines
        assert len(bulk_data_lines) == 2, bulk_data_lines
        assert len(additional_deck_lines) == 2, additional_deck_lines
        ntrash_lines = 2
        assert len(additional_deck_lines[('SUPER', 1, '')]) == 1 + ntrash_lines, additional_deck_lines[('SUPER', 1, '')]
        assert len(additional_deck_lines[('SUPER', 2, '')]) == 3 + ntrash_lines, additional_deck_lines[('SUPER', 2, '')]

        #-----------------------------------------------
        iline = None
        punch = None
        decks = lines_to_decks2(main_lines, iline, punch, log)
        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks
        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 2, case_control_lines
        assert len(bulk_data_lines) == 2, bulk_data_lines
        assert len(additional_deck_lines) == 2, additional_deck_lines
        assert len(additional_deck_lines[('SUPER', 1, '')]) == 1
        assert len(additional_deck_lines[('SUPER', 2, '')]) == 3

    def test_line_to_decks_super_2(self):
        #read_includes = True
        #dumplines = False
        #encoding = None
        #consider_superelements = True
        #pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           #consider_superelements=True, log=None, debug=False)
        main_lines = [
            # system
            'SOL 101',       # executive
            'CEND',          # tag
            'SUBCASE 1',     # case
            '  LOAD = 1',    # case
            "begin super=1", # tag
            'CONM2,1',       # SUPER=1
            "begin super",   # tag
            'CONM2,2',       # SUPER=0
            'GRID,1',        # SUPER=0
            'ENDDATA',       # SUPER=0
        ]
        log = get_logger(log=None, level='debug', encoding='utf-8')
        iline = None
        punch = None
        decks = lines_to_decks2(main_lines, iline, punch, log)
        system_lines, executive_control_lines, case_control_lines, \
            bulk_data_lines, bulk_data_ilines, \
            additional_deck_lines = decks
        assert len(system_lines) == 0, system_lines
        assert len(executive_control_lines) == 1, executive_control_lines
        assert len(case_control_lines) == 2, case_control_lines
        assert len(bulk_data_lines) == 0, bulk_data_lines
        assert len(additional_deck_lines) == 2, additional_deck_lines
        assert len(additional_deck_lines[('SUPER', 1, '')]) == 1
        assert len(additional_deck_lines[('SUPER', 0, '')]) == 3

    def test_lines_to_decks_super_mass(self):
        main_lines = [
            # system
            'SOL 101',       # executive
            'CEND',          # tag
            'SUBCASE 1',     # case
            '  LOAD = 1',    # case
            #"begin super=1", # tag
            'CONM2,1',       # SUPER=1
            #"BEGIN MASSID=300 LABEL = 'cat dog'",   # tag
            'BEGIN   SUPER=7 MASSID=300 LABEL=dog',   # tag
            'CONM2,2',       # SUPER=0
            'GRID,1',        # SUPER=0
            'ENDDATA',       # SUPER=0
        ]
        log = get_logger(log=None, level='debug', encoding='utf-8')
        iline = None
        punch = None
        with self.assertRaises(AssertionError):
            decks = lines_to_decks2(main_lines, iline, punch, log)
            #system_lines, executive_control_lines, case_control_lines, \
                #bulk_data_lines, bulk_data_ilines, \
                #additional_deck_lines = decks
        #assert len(system_lines) == 0, system_lines
        #assert len(executive_control_lines) == 1, executive_control_lines
        #assert len(case_control_lines) == 2, case_control_lines
        #assert len(bulk_data_lines) == 0, bulk_data_lines
        #assert len(additional_deck_lines) == 2, additional_deck_lines
        #assert len(additional_deck_lines[('SUPER', 1, '')]) == 1
        #assert len(additional_deck_lines[('SUPER', 0, '')]) == 3

    def test_lines_to_decks_module(self):
        main_lines = [
            # system
            'SOL 101',       # executive
            'CEND',          # tag
            'SUBCASE 1',     # case
            '  LOAD = 1',    # case
            #"begin super=1", # tag
            'CONM2,1',       # SUPER=1
            #"BEGIN MASSID=300 LABEL = 'cat dog'",   # tag
            "BEGIN MODULE=1 APPEND LABEL='MODULE1'",
            #'BEGIN   SUPER=7 MASSID=300 LABEL=dog',   # tag
            'CONM2,2',       # SUPER=0
            'GRID,1',        # SUPER=0
            'ENDDATA',       # SUPER=0
        ]
        log = get_logger(log=None, level='debug', encoding='utf-8')
        iline = None
        punch = None
        with self.assertRaises(AssertionError):
            decks = lines_to_decks2(main_lines, iline, punch, log)
            #system_lines, executive_control_lines, case_control_lines, \
                #bulk_data_lines, bulk_data_ilines, \
                #additional_deck_lines = decks
        #assert len(system_lines) == 0, system_lines
        #assert len(executive_control_lines) == 1, executive_control_lines
        #assert len(case_control_lines) == 2, case_control_lines
        #assert len(bulk_data_lines) == 0, bulk_data_lines
        #assert len(additional_deck_lines) == 2, additional_deck_lines
        #assert len(additional_deck_lines[('SUPER', 1, '')]) == 1
        #assert len(additional_deck_lines[('SUPER', 0, '')]) == 3

def _write_unicode_deck(bdf_filename: str) -> None:
    with open(bdf_filename, 'w', encoding='utf8') as bdf_file:
        bdf_file.write(
            'CEND\n'
            'SUBCASE 1\n'
            '  DISP = ALL\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID.2,,1.,0.,0.\n'
            'GRID,3,,1.,1.,0.\n'
            'GRID,4,,0.,1.,0.\n'
            '$ helló wörld from two\n'
            'CQUAD4,1,2,3,4,5',
            #'ENDDATA'
        )


if __name__ == '__main__':   # pragma: no cover
    unittest.main()
