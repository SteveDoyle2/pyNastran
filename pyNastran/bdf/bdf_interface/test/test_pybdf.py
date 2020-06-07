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
    BDFInputPy, _show_bad_file, _lines_to_decks, MissingDeckSections)


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
            pybdf._validate_open_file(bdf_filename, op2_filename, check=True)

        with self.assertRaises(RuntimeError):
            _show_bad_file(pybdf, op2_filename, encoding, nlines_previous=10)

        #------------------------------------------------------
        with self.assertRaises(IOError):
            pybdf._validate_open_file(bdf_filename, op2_filename, check=True)
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
            pybdf._validate_open_file(bdf_filename, bdf_dir, check=True)

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

    def test_no_punch(self):
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
        bulk_data_lines = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)[3]
        assert len(bulk_data_lines) == 5, bulk_data_lines
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
        bulk_data_lines = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)[3]
        assert len(bulk_data_lines) == 5, bulk_data_lines

    def test_unicode_errors1(self):
        """tests some error handling"""
        bdf_filename = 'unicode.bdf'
        with open(bdf_filename, 'w') as bdf_file:
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

        read_includes = True
        dumplines = True
        encoding = 'utf8'
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='msc',
                           consider_superelements=False, log=None, debug=False)
        if sys.platform == 'win32':
            with self.assertRaises(RuntimeError):
                pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)
        else:
            pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)

        #with self.assertRaises(RuntimeError):
        #with self.assertRaises(UnicodeDecodeError):
        #unused_bulk_data_lines1 = pybdf.get_lines(bdf_filename, punch=None, make_ilines=True)[3]

    def test_unicode_errors2(self):
        """tests some error handling"""
        bdf_filename = 'unicode.bdf'
        with open(bdf_filename, 'w') as bdf_file:
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

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
