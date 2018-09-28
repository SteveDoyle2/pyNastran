# encoding: utf8
# pylint: disable=W0212
from __future__ import print_function, unicode_literals
import os
from io import open
import unittest
from six import StringIO

from pyNastran.bdf.bdf_interface.pybdf import BDFInputPy, _show_bad_file, _lines_to_decks

class TestPyBDF(unittest.TestCase):

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

    def test_get_lines(self):
        with open('junk.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1\n')

        bdf_filename = StringIO()
        bdf_filename.write(
            'ASSIGN FEM="junk.bdf"\n'
            'ASSIGN MATRIX=junk.mat\n'
            'ASSIGN OUTPUT4=junk.op4\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           log=None, debug=False)
        pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)
        #----------------------------------------------------------------
        bdf_filename = StringIO()
        bdf_filename.write(
            'ASSIGN FEM=junk.f06\n'
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           log=None, debug=False)
        pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)

        #----------------------------------------------------------------
        bdf_filename = StringIO()
        bdf_filename.write(
            "ASSIGN FEM='junk.f06'\n"
            'CEND\n'
            'BEGIN BULK\n'
        )
        bdf_filename.seek(0)

        read_includes = True
        dumplines = True
        encoding = None
        pybdf = BDFInputPy(read_includes, dumplines, encoding, nastran_format='zona',
                           log=None, debug=False)
        pybdf.get_lines(bdf_filename, punch=False, make_ilines=True)

if __name__ == '__main__':
    unittest.main()
