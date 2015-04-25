from __future__ import unicode_literals
import unittest
from six import PY2

import os
import pyNastran
from pyNastran.bdf.bdf import BDF

root_path = pyNastran.__path__[0]
test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

log = None
class TestReadWrite(unittest.TestCase):

    def test_write_1(self):
        """
        Tests 1 read method and various write methods
        """
        model = BDF(log=log, debug=False)

        bdf_name = os.path.join(test_path, 'test_mass.dat')
        model.read_bdf(bdf_name)
        model.write_bdf(os.path.join(test_path, 'test_mass1a.out'), size=8)
        model.write_bdf(os.path.join(test_path, 'test_mass2a.out'), size=8)

        model.write_bdf(os.path.join(test_path, 'test_mass1b.out'), size=8, interspersed=False)
        model.write_bdf(os.path.join(test_path, 'test_mass2b.out'), size=8, interspersed=True)
        os.remove(os.path.join(test_path, 'test_mass1a.out'))
        os.remove(os.path.join(test_path, 'test_mass2a.out'))
        os.remove(os.path.join(test_path, 'test_mass1b.out'))
        os.remove(os.path.join(test_path, 'test_mass2b.out'))

    def test_punch_1(self):
        """
        Tests punch file reading
        """
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        model.read_bdf(bdf_name, xref=False, punch=True)

        model2 = BDF(debug=False)
        #bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        model2.read_bdf(bdf_name, xref=False, punch=True)

    def test_read_include_dir_1(self):
        """
        Tests various read methods using various include files
        """
        # fails correctly
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'test_include.bdf')
        self.assertRaises(IOError, model.read_bdf, bdf_name, include_dir=None, xref=True, punch=False)

        # passes
        full_path = os.path.join(test_path, 'include_dir')
        model2 = BDF(debug=False)
        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(test_path, 'test_include.bdf')
        model2.read_bdf(bdf_filename, include_dir=full_path, xref=True, punch=False)

    def test_enddata_1(self):
        """
        There is an ENDDATA is in the baseline BDF, so None -> ENDDATA
        """
        model = BDF(debug=False)
        full_path = os.path.join(test_path, 'include_dir')
        model2 = BDF(debug=False)

        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(test_path, bdf_filename)
        model2.read_bdf(bdf_filename, include_dir=full_path, xref=True, punch=False)
        for out_filename, is_enddata, write_flag in [
            ('enddata1.bdf', True, None),
            ('enddata2.bdf', True, True),
            ('enddata3.bdf', False, False)]:
            out_filename = os.path.join(test_path, out_filename)
            model2.write_bdf(out_filename=out_filename+'.out', interspersed=True, size=8,
                             is_double=False, enddata=write_flag)
            data = open(out_filename + '.out', 'r').read()
            if is_enddata:
                self.assertTrue('ENDDATA' in data)
            else:
                self.assertFalse('ENDDATA' in data)
            os.remove(out_filename + '.out')

    def test_enddata_2(self):
        """
        There is no ENDDATA is in the baseline BDF, so None -> no ENDDATA
        """
        model = BDF(debug=False)
        full_path = os.path.join(test_path, 'include_dir')
        model2 = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'test_mass.dat')
        model2.read_bdf(bdf_name, include_dir=full_path, xref=True, punch=False)
        for out_filename, is_enddata, write_flag in [
            ('test_mass1.dat', False, None),
            ('test_mass2.dat', True, True),
            ('test_mass3.dat', False, False)]:
            model2.write_bdf(out_filename=out_filename, interspersed=True, size=8,
                             is_double=False, enddata=write_flag)
            data = open(out_filename, 'r').read()
            msg = 'outfilename=%r expected=%r write_flag=%s card_count=%r' % (out_filename, is_enddata, write_flag, model2.card_count.keys())
            if is_enddata:
                self.assertTrue('ENDDATA' in data, msg)
            else:
                self.assertFalse('ENDDATA' in data, msg)
            os.remove(out_filename)

    def test_include_end(self):
        """this test fails incorrectly"""
        if PY2:
            f = open('a.bdf', 'wb')
        else:
            f = open('a.bdf', 'w')
        f.write('CEND\n')
        f.write('BEGIN BULK\n')
        f.write('GRID,1\n')
        f.write("INCLUDE 'b.bdf'\n\n")

        if PY2:
            f = open('b.bdf', 'wb')
        else:
            f = open('b.bdf', 'w')
        f.write('GRID,2\n')
        f.write("INCLUDE 'c.bdf'\n\n")

        if PY2:
            f = open('c.bdf', 'wb')
        else:
            f = open('c.bdf', 'w')
        f.write('GRID,3\n\n')
        f.write("ENDDATA\n")
        f.close()

        model = BDF(log=log, debug=False)
        model.read_bdf('a.bdf')
        model.write_bdf('a.out.bdf')

        os.remove('a.bdf')
        os.remove('b.bdf')
        os.remove('c.bdf')
        os.remove('a.out.bdf')
        self.assertEqual(len(model.nodes), 3)

    def test_read_bad_01(self):
        model = BDF()
        model.active_filenames = ['fake.file']
        with self.assertRaises(IOError):
            model._open_file('fake.file')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()