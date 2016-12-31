from __future__ import unicode_literals, print_function
import os
from codecs import open as codec_open
import unittest
from six import PY2, StringIO

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, get_logger2
from pyNastran.bdf.test.test_case_control_deck import compare_lines

root_path = pyNastran.__path__[0]
test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

log = get_logger2(debug=None)

class TestReadWrite(unittest.TestCase):
    """various BDF I/O tests"""

    def test_write_1(self):
        """
        Tests 1 read method and various write methods
        """
        model = BDF(log=log, debug=False)

        bdf_name = os.path.join(test_path, 'test_mass.dat')
        model.read_bdf(bdf_name)
        model.write_bdf(os.path.join(test_path, 'test_mass1a.out'), size=8)
        model.write_bdf(os.path.join(test_path, 'test_mass2a.out'), size=8)
        msg = model.get_bdf_stats(return_type='list')
        str('\n'.join(msg))

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
        model = BDF(log=log, debug=False)
        bdf_name = os.path.join(test_path, 'include_dir', 'include_alt.inc')
        model.read_bdf(bdf_name, xref=False, punch=True)

        model2 = BDF(log=log, debug=False)
        #bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        model2.read_bdf(bdf_name, xref=False, punch=True)

    def test_read_include_dir_1(self):
        """
        Tests various read methods using various include files
        """
        # fails correctly
        model = BDF(log=log, debug=False)
        bdf_name = os.path.join(test_path, 'test_include.bdf')
        model.read_bdf(bdf_name, xref=True, punch=False)
        #self.assertRaises(IOError, model.read_bdf, bdf_name, xref=True, punch=False)

        # passes
        #full_path = os.path.join(test_path, 'include_dir')
        model2 = BDF(log=log, debug=False)
        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(test_path, 'test_include.bdf')
        model2.read_bdf(bdf_filename, xref=True, punch=False)

    def test_read_include_dir_2(self):
        full_path = os.path.join(test_path)
        model = BDF(log=log, debug=False)
        bdf_filename = 'test_include2.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(full_path, 'test_include2.bdf')
            #print(full_path)
        #print(bdf_filename)
        model.read_bdf(bdf_filename, xref=True, punch=False)
        #model.write_bdf('junk.bdf')


    def test_enddata_1(self):
        """
        There is an ENDDATA is in the baseline BDF, so None -> ENDDATA
        """
        model2 = BDF(log=log, debug=False)

        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(test_path, bdf_filename)
        model2.read_bdf(bdf_filename, xref=True, punch=False)

        cases = [
            ('enddata1.bdf', True, None),
            ('enddata2.bdf', True, True),
            ('enddata3.bdf', False, False),
        ]
        for out_filename, is_enddata, write_flag in cases:
            out_filename = os.path.join(test_path, out_filename)
            model2.write_bdf(out_filename=out_filename+'.out', interspersed=True, size=8,
                             is_double=False, enddata=write_flag)

            with codec_open(out_filename + '.out', 'r') as bdf_file:
                data = bdf_file.read()

            if is_enddata:
                self.assertTrue('ENDDATA' in data)
            else:
                self.assertFalse('ENDDATA' in data)
            os.remove(out_filename + '.out')

    def test_enddata_2(self):
        """
        There is no ENDDATA is in the baseline BDF, so None -> no ENDDATA
        """
        model2 = BDF(log=log, debug=False)
        bdf_name = os.path.join(test_path, 'test_mass.dat')
        model2.read_bdf(bdf_name, xref=True, punch=False)

        cases = [
            ('test_mass1.dat', False, None),
            ('test_mass2.dat', True, True),
            ('test_mass3.dat', False, False)
        ]
        for out_filename, is_enddata, write_flag in cases:
            model2.write_bdf(out_filename=out_filename, interspersed=True, size=8,
                             is_double=False, enddata=write_flag)

            with codec_open(out_filename, 'r') as bdf_file:
                data = bdf_file.read()

            msg = 'outfilename=%r expected=%r write_flag=%s card_count=%r' % (
                out_filename, is_enddata, write_flag, model2.card_count.keys())
            if is_enddata:
                self.assertTrue('ENDDATA' in data, msg)
            else:
                self.assertFalse('ENDDATA' in data, msg)
            os.remove(out_filename)

    def test_add_card_skip(self):
        model = BDF(log=log, debug=False)

        card_name = 'JUNK'
        card_lines1 = ['JUNK', 1, 2, 3]
        card_lines2 = ['JUNK,a,b,c']
        model.add_card(card_lines1, card_name)
        model.add_card(card_lines2, card_name, comment='', is_list=False,
                       has_none=True)
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        msg = bdf_file.getvalue()
        bdf_file.close()
        assert 'JUNK           1       2       3' in msg, msg
        assert 'JUNK           a       b       c' in msg, msg

    def test_add_card_fail(self):
        model = BDF(log=log, debug=False)
        card_lines1 = ['GRID', 1, 'a', 'b', 'c']
        card_lines2 = ['GRID', 1, 'd', 'e', 'f']
        with self.assertRaises(SyntaxError):
            model.add_card(card_lines1, 'GRID')

        model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                                nxref_errors=100,
                                stop_on_xref_error=True)

        card_lines3 = ['GMSPC', 1, 'd', 'e', 'f']
        model.add_card(card_lines3, 'GMSPC')

        card_lines4 = ['GRDSET', 1, 'd2', 'e2', 'f2']
        model.add_card(card_lines4, 'GRDSET')

    def test_include_end(self):
        with codec_open('a.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")

        with codec_open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")

        with codec_open('c.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,3,,3.0\n\n')
            bdf_file.write("ENDDATA\n")

        model = BDF(log=log, debug=False)
        model.read_bdf('a.bdf')
        model.write_bdf('a.out.bdf')

        os.remove('a.bdf')
        os.remove('b.bdf')
        os.remove('c.bdf')
        os.remove('a.out.bdf')
        self.assertEqual(len(model.nodes), 3)
        self.assertEqual(model.nnodes, 3, 'nnodes=%s' % model.nnodes)

    def test_include_end_02(self):
        with codec_open('a.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")
            bdf_file.write('GRID,4,,4.0\n')

        with codec_open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")
            bdf_file.write('GRID,5,,5.0\n')

        with codec_open('c.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,3,,3.0\n\n')

        model = BDF(log=log, debug=False)
        model.read_bdf('a.bdf')
        model.write_bdf('a.out.bdf')

        os.remove('a.bdf')
        os.remove('b.bdf')
        os.remove('c.bdf')
        os.remove('a.out.bdf')
        self.assertEqual(len(model.nodes), 5)
        self.assertEqual(model.nnodes, 5, 'nnodes=%s' % model.nnodes)

    def test_include_03(self):
        with codec_open('a.bdf', 'w') as bdf_file:
            bdf_file.write("INCLUDE 'executive_control.inc'\n\n")
            bdf_file.write('CEND\n')
            bdf_file.write("INCLUDE 'case_control.inc'\n\n")
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")
            bdf_file.write('GRID,4,,4.0\n')

        with codec_open('executive_control.inc', 'w') as bdf_file:
            bdf_file.write('SOL = 103\n')

        with codec_open('case_control.inc', 'w') as bdf_file:
            bdf_file.write('DISP = ALL\n')

        with codec_open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")
            bdf_file.write('GRID,5,,5.0\n')

        with codec_open('c.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,3,,3.0\n\n')

        model = BDF(log=log, debug=False)
        model.read_bdf('a.bdf')
        model.write_bdf('a.out.bdf')

        os.remove('a.bdf')
        os.remove('b.bdf')
        os.remove('c.bdf')
        os.remove('executive_control.inc')
        os.remove('case_control.inc')

        os.remove('a.out.bdf')
        self.assertEqual(len(model.nodes), 5)
        self.assertEqual(model.nnodes, 5, 'nnodes=%s' % model.nnodes)

    def test_include_04(self):
        with codec_open('include4.bdf', 'w') as bdf_file:
            bdf_file.write('$ pyNastran: punch=True\n')
            bdf_file.write('$ pyNastran: dumplines=True\n')
            bdf_file.write("INCLUDE 'include4b.inc'\n\n")

        with codec_open('include4b.inc', 'w') as bdf_file:
            bdf_file.write('$ GRID comment\n')
            bdf_file.write('GRID,2,,2.0\n')

        model = BDF(log=log, debug=False)
        model.read_bdf('include4.bdf')
        model.write_bdf('include4.out.bdf')

        os.remove('include4.out.bdf')
        os.remove('include4b.inc')
        os.remove('include4.bdf')
        #os.remove('include4.inc')
        # os.remove('c.bdf')
        # os.remove('executive_control.inc')
        # os.remove('case_control.inc')

        self.assertEqual(len(model.nodes), 1)
        self.assertEqual(model.nnodes, 1, 'nnodes=%s' % model.nnodes)

    def test_include_05(self):
        with codec_open('include5.bdf', 'w') as bdf_file:
            bdf_file.write('$ pyNastran: punch=True\n')
            bdf_file.write('$ pyNastran: dumplines=True\n')
            bdf_file.write("INCLUDE 'include5b.inc'\n\n")

        with codec_open('include5b.inc', 'w') as bdf_file:
            bdf_file.write('ECHOON\n')
            bdf_file.write('$ GRID comment\n')
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write('ECHOOFF\n')
            bdf_file.write('GRID,3,,3.0\n')
            bdf_file.write('grid,4,,4.0\n')
            bdf_file.write('grid ,5,,5.0\n')

        model = BDF(log=log, debug=False)
        model.read_bdf('include5.bdf')
        assert model.echo is False, model.echo
        #model.write_bdf('include5.out.bdf')

        # os.remove('c.bdf')
        # os.remove('executive_control.inc')
        # os.remove('case_control.inc')

        self.assertEqual(len(model.nodes), 4)
        self.assertEqual(model.nnodes, 4, 'nnodes=%s' % model.nnodes)

        model2 = read_bdf(bdf_filename='include5.bdf', xref=True, punch=False,
                          log=log, encoding=None)
        self.assertEqual(len(model2.nodes), 4)
        self.assertEqual(model2.nnodes, 4, 'nnodes=%s' % model.nnodes)
        os.remove('include5.bdf')
        #os.remove('include5.out.bdf')
        os.remove('include5b.inc')


    def test_encoding_write(self):
        mesh = BDF(log=log, debug=False)
        mesh.add_card(['GRID', 100000, 0, 43.91715, -29., .8712984], 'GRID')
        mesh.write_bdf('out.bdf')
        lines_expected = [
            '$pyNastran: version=msc',
            '$pyNastran: punch=True',
            '$pyNastran: encoding=ascii' if PY2 else '$pyNastran: encoding=utf-8\n',
            '$pyNastran: nnodes=1',
            '$pyNastran: nelements=0',
            '$NODES',
            'GRID      100000        43.91715    -29..8712984',
        ]
        bdf_filename = 'out.bdf'
        with codec_open(bdf_filename, 'r', encoding='ascii') as bdf_file:
            lines = bdf_file.readlines()
            compare_lines(self, lines, lines_expected, has_endline=False)


    def test_include_stop(self):
        with codec_open('a.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")
            bdf_file.write('GRID,1,,1.0\n')
        model = BDF(log=log, debug=False)
        with self.assertRaises(IOError):
            model.read_bdf(bdf_filename='a.bdf', xref=True, punch=False,
                           read_includes=True, encoding=None)
        with self.assertRaises(IOError):
            read_bdf(bdf_filename='a.bdf', xref=True, punch=False,
                     encoding=None, log=log)
        model.read_bdf(bdf_filename='a.bdf', xref=True, punch=False,
                       read_includes=False, encoding=None)
        model.write_bdf('out.bdf')
        os.remove('a.bdf')
        os.remove('out.bdf')

    def test_read_bad_01(self):
        """tests you can't read the same file twice"""
        model = BDF(log=log, debug=False)
        model.active_filenames = ['fake.file']
        with self.assertRaises(IOError):
            model._open_file('fake.file')

    def test_read_bad_02(self):
        """tests when users don't add punch=True to read_bdf(...)"""
        lines = [
            'GRID     1000177       0      1.      0.      0.       0\n',
            'GRID     1000178       0      0.      1.      0.       0\n',
            'GRID     1000186       0      0.      0.      1.       0\n',
            'GRID     1000187       0      1.      1.      1.       0\n',
            'GRID    15000014       0      2.      1.      1.       0\n',
            'RBE2    1500002215000014  123456 1000177 1000178 1000186 1000187\n',
        ]
        bdf_filename = 'xref_test.bdf'
        with codec_open(bdf_filename, 'w') as bdf_file:
            bdf_file.writelines(lines)
        with self.assertRaises(RuntimeError):
            read_bdf(bdf_filename, validate=False, xref=False,
                     punch=False, encoding=None,
                     log=log, debug=True, mode='msc')
        os.remove(bdf_filename)

    def test_disable_cards(self):
        """tests disabling cards"""
        bdf_filename = os.path.join(root_path, '..', 'models',
                                    'solid_bending', 'solid_bending.bdf')
        model = BDF(log=log, debug=False)
        model.disable_cards(['CTETRA'])
        model.read_bdf(bdf_filename)
        assert len(model.elements) == 0, len(model.elements)

    def test_solid_shell_bar_buckling(self):
        bdf_filename = os.path.join(root_path, '..', 'models',
                                    'sol_101_elements', 'buckling_solid_shell_bar.bdf')
        bdf_filename2 = os.path.join(root_path, '..', 'models',
                                     'sol_101_elements', 'buckling_solid_shell_bar2.bdf')
        model = BDF(debug=False)
        model.read_bdf(bdf_filename)
        model.write_bdf(bdf_filename2)

        model2 = BDF(debug=False)
        model2.read_bdf(bdf_filename2)
        #print(model2.get_bdf_stats())

        eigb = model2.methods[42]
        assert eigb.comment == '$ this is a preload buckling case\n', 'comment=%r\n%s' % (eigb.comment, str(eigb))
        os.remove(bdf_filename2)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
