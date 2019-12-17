import os
import unittest
from io import StringIO

from cpylog import get_logger
import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.bdf_interface.pybdf import BDFInputPy
from pyNastran.bdf.bdf_interface.include_file import (
    split_filename_into_tokens, get_include_filename,
    PurePosixPath, PureWindowsPath,
) # ,_split_to_tokens
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf_interface.test.test_case_control_deck import compare_lines

ROOT_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(ROOT_PATH, 'bdf', 'test', 'unit')
MESH_UTILS_PATH = os.path.join(ROOT_PATH, 'bdf', 'mesh_utils', 'test')
MODEL_PATH = os.path.join(ROOT_PATH, '../', 'models')


class TestReadWriteFiles(unittest.TestCase):
    def test_read_include_dir_2_files(self):
        """tests a model that will write to multiple files"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        full_path = os.path.join(TEST_PATH)
        model = BDF(log=log, debug=False)
        bdf_filename = 'test_include2.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(full_path, 'test_include2.bdf')
            #print(full_path)
        #print(bdf_filename)
        model.read_bdf(bdf_filename, xref=True, punch=False, save_file_structure=True)

        all_lines, ilines = model.include_zip(bdf_filename, encoding=None, make_ilines=True)
        #for (ifile, iline), line in zip(ilines, all_lines):
            #print(ifile, iline, line.rstrip())
        #print(ilines)
        #dd

        out_filenames = {}
        for filename in model.active_filenames:
            dirname = os.path.dirname(filename)
            basename = os.path.basename(filename)
            filename2 = os.path.join(dirname, 'out_' + basename)
            abs_path = os.path.abspath(filename)
            out_filenames[abs_path] = filename2

        model.write_bdfs(out_filenames, encoding=None,
                         size=8, is_double=False,
                         enddata=None, close=True, relative_dirname='')
        #read_bdf('out_test_include2.bdf')

        model.write_bdfs(out_filenames, encoding=None,
                         size=8, is_double=False,
                         enddata=None, close=True, relative_dirname=None)
        #read_bdf('out_test_include2.bdf')

        out_filenames2 = {abs_path : filename2}
        model.write_bdfs(out_filenames2, encoding=None,
                         size=8, is_double=False,
                         enddata=None, close=True)
        #read_bdf('out_test_include2.bdf')
        #os.remove('out_test_include2.bdf')


    def test_isat_files(self):
        """read/writes the isat model with the file structure"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'iSat', 'iSat_launch_100Hz.dat')
        model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                         save_file_structure=True, skip_cards=None, read_cards=None,
                         encoding=None, log=log, debug=True, mode='msc')
        assert len(model.include_filenames) == 1, len(model.include_filenames)
        assert len(model.include_filenames[0]) == 2, len(model.include_filenames[0])

        out_filenames = {}
        out_filenames2 = {}
        for i, fname in enumerate(model.active_filenames):
            dirname = os.path.dirname(fname)
            basename = os.path.basename(fname)
            out_filenames[fname] = os.path.join(dirname, 'out_' + basename)
            if 'antenna_pressure' not in fname:
                out_filenames2[fname] = os.path.join(dirname, 'out2_' + basename)
            #print(bdf_filename2)
        #print('out_filenames =', out_filenames)

        all_lines, ilines = model.include_zip(bdf_filename, encoding=None, make_ilines=True)
        #for (ifile, iline), line in zip(ilines, all_lines):
            #if iline > 100:
                #continue
            #print(ifile, iline, line.rstrip())

        assert len(model.include_filenames) == 1, len(model.include_filenames)
        assert len(model.include_filenames[0]) == 2, len(model.include_filenames[0])
        model.log.debug('saving model')
        model.write_bdfs(out_filenames, relative_dirname='')

        model.log.debug('saving new model')
        model.write_bdfs(out_filenames2, relative_dirname='')

        bdf_filename = os.path.abspath(bdf_filename)
        read_bdf(out_filenames[bdf_filename])
        read_bdf(out_filenames2[bdf_filename])

class TestReadWrite(unittest.TestCase):
    """various BDF I/O tests"""

    def test_write_1(self):
        """
        Tests 1 read method and various write methods
        """
        log = get_logger(log=None, level='info', encoding='utf-8')
        model = BDF(log=log, debug=False)

        bdf_name = os.path.join(MESH_UTILS_PATH, 'test_mass.dat')
        model.read_bdf(bdf_name)
        model.write_bdf(os.path.join(MESH_UTILS_PATH, 'test_mass1a.out'), size=8)
        model.write_bdf(os.path.join(MESH_UTILS_PATH, 'test_mass2a.out'), size=8)
        msg = model.get_bdf_stats(return_type='list')
        str('\n'.join(msg))

        model.write_bdf(os.path.join(MESH_UTILS_PATH, 'test_mass1b.out'), size=8, interspersed=False)
        model.write_bdf(os.path.join(MESH_UTILS_PATH, 'test_mass2b.out'), size=8, interspersed=True)
        os.remove(os.path.join(MESH_UTILS_PATH, 'test_mass1a.out'))
        os.remove(os.path.join(MESH_UTILS_PATH, 'test_mass2a.out'))
        os.remove(os.path.join(MESH_UTILS_PATH, 'test_mass1b.out'))
        os.remove(os.path.join(MESH_UTILS_PATH, 'test_mass2b.out'))

    def test_punch_1(self):
        """Tests punch file reading"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        model = BDF(log=log, debug=False)
        bdf_name = os.path.join(TEST_PATH, 'include_dir', 'include_alt.inc')
        model.read_bdf(bdf_name, xref=False, punch=True)

        model2 = BDF(log=log, debug=False)
        #bdf_name = os.path.join(TEST_PATH, 'include_dir', 'include.inc')
        model2.read_bdf(bdf_name, xref=False, punch=True)

    def test_read_include_dir_1(self):
        """Tests various read methods using various include files"""
        # fails correctly
        log = get_logger(log=None, level='info', encoding='utf-8')
        model = BDF(log=log, debug=False)
        bdf_name = os.path.join(TEST_PATH, 'test_include.bdf')
        model.read_bdf(bdf_name, xref=True, punch=False)
        #self.assertRaises(IOError, model.read_bdf, bdf_name, xref=True, punch=False)

        # passes
        #full_path = os.path.join(TEST_PATH, 'include_dir')
        model2 = BDF(log=log, debug=False)
        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(TEST_PATH, 'test_include.bdf')
        model2.read_bdf(bdf_filename, xref=True, punch=False)

    def test_read_include_dir_2(self):
        full_path = os.path.join(TEST_PATH)
        log = get_logger(log=None, level='info', encoding='utf-8')
        model = BDF(log=log, debug=False)
        bdf_filename = 'test_include2.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(full_path, 'test_include2.bdf')
            #print(full_path)
        #print(bdf_filename)
        model.read_bdf(bdf_filename, xref=True, punch=False)
        #model.write_bdf('junk.bdf')

    def test_enddata_1(self):
        """There is an ENDDATA is in the baseline BDF, so None -> ENDDATA"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        model2 = BDF(log=log, debug=False)

        bdf_filename = 'test_include.bdf'
        if not os.path.exists(bdf_filename):
            bdf_filename = os.path.join(TEST_PATH, bdf_filename)
        model2.read_bdf(bdf_filename, xref=True, punch=False)

        cases = [
            ('enddata1.bdf', True, None),
            ('enddata2.bdf', True, True),
            ('enddata3.bdf', False, False),
        ]
        for out_filename, is_enddata, write_flag in cases:
            out_filename = os.path.join(TEST_PATH, out_filename)
            bdf_filename_out = out_filename + '.out'
            model2.write_bdf(out_filename=bdf_filename_out, interspersed=True, size=8,
                             is_double=False, enddata=write_flag)

            with open(out_filename + '.out', 'r') as bdf_file:
                data = bdf_file.read()

            if is_enddata:
                self.assertTrue('ENDDATA' in data)
            else:
                self.assertFalse('ENDDATA' in data)
            os.remove(bdf_filename_out)

    def test_enddata_2(self):
        """
        There is no ENDDATA is in the baseline BDF, so None -> no ENDDATA
        """
        log = get_logger(log=None, level='info', encoding='utf-8')
        model2 = BDF(log=log, debug=False)
        bdf_name = os.path.join(MESH_UTILS_PATH, 'test_mass.dat')
        model2.read_bdf(bdf_name, xref=True, punch=False)

        cases = [
            ('test_mass1.dat', False, None),
            ('test_mass2.dat', True, True),
            ('test_mass3.dat', False, False)
        ]
        for out_filename, is_enddata, write_flag in cases:
            model2.write_bdf(out_filename=out_filename, interspersed=True, size=8,
                             is_double=False, enddata=write_flag)

            with open(out_filename, 'r') as bdf_file:
                data = bdf_file.read()

            msg = 'outfilename=%r expected=%r write_flag=%s card_count=%r' % (
                out_filename, is_enddata, write_flag, model2.card_count.keys())
            if is_enddata:
                self.assertTrue('ENDDATA' in data, msg)
            else:
                self.assertFalse('ENDDATA' in data, msg)
            os.remove(out_filename)

    def test_add_card_skip(self):
        """tests that a fake card 'JUNK' is skipped"""
        log = get_logger(log=None, level='info', encoding='utf-8')
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
        log = get_logger(log=None, level='info', encoding='utf-8')
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
        """tests multiple levels of includes"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('a.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")

        with open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")

        with open('c.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,3,,3.0\n\n')
            bdf_file.write("ENDDATA\n")

        model = BDF(log=log, debug=False)
        model.read_bdf('a.bdf')
        model.write_bdf('a.out.bdf')

        self.assertEqual(len(model.nodes), 3)
        self.assertEqual(model.nnodes, 3, 'nnodes=%s' % model.nnodes)

        model = BDF(log=log, debug=False)
        lines, ilines = model.include_zip(bdf_filename='a.bdf', encoding=None)
        assert len(lines) == 11, len(lines)

        os.remove('a.bdf')
        os.remove('b.bdf')
        os.remove('c.bdf')
        os.remove('a.out.bdf')

    def test_include_end_02(self):
        """tests multiple levels of includes"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('a.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")
            bdf_file.write('GRID,4,,4.0\n')

        with open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")
            bdf_file.write('GRID,5,,5.0\n')

        with open('c.bdf', 'w') as bdf_file:
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
        """tests executive/case control includes"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('a.bdf', 'w') as bdf_file:
            bdf_file.write("INCLUDE 'executive_control.inc'\n\n")
            bdf_file.write('CEND\n')
            bdf_file.write("INCLUDE 'case_control.inc'\n\n")
            bdf_file.write('BEGIN BULK\n')
            bdf_file.write('GRID,1,,1.0\n')
            bdf_file.write("INCLUDE 'b.bdf'\n\n")
            bdf_file.write('GRID,4,,4.0\n')

        with open('executive_control.inc', 'w') as bdf_file:
            bdf_file.write('SOL = 103\n')

        with open('case_control.inc', 'w') as bdf_file:
            bdf_file.write('DISP = ALL\n')

        with open('b.bdf', 'w') as bdf_file:
            bdf_file.write('GRID,2,,2.0\n')
            bdf_file.write("INCLUDE 'c.bdf'\n\n")
            bdf_file.write('GRID,5,,5.0\n')

        with open('c.bdf', 'w') as bdf_file:
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
        """tests pyNastran: punch=True with includes"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('include4.bdf', 'w') as bdf_file:
            bdf_file.write('$ pyNastran: punch=True\n')
            bdf_file.write('$ pyNastran: dumplines=True\n')
            bdf_file.write("INCLUDE 'include4b.inc'\n\n")

        with open('include4b.inc', 'w') as bdf_file:
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
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('include5.bdf', 'w') as bdf_file:
            bdf_file.write('$ pyNastran: punch=True\n')
            bdf_file.write('$ pyNastran: dumplines=True\n')
            bdf_file.write("INCLUDE 'include5b.inc'\n\n")

        with open('include5b.inc', 'w') as bdf_file:
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
        os.remove('pyNastran_dump.bdf')


    def test_encoding_write(self):
        """tests encodings in BDF header"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        mesh = BDF(log=log, debug=False)
        mesh.add_card(['GRID', 100000, 0, 43.91715, -29., .8712984], 'GRID')
        mesh.write_bdf('out.bdf')
        lines_expected = [
            '$pyNastran: version=msc',
            '$pyNastran: punch=True',
            '$pyNastran: encoding=utf-8\n',
            '$pyNastran: nnodes=1',
            '$pyNastran: nelements=0',
            '$NODES',
            'GRID      100000        43.91715    -29..8712984',
        ]
        bdf_filename = 'out.bdf'
        with open(bdf_filename, 'r', encoding='ascii') as bdf_file:
            lines = bdf_file.readlines()
            compare_lines(self, lines, lines_expected, has_endline=False)


    def test_include_stop(self):
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open('a.bdf', 'w') as bdf_file:
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
        os.remove('pyNastran_crash.bdf')

    def test_read_bad_01(self):
        """tests you can't read the same file twice"""
        read_includes = True
        dumplines = False
        encoding = 'ascii'
        log = get_logger(log=None, level='info', encoding='utf-8')
        model = BDFInputPy(read_includes, dumplines, encoding, log=log, debug=False)
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
        log = get_logger(log=None, level='info', encoding='utf-8')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.writelines(lines)
        with self.assertRaises(RuntimeError):
            read_bdf(bdf_filename, validate=False, xref=False,
                     punch=False, encoding=None,
                     log=log, debug=True, mode='msc')
        os.remove(bdf_filename)

    def test_disable_cards(self):
        """tests disabling cards"""
        log = get_logger(log=None, level='info', encoding='utf-8')
        bdf_filename = os.path.join(ROOT_PATH, '..', 'models',
                                    'solid_bending', 'solid_bending.bdf')
        model = BDF(log=log, debug=False)
        model.disable_cards(['CTETRA'])
        model.read_bdf(bdf_filename)
        assert len(model.elements) == 0, len(model.elements)

    def test_solid_shell_bar_buckling(self):
        bdf_filename = os.path.join(ROOT_PATH, '..', 'models',
                                    'sol_101_elements', 'buckling_solid_shell_bar.bdf')
        bdf_filename2 = os.path.join(ROOT_PATH, '..', 'models',
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

    def test_paths(self):
        """tests parsing paths"""
        include_dir = ''
        filename = 'model.bdf'
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        assert filename_out == PureWindowsPath('model.bdf'), 'filename_out=%r windows' % filename_out

        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        assert filename_out == PurePosixPath('model.bdf'), 'filename_out=%r linux/mac' % filename_out
        #-----------------------------------------------------------------------

        include_dir = 'dir'
        filename = 'model.bdf'
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        assert filename_out == PureWindowsPath('dir/model.bdf'), 'filename_out=%r windows' % filename_out

        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        assert filename_out == PurePosixPath('dir/model.bdf'), 'filename_out=%r linux/mac' % filename_out
        #-----------------------------------------------------------------------
        include_dir = 'dir1'
        filename = 'dir2/model.bdf'
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        assert filename_out == PureWindowsPath('dir1', 'dir2', 'model.bdf'), 'filename_out=%r windows' % filename_out

        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        assert filename_out == PurePosixPath('dir1', 'dir2', 'model.bdf'), 'filename_out=%r linux/mac' % filename_out
        #-----------------------------------------------------------------------
        include_dir = 'dir1/'
        filename = '/dir2/model.bdf'
        with self.assertRaises(SyntaxError):
            filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        #assert filename_out == PureWindowsPath('dir1', 'dir2', 'model.bdf'), 'filename_out=%r windows' % filename_out

        #with self.assertRaises(SyntaxError):
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        assert filename_out == PurePosixPath('/dir2', 'model.bdf'), 'filename_out=%r linux/mac' % filename_out

        #-----------------------------------------------------------------------
        include_dir = ''
        filename = 'C:/dir2/model.bdf'
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        assert filename_out == PureWindowsPath('C:/dir2/model.bdf'), 'filename_out=%r linux/mac' % filename_out

        with self.assertRaises(SyntaxError):
            filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        #assert filename_out == PurePosixPath('C:', 'dir2', 'model.bdf'), 'filename_out=%r linux/mac' % filename_out

        #-----------------------------------------------------------------------
        include_dir = ''
        filename = '/dir2/model.bdf'
        with self.assertRaises(SyntaxError):
            filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        #assert filename_out == PureWindowsPath('dir2', 'model.bdf'), 'filename_out=%r linux/mac' % filename_out

        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        assert filename_out == PurePosixPath('/dir2', 'model.bdf'), 'filename_out=%r linux/mac' % filename_out
        #-----------------------------------------------------------------------
        include_dir = ''
        filename = '\\\\nas\\dir2\\model.bdf'
        filename_out = split_filename_into_tokens(include_dir, filename, is_windows=True)
        assert filename_out == PureWindowsPath('\\\\nas\\dir2\\model.bdf'), 'filename_out=%r linux/mac' % filename_out

        with self.assertRaises(SyntaxError):
            filename_out = split_filename_into_tokens(include_dir, filename, is_windows=False)
        #assert filename_out == PurePosixPath('\\\\nas\\dir2\\model.bdf'), 'filename_out=%r linux/mac' % filename_out
        #-----------------------------------------------------------------------


    def test_paths_sat(self):
        """runs through the various satellite includes on windows and linux"""
        sat_path = os.path.abspath(os.path.join(MODEL_PATH, 'Satellite_V02'))
        os.environ['Satellite_V02_base'] = sat_path
        os.environ['Satellite_V02_bddm'] = os.path.join(sat_path, 'BULK', 'MATERIAUX')
        os.environ['Satellite_V02_BULK'] = os.path.join(sat_path, 'BULK')
        os.environ['Satellite_V02_INCLUDE'] = os.path.join(sat_path, 'INCLUDE')

        assert os.path.exists(sat_path), print_bad_path(sat_path)
        assert os.path.exists(os.path.join(sat_path, 'BULK')
                              ), print_bad_path(os.path.join(sat_path, 'BULK'))
        assert os.path.exists(os.path.join(sat_path, 'BULK', 'MATERIAUX')
                              ), print_bad_path(os.path.join(sat_path, 'BULK', 'MATERIAUX'))
        assert os.path.exists(os.path.join(sat_path, 'INCLUDE')
                              ), print_bad_path(os.path.join(sat_path, 'INCLUDE'))

        pths = [
            "INCLUDE 'Satellite_V02_bddm:Satellite_V02_Materiaux.blk'",
            "INCLUDE 'Satellite_V02_BULK:CONM2/Satellite_V02_CONM2.blk'",
            "INCLUDE 'Satellite_V02_BULK:RBE2/Satellite_V02_RBE2.blk'",
            "INCLUDE 'Satellite_V02_BULK:TUBE/Satellite_V02_TubeCentral.blk'",
            "INCLUDE 'Satellite_V02_BULK:TUBE/Satellite_V02_Barre_TubeCentral.blk'",
            "INCLUDE 'Satellite_V02_INCLUDE:Satellite_V02_Panneau_Etoile.dat'",
            "INCLUDE 'Satellite_V02_BULK:ETOILE/Satellite_V02_Barre_Panneau_Etoile.blk'",
            "INCLUDE 'Satellite_V02_BULK:TOP/Satellite_V02_Panneau_PZ.blk'",
            "INCLUDE 'Satellite_V02_BULK:TOP/Satellite_V02_Barre_Panneau_PZ.blk'",
            "INCLUDE 'Satellite_V02_BULK:BOTTOM/Satellite_V02_Panneau_MZ.blk'",
            "INCLUDE 'Satellite_V02_INCLUDE:Satellite_V02_Tube_Cone.dat'",
            "INCLUDE 'Satellite_V02_BULK:COORDS/satellite_V02_Coord.blk'",
            "INCLUDE 'Satellite_V02_INCLUDE:Satellite_V02_Panneau_Externe.dat'",
        ]
        for pth in pths:
            pth2 = get_include_filename([pth], include_dir='', is_windows=True)
            #if not os.path.exists(pth2):
                #msg = 'Invalid Path\nold:  %r\nnew:  %r' % (pth, pth2)
                #msg += print_bad_path(pth2)
                #raise RuntimeError(msg)
            #print('pth1 =', pth2)

            pth2 = get_include_filename([pth], include_dir='', is_windows=False)
            #print('pth2 =', pth2, '\n')
        #filename_tokens = _split_to_tokens(r'\\nas3\dir1\dir2', is_windows=True)

        #Satellite_V02_base = M:\ACA\Satellite_V02
        #Satellite_V02_bddm = M:\ACA\Satellite_V02/BULK/MATERIAUX
        #Satellite_V02_BULK = M:\ACA\Satellite_V02/INCLUDE
        #Satellite_V02_INCLUDE = M:\ACA\Satellite_V02/BULK
        #split_filename_into_tokens

    def test_paths_sat_02(self):
        """the satellite model should work"""
        sat_path = os.path.abspath(os.path.join(MODEL_PATH, 'Satellite_V02'))
        os.environ['Satellite_V02_base'] = sat_path
        os.environ['Satellite_V02_bddm'] = os.path.join(sat_path, 'BULK', 'MATERIAUX')
        os.environ['Satellite_V02_BULK'] = os.path.join(sat_path, 'BULK')
        os.environ['Satellite_V02_INCLUDE'] = os.path.join(sat_path, 'INCLUDE')
        bdf_filename = os.path.join(sat_path, 'JOBS', 'QS', 'satellite_V02_ACA_QS_SOL101_VarEnv.dat')
        read_bdf(bdf_filename, debug=False)

    def test_two_envs(self):
        """fails for two environment variables"""
        sat_path = os.path.abspath(os.path.join(MODEL_PATH, 'Satellite_V02'))
        os.environ['Satellite_V02_base'] = sat_path
        os.environ['Satellite_V02_bddm'] = os.path.join('BULK', 'MATERIAUX')
        #os.environ['Satellite_V02_INCLUDE'] = os.path.join(sat_path, 'INCLUDE')
        #os.environ['Satellite_V02_BULK'] = os.path.join(sat_path, 'BULK')

        pth = "INCLUDE 'Satellite_V02_base:Satellite_V02_bddm:Satellite_V02_Materiaux.blk'"
        with self.assertRaises(SyntaxError):
            pth2 = get_include_filename([pth], include_dir=r'C:\dir\dir2', is_windows=True)
        with self.assertRaises(SyntaxError):
            pth2 = get_include_filename([pth], include_dir=r'C:\dir\dir2', is_windows=False)

        #print('Path:\nold:  %r\nnew:  %r' % (pth, pth2))

    def test_dollar_envs(self):
        """tests sane environment variables"""
        sat_path = os.path.abspath(os.path.join(MODEL_PATH, 'Satellite_V02'))
        os.environ['Satellite_V02_base'] = sat_path
        os.environ['Satellite_V02_bddm'] = os.path.join('BULK', 'MATERIAUX')
        #os.environ['Satellite_V02_INCLUDE'] = os.path.join(sat_path, 'INCLUDE')
        #os.environ['Satellite_V02_BULK'] = os.path.join(sat_path, 'BULK')

        pth = "INCLUDE '%Satellite_V02_bddm%:Satellite_V02_Materiaux.blk'"
        pth2 = get_include_filename([pth], include_dir='', is_windows=True)
        #print(pth2)

        #pth = "INCLUDE '$Satellite_V02_bddm:Satellite_V02_Materiaux.blk'"
        #pth2 = get_include_filename([pth], include_dir='', is_windows=False)
        #print(pth2)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
