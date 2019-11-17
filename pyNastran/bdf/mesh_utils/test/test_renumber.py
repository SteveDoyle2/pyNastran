"""tests bdf_renumber"""
import os
import unittest
from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
#from pyNastran.utils.dev import get_files_of_type

import pyNastran
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'test')
UNIT_PATH = os.path.join(TEST_PATH, 'unit')

class TestRenumber(unittest.TestCase):
    def test_renumber_01(self):
        msg = 'CEND\n'
        msg += 'BEGIN BULK\n'
        msg += 'GRID,10,,1.0,1.0\n'
        msg += 'GRID,30,,3.0,2.0\n'
        msg += 'GRID,20,,2.0,3.0\n'
        msg += 'GRID,33,,3.3,4.0\n'
        msg += 'GRID,34,,3.4,5.0\n'
        msg += 'GRID,35,,3.5,6.0\n'
        msg += 'GRID,36,,3.6,7.0\n'
        msg += 'SPOINT,4,THRU,8\n'
        msg += 'SPOINT,11\n'
        msg += 'CTRIA3,10,8,30,20,10\n'
        msg += 'PSHELL,8,4,0.1\n'
        msg += 'MAT1,4,3.0e7,,0.3\n'

        msg += 'MPC,10,20,1,1.0,10,2,1.0\n'
        msg += 'SPC,2,30,3,-2.6\n'
        msg += 'SPC1,313,12456,33,THRU,34\n'
        msg += 'SPC,314,30,3,-2.6,36,3,-2.6\n'

        msg += '$SPCD,SID,G1,C1, D1,  G2,C2,D2\n'
        msg += 'SPCD, 100,33,436,-2.6,10, 2,.9\n'
        msg += 'SPCD, 101,34,436,-2.6\n'

        msg += '$RBAR, EID, GA, GB, CNA\n'
        msg += 'RBAR,    5, 10, 20, 123456\n'

        msg += '$RBAR1, EID, GA, GB, CB, ALPHA\n'
        msg += 'RBAR1,    9, 20, 10, 123, 6.5-7\n'

        msg += 'RBE1        1001    33    123456\n'
        msg += '              UM    20       123    35       123    34       123\n'
        msg += '                    10       123\n'
        msg += '$[RBE3, eid, None, refgrid, refc]\n'
        msg += 'RBE3       12225           33     123456      1.     123    34      36\n'
        msg += '            20      10\n'
        msg += 'ENDDATA\n'

        bdf_filename = 'renumber_in.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        msg_expected = 'CEND\n'
        msg_expected += 'BEGIN BULK\n'
        msg_expected += 'GRID,6\n'
        msg_expected += 'GRID,7\n'
        msg_expected += 'GRID,8\n'
        msg_expected += 'SPOINT,1,THRU,5\n'
        msg_expected += 'CTRIA3,1,1,8,7,6\n'
        msg_expected += 'PSHELL,1,1,0.1\n'
        msg_expected += 'MAT1,1,3.0e7,,0.3\n'
        msg_expected += 'ENDDATA\n'

        # for now we're testing things don't crash
        bdf_filename_renumber = 'renumber_out.bdf'
        bdf_renumber(bdf_filename, bdf_filename_renumber)

        #model = BDF(debug=False)
        #model.read_bdf(bdf_filename)
        #model.write_bdf(bdf_filename_check)
        model = BDF(debug=False)
        model.read_bdf(bdf_filename_renumber)

        os.remove(bdf_filename)
        os.remove(bdf_filename_renumber)

    def test_renumber_02(self):
        bdf_filename = os.path.join(MODEL_PATH, 'iSat', 'ISat_Dploy_Sm.dat')
        bdf_filename_renumber = os.path.join(MODEL_PATH, 'iSat', 'ISat_Dploy_Sm_renumber.dat')
        bdf_filename_check = os.path.join(MODEL_PATH, 'iSat', 'ISat_Dploy_Sm_check.dat')
        log = SimpleLogger('error')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check, log=log)

    def test_renumber_03(self):
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'cbush', 'cbush.dat')
        bdf_filename_renumber = os.path.join(MODEL_PATH, 'unit', 'cbush', 'cbush_renumber.dat')
        bdf_filename_check = os.path.join(MODEL_PATH, 'unit', 'cbush', 'cbush_check.dat')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)

    def test_renumber_04(self):
        dirname = os.path.join(MODEL_PATH, 'complex', 'tet10')
        bdf_filename = os.path.join(dirname, 'Simple_Example.bdf')
        bdf_filename_renumber = os.path.join(dirname, 'Simple_Example_renumber.bdf')
        bdf_filename_check = os.path.join(dirname, 'Simple_Example_check.bdf')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)

    def test_renumber_05(self):
        """renumbers a deck in a couple ways"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        bdf_filename_out1 = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero1.out')
        bdf_filename_out2 = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero2.out')
        bdf_filename_out3 = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero3.out')
        model = bdf_renumber(bdf_filename, bdf_filename_out1, size=8,
                             is_double=False, starting_id_dict=None,
                             round_ids=False, cards_to_skip=None, debug=False)

        model = read_bdf(bdf_filename, log=log)
        bdf_renumber(model, bdf_filename_out2, size=16, is_double=False,
                     starting_id_dict={
                         'eid' : 1000, 'pid':2000, 'mid':3000,
                         'spc_id' : 4000,},
                     round_ids=False, cards_to_skip=None)
        bdf_renumber(bdf_filename, bdf_filename_out3, size=8,
                     is_double=False, starting_id_dict=None,
                     round_ids=True, cards_to_skip=None)
        read_bdf(bdf_filename_out1, log=log)
        read_bdf(bdf_filename_out2, log=log)
        read_bdf(bdf_filename_out3, log=log)

    #def test_renumber_06(self):
        #dirname = os.path.join(UNIT_PATH, 'obscure')
        #bdf_filenames = get_files_of_type(dirname, extension='.bdf')
        #for bdf_filename in bdf_filenames:
            #print('bdf_filename = %s' % (bdf_filename))
            #basename = os.path.basename(bdf_filename)
            #base, ext = os.path.splitext(basename)
            #bdf_filename_renumber = os.path.join(dirname, base + '_renumber.bdf_test')
            #bdf_filename_check = os.path.join(dirname, base + '_check.bdf_test')
            #check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)


def check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check,
                   log=None):
    """renumbers the file, then reloads both it and the renumbered deck"""
    bdf_renumber(bdf_filename, bdf_filename_renumber)

    model = BDF(debug=False, log=log)
    model.read_bdf(bdf_filename)
    model.write_bdf(bdf_filename_check, interspersed=False)

    model = BDF(debug=False, log=log)
    model.read_bdf(bdf_filename_renumber)

    os.remove(bdf_filename_renumber)
    os.remove(bdf_filename_check)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
