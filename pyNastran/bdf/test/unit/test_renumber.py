from __future__ import print_function
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_interface.dev_utils import bdf_renumber
from pyNastran.utils.dev import get_files_of_type
import unittest
import os

import pyNastran
pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, '..', 'models')
test_path = os.path.join(pkg_path, 'bdf', 'test')
unit_path = os.path.join(test_path, 'unit')

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
        msg += '$[RBE3, self.eid, None, self.refgrid, self.refc]\n'
        msg += 'RBE3       12225           33     123456      1.     123    34      36\n'
        msg += '            20      10\n'

        msg += 'ENDDATA\n'
        with open('renumber_in.bdf', 'w') as f:
            f.write(msg)

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
        bdf_filename = 'renumber_in.bdf'
        bdf_filename_renumber = 'renumber_out.bdf'
        bdf_renumber(bdf_filename, bdf_filename_renumber)

        #model = BDF(debug=False)
        #model.read_bdf(bdf_filename)
        #model.write_bdf(bdf_filename_check)
        model = BDF(debug=False)
        model.read_bdf(bdf_filename_renumber)

    def test_renumber_02(self):
        bdf_filename = os.path.join(model_path, 'iSat', 'ISat_Dploy_Sm.dat')
        bdf_filename_renumber = os.path.join(model_path, 'iSat', 'ISat_Dploy_Sm_renumber.dat')
        bdf_filename_check = os.path.join(model_path, 'iSat', 'ISat_Dploy_Sm_check.dat')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)


    def test_renumber_03(self):
        bdf_filename = os.path.join(model_path, 'cbush', 'cbush.dat')
        bdf_filename_renumber = os.path.join(model_path, 'cbush', 'cbush_renumber.dat')
        bdf_filename_check = os.path.join(model_path, 'cbush', 'cbush_check.dat')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)

    def test_renumber_04(self):
        bdf_filename = os.path.join(model_path, 'complex', 'tet10', 'Simple_Example.bdf')
        bdf_filename_renumber = os.path.join(model_path, 'complex', 'tet10', 'Simple_Example_renumber.bdf')
        bdf_filename_check = os.path.join(model_path, 'complex', 'tet10', 'Simple_Example_check.bdf')
        check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)

    def test_renumber_05(self):
        dirname = os.path.join(unit_path, 'obscure')
        bdf_filenames = get_files_of_type(dirname, extension='.bdf')
        for bdf_filename in bdf_filenames:
            print('bdf_filename = %s' % (bdf_filename))
            basename = os.path.basename(bdf_filename)
            base, ext = os.path.splitext(basename)
            bdf_filename_renumber = os.path.join(dirname, base + '_renumber.bdf_test')
            bdf_filename_check = os.path.join(dirname, base + '_check.bdf_test')
            check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check)


def check_renumber(bdf_filename, bdf_filename_renumber, bdf_filename_check):
    bdf_renumber(bdf_filename, bdf_filename_renumber)

    model = BDF(debug=False)
    model.read_bdf(bdf_filename)
    model.write_bdf(bdf_filename_check, interspersed=False)

    model = BDF(debug=False)
    model.read_bdf(bdf_filename_renumber)

if __name__ == '__main__':
    unittest.main()
