import os
import unittest
from pyNastran.op2.dev.op2 import read_post_op2
from pyNastran.op2.op2 import OP2
import timeit


import pyNastran
test_path = pyNastran.__path__[0]

class TestOP2New(unittest.TestCase):
    def _test_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        op2_filename = os.path.join(folder, 'solid_bending', 'solid_bending.op2')
        op2_filename = os.path.join(folder, 'solid_bending', 'solid_bending.op2')

        #op2 = OP2(op2_filename)
        #op2.rdn2cop2()
        #print(op2.object_attributes())
        #print('dbnames =', op2.dbnames)
        #print('dblist =', op2.dblist)

        o2 = read_post_op2(op2_filename, verbose=True, getougv1=True)
        print('a =', o2.keys())
        print('b =', o2['mats'].keys())
        print('c =', o2['mats']['ougv1'])

    def test_02(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        op2_filename = os.path.join(folder, 'iSat', 'isat_dploy_sm_orig.op2')

        #op2 = OP2(op2_filename)
        #op2.rdn2cop2()
        #print(op2.object_attributes())
        #print('dbnames =', op2.dbnames)
        #print('dblist =', op2.dblist)
        o2 = read_post_op2(op2_filename, verbose=True, getougv1=True)
        #print('a =', o2.keys())
        #print('b =', o2['mats'].keys())
        #print('c =', len(o2['mats']['ougv1']))
        #print('d =', o2['mats']['ougv1'][0].keys())
        #print('e =', o2['mats']['ougv1'][0]['lambda'])
        #print('e.shape =', o2['mats']['ougv1'][0]['lambda'].shape)
        #print('e.shape =', o2['mats']['ougv1'][0]['ougv1'].shape)

        #s = 'o2 = rdpostop2(r"F:\work\pyNastran\pyNastran\master2\models\iSat\isat_dploy_sm.op2", verbose=True, getougv1=True)'
        #timeit.timeit(s, setup="from pyNastran.op2.dev.op2 import rdpostop2", number=10)

    def _test_02b(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        op2_filename = os.path.join(folder, 'iSat', 'isat_dploy_sm.op2')

        #op2 = OP2()
        #op2.read_op2(op2_filename, combine=False)

        s = '\n\nop2 = OP2()\nop2.read_op2(r"F:\work\pyNastran\pyNastran\master2\models\iSat\isat_dploy_sm.op2", combine=False)'
        timeit.timeit(s, setup="from pyNastran.op2.op2 import OP2", number=10)

    def test_03(self):
        op2_filename = os.path.join('sol_101_elements', 'mode_solid_shell_bar.op2')

        # op2 = OP2(op2_filename)
        # op2.rdn2cop2()
        # print(op2.object_attributes())
        #print('dbnames =', op2.dbnames)
        #print('dblist =', op2.dblist)
        o2 = read_post_op2(op2_filename, verbose=True, getougv1=True)

    def test_04(self):
        op2_filename = 'mat_b_dn.op2'
        o2 = read_post_op2(op2_filename, verbose=True, getougv1=True)
        #print('a =', o2.keys())
        #print('b =', o2['mats'].keys())
        #print('c =', o2['mats']['ougv1'][0].keys())
        #print('d =', o2['mats']['ougv1'][0]['lambda'])


if __name__ == '__main__':
    unittest.main()
