import os
import unittest
import pyNastran
testPath = pyNastran.__path__[0]
#print "testPath = ",testPath

from pyNastran.op2.op2 import OP2
from pyNastran.op2.test.test_op2 import run_op2
from pyNastran.bdf.test.bdf_unit_tests import Tester


class TestOP2(Tester):
    def _spike(self):
        op2 = OP2()
        op2.set_results('solidStress.oxx')
        op2.read_op2(op2Filename)

    def test_set_results(self):
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        op2Filename = os.path.join(folder, 'solid_bending', 'solid_bending.op2')
        op2 = OP2()
        op2.set_results('stress')
        op2.read_op2(op2Filename)
        self.assertEqual(len(op2.solidStress), 1), len(op2.solidStress)
        self.assertEqual(len(op2.displacements), 0), len(op2.displacements)

        op2 = OP2()
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2Filename)
        self.assertEqual(len(op2.solidStress), 1), len(op2.solidStress)
        self.assertEqual(len(op2.displacements), 1), len(op2.displacements)

    def test_op2_01(self):
        op2Filename = os.path.join('solid_bending', 'solid_bending.op2')
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        make_geom = True
        write_bdf = True
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, debug=debug, stopOnFailure=True)

    def test_op2_02(self):
        op2Filename = os.path.join('plate_py', 'plate_py.op2')
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        make_geom = True
        write_bdf = True
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, debug=debug, stopOnFailure=True)

if __name__ == '__main__':
    unittest.main()
