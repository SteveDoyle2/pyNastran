import os
import unittest
import pyNastran
testPath = pyNastran.__path__[0]
#print "testPath = ",testPath

from pyNastran.op2.test.test_op2 import run_op2
from pyNastran.bdf.test.bdf_unitTests import Tester


class OP2_Test(Tester):
    def test_op2_01(self):
        op2Filename = os.path.join('solid_bending', 'solid_bending.op2')
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        makeGeom = True
        writeBDF = True
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        run_op2(op2file, makeGeom=makeGeom, writeBDF=writeBDF, iSubcases=[],
                write_f06=write_f06, debug=debug, stopOnFailure=True)

    def test_op2_02(self):
        op2Filename = os.path.join('plate_py', 'plate_py.op2')
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        makeGeom = True
        writeBDF = True
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        run_op2(op2file, makeGeom=makeGeom, writeBDF=writeBDF, iSubcases=[],
                write_f06=write_f06, debug=debug, stopOnFailure=True)

if __name__ == '__main__':
    unittest.main()
