import os
import sys
import unittest
import pyNastran
testPath = pyNastran.__path__[0]
#print "testPath = ",testPath

from pyNastran.op2.test.test_op2 import runOP2
from pyNastran.bdf.test.bdf_unitTests import Tester


class OP2_Test(Tester):
    def test_op2_01(self):
        op2Filename = 'solidBending.op2'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        makeGeom = True
        writeBDF = True
        writeF06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        runOP2(op2file, makeGeom=makeGeom, writeBDF=writeBDF, iSubcases=[],
               writeF06=writeF06, debug=debug, stopOnFailure=True)

    def test_op2_02(self):
        op2Filename = 'plate_py.op2'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        makeGeom = True
        writeBDF = True
        writeF06 = True
        debug = False
        op2file = os.path.join(folder, op2Filename)
        runOP2(op2file, makeGeom=makeGeom, writeBDF=writeBDF, iSubcases=[],
               writeF06=writeF06, debug=debug, stopOnFailure=True)

if __name__ == '__main__':
    unittest.main()
