from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
import pyNastran
testPath = pyNastran.__path__[0]
#print("testPath = %s" %(testPath))
from pyNastran.bdf.test.test_bdf import runBDF, run_all_files_in_folder


class Tester(unittest.TestCase):

    def runBDF(self, folder, bdfFilename, xref=False, cid=None,
               meshForm='combined', debug=False):
        cid = 0
        #xref = False
        runBDF(folder, bdfFilename, xref=xref, cid=cid, isFolder=True,
               meshForm=meshForm, debug=debug)

    def runAllFilesInFolder(self, folder, xref=False, cid=None, debug=False):
        run_all_files_in_folder(folder, xref=xref, cid=cid, debug=debug)


class BDF_Test(Tester):
    def test_bdf_01(self):
        bdfFilename = 'solidBending.bdf'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        self.runBDF(folder, bdfFilename, xref=True)

    def test_bdf_02(self):
        bdfFilename = 'plate_py.dat'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        self.runBDF(folder, bdfFilename, xref=True)

    def test_bdf_03(self):
        bdfFilename = 'beam_modes.dat'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        #self.runBDF(folder, bdfFilename, xref=True) ## PBEAML is not supported

    def test_bdf_04(self):
        bdfFilename = 'testA.bdf'
        folder = os.path.abspath(os.path.join(testPath, 'bdf', 'test', 'unit'))
        self.runBDF(folder, bdfFilename)
        #self.runBDF(folder, bdfFilename, xref=True) ## PBEAML is not supported

if __name__ == '__main__':
    unittest.main()
