## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 

from __future__ import division, print_function
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
    
if __name__=='__main__':
    unittest.main()
