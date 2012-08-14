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
        folder = os.path.abspath(os.path.join(testPath,'..','models'))
        makeGeom = True
        writeBDF = True
        writeF06 = True
        debug = False
        op2file = os.path.join(folder,op2Filename)
        runOP2(op2file,makeGeom=makeGeom,writeBDF=writeBDF,iSubcases=[],writeF06=writeF06,debug=debug,stopOnFailure=True)

    def test_op2_02(self):
        op2Filename = 'plate_py.op2'
        folder = os.path.abspath(os.path.join(testPath,'..','models'))
        makeGeom = True
        writeBDF = True
        writeF06 = True
        debug = False
        op2file = os.path.join(folder,op2Filename)
        runOP2(op2file,makeGeom=makeGeom,writeBDF=writeBDF,iSubcases=[],writeF06=writeF06,debug=debug,stopOnFailure=True)
    
if __name__=='__main__':
    unittest.main()
