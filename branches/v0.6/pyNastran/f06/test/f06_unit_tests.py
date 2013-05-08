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
import unittest
from math import sqrt

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2
from pyNastran.op4.op4 import OP4

testpath = os.path.join(pyNastran.__path__[0], '..', 'models')
#print "testpath =", testpath

class TestF06(unittest.TestCase):
    def run_model(self, bdf_name=None, f06_name=None, op2_name=None,
                  op4_name=None, dynamic_vars=None):
        outputs = []
        if bdf_name:
            bdf = BDF(debug=False, log=None)
            if dynamic_vars is not None:
                print('dynamic')
                bdf.set_dynamic_syntax(dynamic_vars)
            bdf.read_bdf(bdf_name)
            bdf.write_bdf(bdf_name+'.out', interspersed=False)
            bdf.write_bdf(bdf_name+'.out', interspersed=True)
            outputs.append(bdf)

        if f06_name:
            f06 = F06(f06_name, debug=False, log=None)
            f06.read_f06()
            f06.write_f06(f06_name+'.out')
            outputs.append(f06)

        if op2_name:
            op2 = OP2(op2_name, debug=False)
            op2.read_op2()
            op2.write_f06(op2_name+'.out')
            outputs.append(op2)
        
        if op4_name:
            op4 = OP4(op4_name)
            op4.read_op4(op4Name, matrixNames=None, precision='default')

        assert len(outputs) > 0
        if len(outputs) == 1: return outputs[0]
        return outputs
    
    def test_plate_openmdao(self):
        bdfname2 = os.path.join(testpath, 'plate', 'plate_openmdao.bdf')
        dynamic_vars = {'t' : 42.}
        bdf2 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)
        self.assertEquals(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)
        
        dynamic_vars = {'t' : 42}
        with self.assertRaises(SyntaxError):
            bdf3 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)

        dynamic_vars = {'t' : 'asdddddddf'}
        with self.assertRaises(SyntaxError):
            bdf3 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)

    def test_plate_vonmises(self):
        bdfname = os.path.join(testpath, 'plate', 'plate.bdf')
        f06name = os.path.join(testpath, 'plate', 'plate.f06')
        op2name = os.path.join(testpath, 'plate', 'plate.op2')

        (bdf, f06, op2) = self.run_model(bdfname, f06name, op2name)
        self.assertEquals(bdf.properties[1].t,  0.3, 't=%s' % bdf.properties[1].t)

        self.assertEquals(len(bdf.nodes), 36, bdf.nodes)
        self.assertEquals(len(bdf.elements), 25, bdf.elements)
        self.assertEquals(len(bdf.properties), 1, bdf.properties)
        self.assertEquals(len(bdf.materials), 1, bdf.materials)
        self.assertEquals(len(bdf.loads), 2, bdf.loads)  # FORCE, LOAD
        self.assertEquals(len(bdf.params), 2, bdf.params)
        self.assertEquals(bdf.sol, 101, bdf.sol)
        
        for (loadcase, stress) in f06.plateStress.iteritems():
            #print("%3s %3s %6s %8s" % ('EID', 'NID', 'iLayer', 'VM_Stress'))

            # stress is a PlateStressObject
            if stress.isVonMises():
                #vonMises = 'VON MISES'
                for eid,ovm in sorted(stress.ovmShear.iteritems()):
                    ovmkeys = ovm.keys()
                    ovmkeys.remove('C')
                    ovmkeys.sort()
                    ovmkeys = ['C'] + ovmkeys
                    for nid in ovmkeys:
                        ovmi = ovm[nid]
                        for ilayer, ovmii in enumerate(ovmi):
                            #print("%8s %8s %6s %8s" % (eid, nid, ilayer, ovmii))
                            pass
            else:
                #vonMises = 'MAX SHEAR'
                for eid,ovm in sorted(stress.ovmShear.iteritems()):
                    ovmkeys = ovm.keys()
                    ovmkeys.remove('C')
                    ovmkeys.sort()
                    ovmkeys = ['C'] + ovmkeys
                    for nid in ovmkeys:
                        ovmi = ovm[nid]
                        for ilayer, ovmii in enumerate(ovmi):
                            o1 = stress.oxx[eid][nid][ilayer]
                            o2 = stress.oyy[eid][nid][ilayer]
                            t12 = stress.txy[eid][nid][ilayer]
                            ovmii2 = sqrt(o1**2 - o1*o2 + o2**2 + 3*t12**2)
                            self.assertAlmostEqual(ovmii, ovmii2)
                            #print("%3s %3s %6s %8s" % (eid, nid, ilayer, ovmii2))
        self.assertEquals(op2.plateStress[1].ovmShear[25]['C'][0], 276.8023376464844)
        self.assertEquals(f06.plateStress[1].ovmShear[25]['C'][0], 276.8023)
        #f06.print_stats()
        #op2.print_stats()

if __name__ == '__main__':
    unittest.main()