import os
import unittest
from math import sqrt

import pyNastran
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2
testpath = os.path.join(pyNastran.__path__[0], 'f06', 'test', 'tests')
print "testpath =", testpath

class TestF06(unittest.TestCase):
    def runF06(self, f06_name, op2_name=None):
        f06 = F06(f06_name)
        f06.readF06()
        f06.writeF06(f06_name+'.out')
        
        if op2_name:
            op2 = OP2(op2_name)
            op2.readOP2()
            op2.writeF06(op2_name+'.out')
            return f06, op2
        
        return f06
    
    def test_plate_vonmises(self):
        f06name = os.path.join(testpath, 'plate.f06')
        op2name = os.path.join(testpath, 'plate.op2')
        print("f06name = ", f06name)
        f06, op2 = self.runF06(f06name, op2name)
        
        for (loadcase, stress) in f06.plateStress.iteritems():
            print "%3s %3s %6s %8s" % ('EID', 'NID', 'iLayer', 'VM_Stress')

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
                            print "%8s %8s %6s %8s" % (eid, nid, ilayer, ovmii)
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
                            ovmii = sqrt(o1**2 - o1*o2 + o2**2 + 3*t12**2)
                            print "%3s %3s %6s %8s" % (eid, nid, ilayer, ovmii)
        self.assertEquals(op2.plateStress[1].ovmShear[25]['C'][0], 276.8023376464844)
        self.assertEquals(f06.plateStress[1].ovmShear[25]['C'][0], 276.8023)

if __name__ == '__main__':
    unittest.main()
