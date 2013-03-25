import os
import unittest
from math import sqrt

import pyNastran
from pyNastran.bdf.bdf2 import BDF
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2
testpath = os.path.join(pyNastran.__path__[0], 'f06', 'test', 'tests')
print "testpath =", testpath

class TestF06(unittest.TestCase):
    def run_model(self, bdf_name=None, f06_name=None, op2_name=None,
                  dynamic_vars=None):
        outputs = []
        if bdf_name:
            bdf = BDF(debug=True, log=None)
            if dynamic_vars is not None:
                print('dynamic')
                bdf.set_dynamic_syntax(dynamic_vars)
            bdf.read_bdf(bdf_name)
            bdf.write_bdf(bdf_name+'.out')
            bdf.write_bdf_as_patran(bdf_name+'.out')
            outputs.append(bdf)

        if f06_name:
            f06 = F06(f06_name, debug=False, log=None)
            f06.readF06()
            f06.writeF06(f06_name+'.out')
            outputs.append(f06)

        if op2_name:
            op2 = OP2(op2_name)
            op2.readOP2()
            op2.writeF06(op2_name+'.out')
            outputs.append(op2)

        assert len(outputs) > 0
        if len(outputs) == 1: return outputs[0]
        return outputs
    
    def test_plate_vonmises(self):
        bdfname = os.path.join(testpath, 'plate.bdf')
        bdfname2 = os.path.join(testpath, 'plate_openmdao.bdf')
        f06name = os.path.join(testpath, 'plate.f06')
        op2name = os.path.join(testpath, 'plate.op2')

        (bdf, f06, op2) = self.run_model(bdfname, f06name, op2name)
        
        dynamic_vars = {'t' : 42.}
        bdf2 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)
        self.assertEquals(bdf.properties[1].t,  0.3, 't=%s' % bdf.properties[1].t)
        self.assertEquals(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)
        
        dynamic_vars = {'t' : 42}
        with self.assertRaises(SyntaxError):
            bdf3 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)

        dynamic_vars = {'t' : 'asdddddddf'}
        with self.assertRaises(SyntaxError):
            bdf3 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)

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
                            self.assertAlmostEquals(ovmii, ovmii2)
                            #print("%3s %3s %6s %8s" % (eid, nid, ilayer, ovmii2))
        self.assertEquals(op2.plateStress[1].ovmShear[25]['C'][0], 276.8023376464844)
        self.assertEquals(f06.plateStress[1].ovmShear[25]['C'][0], 276.8023)
        #f06.print_stats()
        #op2.print_stats()

if __name__ == '__main__':
    unittest.main()