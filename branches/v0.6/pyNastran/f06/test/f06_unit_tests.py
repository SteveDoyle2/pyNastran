import os
import unittest
from math import sqrt

from numpy import array, array_equiv

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06, FatalError
from pyNastran.op2.op2 import OP2
from pyNastran.op4.op4 import OP4

testpath = os.path.join(pyNastran.__path__[0], '..', 'models')
#print "testpath =", testpath

class TestF06(unittest.TestCase):
    def run_model(self, bdf_name=None, f06_name=None, op2_name=None,
                  op4_name=None, dynamic_vars=None, f06_has_weight=True):
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
            if f06_has_weight:
                assert f06.grid_point_weight.reference_point is not None
            else:
                assert f06.grid_point_weight.reference_point is None


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

    def test_blade2dv_fatal_1(self):
        f06_name = os.path.join(testpath, 'blade_2dv', 'blade_2dv.f06_fatal')
        f06 = F06(f06_name, debug=False, log=None)
        with self.assertRaises(FatalError):
            f06.readF06()

    def test_blade2dv_fatal_2(self):
        f06_name = os.path.join(testpath, 'blade_2dv', 'blade_2dv.f06_fatal')
        bdf_name = os.path.join(testpath, 'blade_2dv', 'blade_2dv.bdf')
        #bdf2 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)
        #self.assertEquals(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)

        f06 = F06(f06_name, debug=False, log=None)

        try:  # this is supposed to raise a FatalError
            f06.read_f06()
            assert False, 'a FATAL should have been raised'
        except FatalError:
            pass

        f06.write_f06(f06_name+'.out')

        ref_point = f06.grid_point_weight.reference_point
        #print "ref_point", ref_point
        ref_point_exact = 0
        self.assertEqual(ref_point, ref_point_exact)

        MO = f06.grid_point_weight.MO
        #print "MO", MO
        MO_exact = array(
            [[  1.22085800e-01,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  5.33146300e-01,  -1.22767700e-05],
             [  0.00000000e+00,   1.22085800e-01,   0.00000000e+00,  -5.33146300e-01,  0.00000000e+00,   1.57186600e-01],
             [  0.00000000e+00,   0.00000000e+00,   1.22085800e-01,   1.22767700e-05, -1.57186600e-01,   0.00000000e+00],
             [  0.00000000e+00,  -5.33146300e-01,   1.22767700e-05,   3.20227600e+00, -7.13340800e-06,  -6.83890800e-01],
             [  5.33146300e-01,   0.00000000e+00,  -1.57186600e-01,  -7.13340800e-06,  3.45033400e+00,  -7.35886500e-05],
             [ -1.22767700e-05,   1.57186600e-01,   0.00000000e+00,  -6.83890800e-01, -7.35886500e-05,   2.50287600e-01]])


        S = f06.grid_point_weight.S
        S_exact = array([[ 1.,  0.,  0.],
                         [ 0.,  1.,  0.],
                         [ 0.,  0.,  1.]])
        #print "S", S

        mass = f06.grid_point_weight.mass
        mass_exact = array([ 0.1220858,  0.1220858,  0.1220858])
        self.assertTrue(array_equiv(mass, mass_exact))
        #print "mass =", mass

        cg = f06.grid_point_weight.cg
        cg_exact = array(
            [[  0.00000000e+00,   1.00558600e-04,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   0.00000000e+00,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   1.00558600e-04,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00]])
        #print "cg =", cg
        self.assertTrue(array_equiv(cg, cg_exact))

        IS = f06.grid_point_weight.IS
        #print "IS", IS
        IS_exact = array([[  8.74036600e-01,  -8.67305300e-06,  -2.54028500e-03],
                          [ -8.67305300e-06,   9.19714300e-01,   1.99762300e-05],
                          [ -2.54028500e-03,   1.99762300e-05,   4.79082500e-02]])

        IQ = f06.grid_point_weight.IQ
        #print "IQ", IQ
        IQ_exact = array([[ 0.04790044,  0.,          0.        ],
                          [ 0.,          0.9197143,   0.        ],
                          [ 0.,          0.,          0.8740444 ]])
        self.assertTrue(array_equiv(IQ, IQ_exact))

    def test_blade2dv_fatal_3(self):
        f06_name = os.path.join(testpath, 'blade_2dv', 'blade_2dv.f06_fatal')
        bdf_name = os.path.join(testpath, 'blade_2dv', 'blade_2dv.bdf')
        #bdf2 = self.run_model(bdfname2, dynamic_vars=dynamic_vars)
        #self.assertEquals(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)

        f06 = F06(f06_name, debug=False, log=None)

        # we skip the fatal by stopping after reading the matrices
        f06.stop_after_reading_grid_point_weight()
        f06.read_f06()

        f06.write_f06(f06_name+'.out')

        ref_point = f06.grid_point_weight.reference_point
        #print "ref_point", ref_point
        ref_point_exact = 0
        self.assertEqual(ref_point, ref_point_exact)

        MO = f06.grid_point_weight.MO
        #print "MO", MO
        MO_exact = array(
            [[  1.22085800e-01,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  5.33146300e-01,  -1.22767700e-05],
             [  0.00000000e+00,   1.22085800e-01,   0.00000000e+00,  -5.33146300e-01,  0.00000000e+00,   1.57186600e-01],
             [  0.00000000e+00,   0.00000000e+00,   1.22085800e-01,   1.22767700e-05, -1.57186600e-01,   0.00000000e+00],
             [  0.00000000e+00,  -5.33146300e-01,   1.22767700e-05,   3.20227600e+00, -7.13340800e-06,  -6.83890800e-01],
             [  5.33146300e-01,   0.00000000e+00,  -1.57186600e-01,  -7.13340800e-06,  3.45033400e+00,  -7.35886500e-05],
             [ -1.22767700e-05,   1.57186600e-01,   0.00000000e+00,  -6.83890800e-01, -7.35886500e-05,   2.50287600e-01]])


        S = f06.grid_point_weight.S
        S_exact = array([[ 1.,  0.,  0.],
                         [ 0.,  1.,  0.],
                         [ 0.,  0.,  1.]])
        #print "S", S

        mass = f06.grid_point_weight.mass
        mass_exact = array([ 0.1220858,  0.1220858,  0.1220858])
        self.assertTrue(array_equiv(mass, mass_exact))
        #print "mass =", mass

        cg = f06.grid_point_weight.cg
        cg_exact = array(
            [[  0.00000000e+00,   1.00558600e-04,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   0.00000000e+00,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   1.00558600e-04,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00]])
        #print "cg =", cg
        self.assertTrue(array_equiv(cg, cg_exact))

        IS = f06.grid_point_weight.IS
        #print "IS", IS
        IS_exact = array([[  8.74036600e-01,  -8.67305300e-06,  -2.54028500e-03],
                          [ -8.67305300e-06,   9.19714300e-01,   1.99762300e-05],
                          [ -2.54028500e-03,   1.99762300e-05,   4.79082500e-02]])

        IQ = f06.grid_point_weight.IQ
        #print "IQ", IQ
        IQ_exact = array([[ 0.04790044,  0.,          0.        ],
                          [ 0.,          0.9197143,   0.        ],
                          [ 0.,          0.,          0.8740444 ]])
        self.assertTrue(array_equiv(IQ, IQ_exact))

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

        (bdf, f06, op2) = self.run_model(bdfname, f06name, op2name, f06_has_weight=False)
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