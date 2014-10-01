import os
import unittest
from numpy import array, allclose

import pyNastran
from pyNastran.bdf.bdf import BDF
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')

log = None
class TestLoadSum(unittest.TestCase):
    def test_loads_sum_01(self):
        model = BDF(log=log)
        bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
        model.read_bdf(bdf_filename)
        loadcase_id = 1
        #print "keys1", model.loads.keys()

        p0 = array([0., 0., 0.])
        F_expected = array([23000., 0., 0.])
        M_expected = array([0., 33209.869, -22803.951])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

    def test_loads_sum_02(self):
        model = BDF(log=log)
        bdf_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        model.read_bdf(bdf_filename)
        loadcase_id = 10000
        #print "keys2", model.loads.keys()

        p0 = array([0., 0., 0.])
        F_expected = array([0., 0., 10000.])
        M_expected = array([5000., -5000., 0.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        loadcase_id = 123458
        p0 = array([0., 0., 0.])
        F_expected = array([0., 0., 10000.])
        M_expected = array([5000., -5000., 0.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

    def test_loads_sum_03(self):
        if 0:
            model = BDF(log=log)
            bdf_filename = os.path.join(model_path, 'iSat', 'ISat_Launch_Sm_4pt.dat')
            model.read_bdf(bdf_filename)
            loadcase_id = 1
            #print "keys3", model.loads.keys()

            p0 = array([0., 0., 0.])
            F_expected = array([0., 0., 1.])
            M_expected = array([0., 0., 0.])
            F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

            self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
            self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

    def test_loads_sum_04(self):
        p0 = array([0., 0., 0.])
        model = BDF(log=log)
        bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
        model.read_bdf(bdf_filename)
        #print "keys4", model.loads.keys()

        loadcase_id = 1
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        loadcase_id = 2
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        #---------
        loadcase_id = 3
        A = 0.
        for e, element in model.elements.iteritems():
            A += element.Area()
        A_expected = 100.
        self.assertTrue(allclose(A, A_expected), 'loadcase_id=%s A_expected=%s A=%s' % (loadcase_id, A_expected, A))
        p = 3.
        Fi = p * A
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(p*A, F[2]), 'loadcase_id=%s p*A=%s F=%s' % (loadcase_id, p*A, F))

        F_expected = array([0.,        0., 300.])
        M_expected = array([1500., -1500.,   0.])
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        #---
        loadcase_id = 10
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        #---
        loadcase_id = 4
        F_expected = array([0.,        0., 300.])
        M_expected = array([1500., -1500.,   0.])
        F_expected *= 5.
        M_expected *= 5.
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        #print 'F =', F
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        loadcase_id = 5
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        F_expected = array([0.,        0., 300.])
        M_expected = array([1500., -1500.,   0.])
        F_expected *= 7.
        M_expected *= 7.
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        loadcase_id = 6
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        F_expected = array([0.,        0., 300.])
        M_expected = array([1500., -1500.,   0.])
        F_expected *= 7. * 5.
        M_expected *= 7. * 5.
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

        #---------
        loadcase_id = 11
        A_expected = 4.
        A = 4.
        p = 3.
        Fi = p * A
        element = model.elements[1]
        normal = element.Normal()
        normal_expected = array([0., 0., 1.])
        self.assertTrue(allclose(normal_expected, normal), 'loadcase_id=%s normal_expected=%s normal=%s' % (loadcase_id, normal_expected, normal))
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(p*A, F[2]), 'loadcase_id=%s p*A=%s F=%s' % (loadcase_id, p*A, F))

        F_expected = array([0.,    0., 12.])
        M_expected = array([12., -12.,  0.])
        self.assertTrue(allclose(F_expected, F), 'loadcase_id=%s F_expected=%s F=%s' % (loadcase_id, F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'loadcase_id=%s M_expected=%s M=%s' % (loadcase_id, M_expected, M))

if __name__ == '__main__':
    unittest.main()