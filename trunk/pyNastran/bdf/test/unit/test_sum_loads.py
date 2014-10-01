import os
import unittest
from numpy import array, array_equal, allclose

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

        self.assertTrue(allclose(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))

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

        self.assertTrue(allclose(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))


        loadcase_id = 123458
        p0 = array([0., 0., 0.])
        F_expected = array([0., 0., 10000.])
        M_expected = array([5000., -5000., 0.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

        self.assertTrue(allclose(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
        self.assertTrue(allclose(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))

    def test_loads_sum_03(self):
        if 0:
            model = BDF(log=log)
            bdf_filename = os.path.join(model_path, 'iSat', 'ISat_Launch_Sm_4pt.dat')
            model.read_bdf(bdf_filename)
            loadcase_id = 1
            print "keys3", model.loads.keys()

            p0 = array([0., 0., 0.])
            F_expected = array([0., 0., 1.])
            M_expected = array([0., 0., 0.])
            F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)

            self.assertTrue(array_equal(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
            self.assertTrue(array_equal(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))

    def test_loads_sum_04(self):
        model = BDF(log=log)
        bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
        model.read_bdf(bdf_filename)
        #print "keys4", model.loads.keys()

        p0 = array([0., 0., 0.])

        loadcase_id = 1
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        self.assertTrue(array_equal(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
        self.assertTrue(array_equal(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))

        loadcase_id = 2
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        F, M = model.sum_forces_moments(p0, loadcase_id, include_grav=False)
        self.assertTrue(array_equal(F_expected, F), 'F_expected=%s F=%s' % (F_expected, F))
        self.assertTrue(array_equal(M_expected, M), 'M_expected=%s M=%s' % (M_expected, M))

if __name__ == '__main__':
    unittest.main()