import os
import unittest
from io import StringIO
from numpy import array, allclose, cross
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf import GRID
from pyNastran.bdf.mesh_utils.loads import sum_forces_moments, sum_forces_moments_elements
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')


log = None
class TestLoadSum(unittest.TestCase):
    def test_loads_sum_01(self):
        """tests FORCE"""
        model = BDF(log=log, debug=False)
        bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
        model.read_bdf(bdf_filename)
        loadcase_id = 1
        #print("keys1", model.loads.keys())

        p0 = array([0., 0., 0.])
        F_expected = array([23000., 0., 0.])
        M_expected = array([0., 33209.869, -22803.951])
        eids = None
        nids = None
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

    def test_loads_sum_02(self):
        """tests FORCE"""
        model = BDF(log=log, debug=False)
        bdf_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        model.read_bdf(bdf_filename)
        loadcase_id = 10000
        #print("keys2", model.loads.keys())

        p0 = array([0., 0., 0.])
        F_expected = array([0., 0., 10000.])
        M_expected = array([5000., -5000., 0.])
        eids = None
        nids = None
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        loadcase_id = 123458
        p0 = array([0., 0., 0.])
        F_expected = array([0., 0., 10000.])
        M_expected = array([5000., -5000., 0.])
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

    def test_loads_sum_03(self):
        """tests N/A"""
        if 0:  # pragma: no cover
            model = BDF(log=log, debug=False)
            bdf_filename = os.path.join(model_path, 'iSat', 'ISat_Launch_Sm_4pt.dat')
            model.read_bdf(bdf_filename)
            loadcase_id = 1
            #print("keys3", model.loads.keys())

            p0 = array([0., 0., 0.])
            F_expected = array([0., 0., 1.])
            M_expected = array([0., 0., 0.])
            eids = None
            nids = None
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

            self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
            self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

    def test_loads_sum_04(self):
        """
        tests:
          - 1=FORCE
          - 2=LOAD/FORCE
          - 3=LOAD/PLOAD4
          - 4=LOAD/PLOAD4
          - 5=LOAD/PLOAD4
          - 6=LOAD/PLOAD4
          - 10=PLOAD4
          - 11=PLOAD4
        """
        p0 = array([0., 0., 0.])
        model = BDF(log=log, debug=False)
        bdf_filename = os.path.join(model_path, 'plate', 'plate.bdf')
        #print(bdf_filename)
        model.read_bdf(bdf_filename)
        #print("keys4", model.loads.keys())

        loadcase_id = 1
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        eids = None
        nids = None
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        loadcase_id = 2
        F_expected = array([600., 0., 0.])
        M_expected = array([0., 0., -3000.])
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        #---------
        loadcase_id = 3
        A = 0.
        for e, element in model.elements.items():
            A += element.Area()
        A_expected = 100.
        self.assertTrue(allclose(A, A_expected), 'loadcase_id=%s A_expected=%s A=%s' % (loadcase_id, A_expected, A))
        p = 3.
        #Fi = p * A

        eids = None
        nids = None
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        self.assertTrue(allclose(p*A, F1[2]), 'loadcase_id=%s p*A=%s F1=%s' % (loadcase_id, p*A, F1))

        F_expected = array([0., 0., 300.])
        M_expected = array([1500., -1500., 0.])
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        #---
        loadcase_id = 10
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        #---
        loadcase_id = 4
        F_expected = array([0., 0., 300.])
        M_expected = array([1500., -1500., 0.])
        F_expected *= 5.
        M_expected *= 5.
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        #print('F =', F)
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        loadcase_id = 5
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        F_expected = array([0., 0., 300.])
        M_expected = array([1500., -1500., 0.])
        F_expected *= 7.
        M_expected *= 7.
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        loadcase_id = 6
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        F_expected = array([0., 0., 300.])
        M_expected = array([1500., -1500., 0.])
        F_expected *= 7. * 5.
        M_expected *= 7. * 5.
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

        #---------
        loadcase_id = 11
        A_expected = 4.
        A = 4.
        p = 3.
        #Fi = p * A
        element = model.elements[1]
        normal = element.Normal()
        normal_expected = array([0., 0., 1.])
        self.assertTrue(allclose(normal_expected, normal), 'loadcase_id=%s normal_expected=%s normal=%s' % (loadcase_id, normal_expected, normal))
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        self.assertTrue(allclose(p*A, F1[2]), 'loadcase_id=%s p*A=%s F1=%s' % (loadcase_id, p*A, F1))

        F_expected = array([0., 0., 12.])
        M_expected = array([12., -12., 0.])
        self.assertTrue(allclose(F_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F_expected, F1))
        self.assertTrue(allclose(M_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M_expected, M1))

    def test_loads_sum_05(self):
        """
        tests:
         - 1=LOAD/PLOAD4
         - 2=LOAD/PLOAD4/FORCE
         - 5=PLOAD4
         - 6=PLOAD4
         - 1001=PLOAD4
         - 1002=1002
         - 1003=PLOAD
        """
        model = BDF(log=log, debug=False)
        bdf_filename = os.path.join(model_path, 'real', 'loads', 'loads.bdf')
        model.read_bdf(bdf_filename)

        p = 3.
        A = 1.
        n = array([0., 0., 1.])
        F1001_expected = p * A * n
        r = array([0.5, 1.5, 0.])
        p0 = array([0., 0., 0.])
        M1001_expected = cross(r, F1001_expected)

        loadcase_id = 1001
        eids = None
        nids = None
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F1001_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1001_expected, F1))
        self.assertTrue(allclose(M1001_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1001_expected, M1))

        loadcase_id = 1002
        r = array([4., 2., 0.])
        F1002_expected = array([0., 0., 1.])
        M1002_expected = cross(r, F1002_expected)
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        self.assertTrue(allclose(F1002_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1002_expected, F1))
        self.assertTrue(allclose(M1002_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1002_expected, M1))

        loadcase_id = 1
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F1001_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1001_expected, F1))
        self.assertTrue(allclose(M1001_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1001_expected, M1))

        loadcase_id = 2
        F2_expected = F1001_expected + F1002_expected
        M2_expected = M1001_expected + M1002_expected
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F2_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F2_expected, F1))
        self.assertTrue(allclose(M2_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M2_expected, M1))


        F6_expected = 2. * (3. * F1001_expected + 13. * F1002_expected)
        M6_expected = 2. * (3. * M1001_expected + 13. * M1002_expected)

        F7_expected = 7. * 11. * F6_expected
        M7_expected = 7. * 11. * M6_expected
        if 0:  # pragma: no cover
            loadcase_id = 6
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F6_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F6_expected, F1))
            self.assertTrue(allclose(M6_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M6_expected, M1))

            loadcase_id = 7
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F7_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F7_expected, F1))
            self.assertTrue(allclose(M7_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M7_expected, M1))

        loadcase_id = 5
        p = 2.
        A = 1.
        n = array([0., 1., 1.]) / np.sqrt(2.)
        F5_expected = p * A * n
        r = array([0.5, 0.5, 0.])
        M5_expected = cross(r, F5_expected)
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F5_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F5_expected, F1))
        self.assertTrue(allclose(M5_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M5_expected, M1))
        #print('loadcase_id=%s F1=%s M1=%s' % (loadcase_id, F1, M1))

        loadcase_id = 6
        p = 2.
        A = 1.
        n = array([0., 0., 0.5]) / 0.5
        F6_expected = p * A * n
        r = array([0.5, 0.5, 0.])
        M6_expected = cross(r, F6_expected)
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F6_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F6_expected, F1))
        self.assertTrue(allclose(M6_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M6_expected, M1))
        #print('loadcase_id=%s F1=%s M1=%s' % (loadcase_id, F1, M1))

        loadcase_id = 1003
        p = 9.
        A = 1.
        n = array([0., 0., 1.])
        F1003_expected = p * A * n
        r = array([0.5, 0.5, 0.])
        M1003_expected = cross(r, F1003_expected)
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        self.assertTrue(allclose(F1003_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1003_expected, F1))
        self.assertTrue(allclose(M1003_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1003_expected, M1))

        loadcase_id = 8
        F8_expected = 2. * (3. * F7_expected + 2. * F1003_expected)
        M8_expected = 2. * (3. * M7_expected + 2. * M1003_expected)
        if 0:  # pragma: no cover
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F8_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F8_expected, F1))
            self.assertTrue(allclose(M8_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M8_expected, M1))

        loadcase_id = 800
        p = 3.5
        A = 1.
        n = array([0., 0., 1.])
        F800_expected = p * A * n
        r = array([3.5, 1.5, 0.])
        M800_expected = cross(r, F800_expected)
        if 0:  # pragma: no cover
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F800_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F800_expected, F1))
            self.assertTrue(allclose(M800_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M800_expected, M1))

        loadcase_id = 801
        F801_expected = F800_expected
        M801_expected = M800_expected
        if 0:  # pragma: no cover
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F801_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F801_expected, F1))
            self.assertTrue(allclose(M801_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M801_expected, M1))

        loadcase_id = 802
        p = 3.5
        A = 0.5
        n = array([0., 0., 1.])
        F802_expected = p * A * n
        rx = (3. + 4. + 4.) / 3.
        ry = (1. + 1. + 2.) / 3.
        r = array([rx, ry, 0.])
        M802_expected = cross(r, F802_expected)
        if 0:  # pragma: no cover
            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F802_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F802_expected, F1))
            self.assertTrue(allclose(M802_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M802_expected, M1))
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)
        model.write_bdf(bdf_file, size=16)

    def _test_loads_sum_06(self):
        model = BDF(log=log, debug=False)
        bdf_filename = os.path.join(model_path, 'real', 'loads', 'bars.bdf')
        model.read_bdf(bdf_filename)
        p0 = array([0., 0., 0.])

        loadcase_id = 1
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        if 0:  # pragma: no cover
            r = array([0., 0., 0.])
            F1_expected = array([0., 0., 1.])
            M1_expected = cross(r, F1_expected)

            F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False)
            F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids, include_grav=False)
            assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
            assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
            self.assertTrue(allclose(F1_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1_expected, F1))
            self.assertTrue(allclose(M1_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1_expected, M1))

    def test_loads_sum_radial_01(self):
        model = BDF(debug=False)
        model.nodes[1] = GRID(1, cp=1, xyz=[0., 0., 0.], cd=0, ps='', seid=0,
                              comment='')
        cid = 1
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xaxis = [1., 0., 0.]
        model.add_cord2c(cid, rid=0, origin=origin, zaxis=zaxis, xzplane=xaxis,
                         comment='')

        sid = 1
        node = 1
        cid = 1
        mag = 1.1
        xyz = [1., 0., 0.]
        unused_radial_force = model.add_force(sid, node, mag, xyz, cid=cid, comment='')

        sid = 2
        xyz = [1., 90., 0.]
        mag = 2.2
        unused_theta_force = model.add_force(sid, node, mag, xyz, cid=cid, comment='')
        model.cross_reference()

        p0 = 1
        eids = None
        nids = None

        loadcase_id = 1
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False,
                                    xyz_cid0=None)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                             include_grav=False,
                                             xyz_cid0=None)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)

        F1_expected = np.array([1.1, 0., 0.])
        M1_expected = np.array([0., 0., 0.])
        self.assertTrue(allclose(F1_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F1_expected, F1))
        self.assertTrue(allclose(M1_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M1_expected, M1))

        loadcase_id = 2
        F1, M1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False,
                                    xyz_cid0=None)
        F2, M2 = sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                             include_grav=False,
                                             xyz_cid0=None)
        assert np.allclose(F1, F2), 'F1=%s F2=%s' % (F1, F2)
        assert np.allclose(M1, M2), 'M1=%s M2=%s' % (M1, M2)
        F2_expected = np.array([0., 2.2, 0.])
        M2_expected = np.array([0., 0., 0.])
        self.assertTrue(allclose(F2_expected, F1), 'loadcase_id=%s F_expected=%s F1=%s' % (loadcase_id, F2_expected, F1))
        self.assertTrue(allclose(M2_expected, M1), 'loadcase_id=%s M_expected=%s M1=%s' % (loadcase_id, M2_expected, M1))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
