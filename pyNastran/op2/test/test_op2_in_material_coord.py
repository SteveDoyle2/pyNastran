import os
import unittest
import numpy as np

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.data_in_material_coord import (data_in_material_coord,
        get_eids_from_op2_vector, force_vectors, stress_vectors,
        strain_vectors)
pkg_path = pyNastran.__path__[0]


class TestMaterialCoord(unittest.TestCase):
    def test_force(self):
        is_failed = False
        bdf = BDF(debug=False)
        op2 = OP2(debug=False)
        basepath = os.path.join(pkg_path, 'op2', 'test', 'test_flat_plate_metallic')
        bdf.read_bdf(os.path.join(basepath, 'flat_plate_metallic.bdf'))
        op2.read_op2(os.path.join(basepath, 'flat_plate_metallic.op2'))
        op2_new = data_in_material_coord(bdf, op2)
        subcase = 3
        for vecname in force_vectors:
            name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
            ref_result = np.loadtxt(name)
            vector = getattr(op2_new, vecname)[subcase]
            data = vector.data
            eids = get_eids_from_op2_vector(vector)
            check = eids != 0
            if not np.allclose(data[:, check], ref_result, rtol=0.001):
                print('failed %r' % name)
                is_failed = True
        assert is_failed == False

    def test_stress(self):
        is_failed = False
        bdf = BDF(debug=False)
        op2 = OP2(debug=False)
        basepath = os.path.join(pkg_path, 'op2', 'test', 'test_flat_plate_metallic')
        bdf.read_bdf(os.path.join(basepath, 'flat_plate_metallic.bdf'))
        op2.read_op2(os.path.join(basepath, 'flat_plate_metallic.op2'))
        op2_new = data_in_material_coord(bdf, op2)
        subcase = 3
        for vecname in stress_vectors:
            name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
            ref_result = np.loadtxt(name)
            vector = getattr(op2_new, vecname)[subcase]
            data = vector.data
            eids = get_eids_from_op2_vector(vector)
            check = eids != 0
            if not np.allclose(data[:, check], ref_result, rtol=0.001):
                print('failed %r' % name)
                is_failed = True
        assert is_failed == False

    def test_strain(self):
        bdf = BDF(debug=False)
        op2 = OP2(debug=False)
        basepath = os.path.join(pkg_path, 'op2', 'test', 'test_flat_plate_metallic')
        bdf.read_bdf(os.path.join(basepath, 'flat_plate_metallic.bdf'))
        op2.read_op2(os.path.join(basepath, 'flat_plate_metallic.op2'))
        op2_new = data_in_material_coord(bdf, op2)
        subcase = 3
        is_failed = False
        for vecname in strain_vectors:
            name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
            ref_result = np.loadtxt(name)
            vector = getattr(op2_new, vecname)[subcase]
            data = vector.data
            eids = get_eids_from_op2_vector(vector)
            check = eids != 0
            np.savetxt(vecname + '.txt', data[:, check][0])
            if not np.allclose(data[:, check], ref_result, rtol=0.001):
                print('failed %r' % name)
                is_failed = True
        assert is_failed == False

if __name__ == '__main__':
    unittest.main()
