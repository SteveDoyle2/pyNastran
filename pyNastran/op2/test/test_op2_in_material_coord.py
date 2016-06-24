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


CASES = [
    ['test_flat_plate_metallic', 'flat_plate_metallic', 3],
    ['test_dummy_wing_metallic', 'dummy_wing_metallic', 1],
    ]

SUMDIFFTOL = 0.2


class TestMaterialCoord(unittest.TestCase):
    def test_force(self):
        for folder, prefix, subcase in CASES:
            is_failed = False
            bdf = BDF(debug=False)
            op2 = OP2(debug=False)
            basepath = os.path.join(pkg_path, 'op2', 'test', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in force_vectors:
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    continue
                ref_result = np.loadtxt(name)
                vector = getattr(op2_new, vecname)[subcase]
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    test = abs((data[:, check][:, 0::5, :] - ref_result[0::5]).sum())
                else:
                    test = abs((data[:, check] - ref_result).sum())
                if test > SUMDIFFTOL:
                    print('FAILED %r' % name)
                    is_failed = True
            assert is_failed == False

    def test_stress(self):
        for folder, prefix, subcase in CASES:
            is_failed = False
            bdf = BDF(debug=False)
            op2 = OP2(debug=False)
            basepath = os.path.join(pkg_path, 'op2', 'test', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in stress_vectors:
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    continue
                ref_result = np.loadtxt(name)
                vector = getattr(op2_new, vecname)[subcase]
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    test = abs((data[:, check][:, 0::10, :] - ref_result[0::10]).sum())
                    test += abs((data[:, check][:, 1::10, :] - ref_result[1::10]).sum())
                else:
                    test = abs((data[:, check] - ref_result).sum())
                if test > SUMDIFFTOL:
                    print('FAILED %r' % name)
                    is_failed = True
            assert is_failed == False

    def test_strain(self):
        for folder, prefix, subcase in CASES:
            is_failed = False
            bdf = BDF(debug=False)
            op2 = OP2(debug=False)
            basepath = os.path.join(pkg_path, 'op2', 'test', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in strain_vectors:
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    continue
                ref_result = np.loadtxt(name)
                vector = getattr(op2_new, vecname)[subcase]
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    test = abs((data[:, check][:, 0::10, :] - ref_result[0::10]).sum())
                    test += abs((data[:, check][:, 1::10, :] - ref_result[1::10]).sum())
                else:
                    test = abs((data[:, check] - ref_result).sum())
                if test > SUMDIFFTOL:
                    print('FAILED %r' % name)
                    is_failed = True
            assert is_failed == False


if __name__ == '__main__':
    unittest.main()
