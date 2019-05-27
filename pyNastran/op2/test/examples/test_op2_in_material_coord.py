import os
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.data_in_material_coord import (
    data_in_material_coord,
    get_eids_from_op2_vector, force_vectors, stress_vectors,
    strain_vectors)
pkg_path = pyNastran.__path__[0]


CASES = [
    ['test_flat_plate_metallic', 'flat_plate_metallic', 3],
    ['test_dummy_wing_metallic', 'dummy_wing_metallic', 1],
    ['test_flat_plate_composite', 'flat_plate_composite', 1],
]

RTOL = 0.01
ATOL = 0.01


class TestMaterialCoordReal(unittest.TestCase):
    def test_force(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in force_vectors:
                vector = getattr(op2_new, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    raise AssertionError('Not found reference result {0}\n{1}'.format(
                        name, print_bad_path(name)))
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, check][:, 0::5, :], ref_result[0::5], rtol=RTOL, atol=ATOL)
                else:
                    assert np.allclose(data[:, check], ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')

    def test_stress(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in stress_vectors:
                vector = getattr(op2_new, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    raise AssertionError('Not found reference result {0}\n{1}'.format(
                        name, print_bad_path(name)))
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, check][:, 0::10, :], ref_result[0::10], rtol=RTOL, atol=ATOL)
                    assert np.allclose(data[:, check][:, 1::10, :], ref_result[1::10], rtol=RTOL, atol=ATOL)
                else:
                    assert np.allclose(data[:, check], ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')

    def test_strain(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in strain_vectors:
                vector = getattr(op2_new, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, '{0}_subcase_{1:02d}.txt'.format(vecname, subcase))
                if not os.path.isfile(name):
                    raise AssertionError('Not found reference result {0}\n{1}'.format(
                        name, print_bad_path(name)))
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, check][:, 0::10, :], ref_result[0::10], rtol=RTOL, atol=ATOL)
                    assert np.allclose(data[:, check][:, 1::10, :], ref_result[1::10], rtol=RTOL, atol=ATOL)
                else:
                    assert np.allclose(data[:, check], ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
