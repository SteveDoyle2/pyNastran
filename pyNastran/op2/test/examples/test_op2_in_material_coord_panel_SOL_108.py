from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.utils.log import get_logger
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.data_in_material_coord import (
    data_in_material_coord,
    get_eids_from_op2_vector, force_vectors, stress_vectors,
    strain_vectors)
pkg_path = pyNastran.__path__[0]


CASES = [
         ['test_panel_SOL_108', 'panel_SOL_108_center', [1.0, 9.5]],
         ['test_panel_SOL_108', 'panel_SOL_108_corner', [1.0, 9.5]],
        ]

RTOL = 0.01
ATOL = 0.01

def calc_phasedeg(vec):
    return np.round(np.rad2deg(np.arctan2(vec.imag, vec.real)) % 360., 4)


class TestMaterialCoordComplex(unittest.TestCase):
    def test_force(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for freq in freqs:
                for vecname in force_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    name = os.path.join(basepath, '%s_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
                    ref_force = np.loadtxt(name)
                    mag = ref_force[0::2]
                    phase = ref_force[1::2]
                    data = vector.data
                    eids = get_eids_from_op2_vector(vector)
                    check = eids != 0
                    if 'cquad8' in vecname:
                        assert np.allclose(np.abs(data[:, check][:, 0::5, :]), mag[0::5], rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check][:, 0::5, :]), phase[0::5], rtol=RTOL, atol=ATOL)
                    else:
                        assert np.allclose(np.abs(data[:, check]), mag, rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check]), phase, rtol=RTOL, atol=ATOL)

    def test_stress(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for freq in freqs:
                for vecname in stress_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    name = os.path.join(basepath, '%s_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
                    ref_stress = np.loadtxt(name)
                    mag = ref_stress[:, 1::2]
                    phase = ref_stress[:, 2::2]
                    data = vector.data
                    eids = get_eids_from_op2_vector(vector)
                    check = eids != 0
                    if 'cquad8' in vecname:
                        assert np.allclose(np.abs(data[:, check][:, 0::10, :]), mag[0::10], rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check][:, 1::10, :]), phase[1::10], rtol=RTOL, atol=ATOL)
                    else:
                        assert np.allclose(np.abs(data[:, check]), mag, rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check]), phase, rtol=RTOL, atol=ATOL)

    def test_strain(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2.read_op2(os.path.join(basepath, prefix + '.op2'))
            op2_new = data_in_material_coord(bdf, op2)
            for freq in freqs:
                for vecname in strain_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    name = os.path.join(basepath, '%s_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
                    ref_strain = np.loadtxt(name)
                    mag = ref_strain[:, 1::2]
                    phase = ref_strain[:, 2::2]
                    data = vector.data
                    eids = get_eids_from_op2_vector(vector)
                    check = eids != 0
                    if 'cquad8' in vecname:
                        assert np.allclose(np.abs(data[:, check][:, 0::10, :]), mag[0::10], rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check][:, 1::10, :]), phase[1::10], rtol=RTOL, atol=ATOL)
                    else:
                        assert np.allclose(np.abs(data[:, check]), mag, rtol=RTOL, atol=ATOL)
                        assert np.allclose(calc_phasedeg(data[:, check]), phase, rtol=RTOL, atol=ATOL)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
