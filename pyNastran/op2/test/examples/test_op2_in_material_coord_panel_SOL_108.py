import os
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import read_op2
from pyNastran.op2.data_in_material_coord import (
    data_in_material_coord,
    get_eids_from_op2_vector, force_vectors, stress_vectors,
    strain_vectors)
PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'op2', 'test', 'examples', 'coord_transform')


CASES = [
    ['test_panel_SOL_108', 'panel_SOL_108_center', [1.0, 9.5]],
    ['test_panel_SOL_108', 'panel_SOL_108_corner', [1.0, 9.5]],
    ['test_panel_SOL_108', 'panel_SOL_108_center_tria', [1.0]],
    ['test_panel_SOL_108', 'panel_SOL_108_corner_tria', [1.0]],
]

RTOL = 0.001

def calc_phasedeg(vec):
    out = np.round(np.rad2deg(np.arctan2(vec.imag, vec.real)) % 360., 4)
    out[np.isclose(vec.real, 0)] = 0
    return out



class TestMaterialCoordComplex(unittest.TestCase):
    def test_force(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, f'{prefix}.bdf')
            op2_filename = os.path.join(basepath, f'{prefix}.op2')
            bdf.read_bdf(bdf_filename)
            op2 = read_op2(
                op2_filename,
                debug=False, log=log,
                exclude_results=['stress', 'strain'],
            )
            op2_new = data_in_material_coord(bdf, op2)

            for freq in freqs:
                for vecname in force_vectors:
                    vector = getattr(op2_new.op2_results.force, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, f'{vecname}_center_freq_{freq:1.1f}.txt')
                    else:
                        name = os.path.join(basepath, f'{vecname}_corner_freq_{freq:1.1f}.txt')
                    if not os.path.isfile(name):
                        raise AssertionError(f'Not found reference result {name}')
                    ref_force = np.loadtxt(name)
                    mag = ref_force[0::2]
                    phase = ref_force[1::2]
                    if freq == 1.0:
                        data = vector.data[0]
                    elif freq == 9.5:
                        data = vector.data[17]
                    #eids = get_eids_from_op2_vector(vector)
                    #check = eids != 0
                    assert np.allclose(np.abs(data[:, :]), mag, rtol=RTOL)
                    assert np.allclose(calc_phasedeg(data), phase, rtol=RTOL)

    def test_stress(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, f'{prefix}.bdf')
            op2_filename = os.path.join(basepath, f'{prefix}.op2')
            bdf.read_bdf(bdf_filename)
            op2 = read_op2(
                op2_filename,
                debug=False, log=log,
                exclude_results=['element_forces', 'strain'],
            )
            op2_new = data_in_material_coord(bdf, op2)

            for freq in freqs:
                for vecname in stress_vectors:
                    vector = getattr(op2_new.op2_results.stress, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, f'{vecname}_center_freq_{freq:1.1f}.txt')
                    else:
                        name = os.path.join(basepath, f'{vecname}_corner_freq_{freq:1.1f}.txt')
                    if not os.path.isfile(name):
                        raise AssertionError(f'Not found reference result {name}')
                    ref_stress = np.loadtxt(name)
                    mag = ref_stress[:, 1::2]
                    phase = ref_stress[:, 2::2]
                    if freq == 1.0:
                        data = vector.data[0]
                    elif freq == 9.5:
                        data = vector.data[17]
                    eids = get_eids_from_op2_vector(vector)
                    check = eids != 0
                    if 'cquad8' in vecname:
                        assert np.allclose(np.abs(data[check][0::10, :]), mag[0::10], rtol=RTOL)
                        assert np.allclose(calc_phasedeg(data[check][1::10, :]), phase[1::10], rtol=RTOL)
                    else:
                        assert np.allclose(np.abs(data[check]), mag, rtol=RTOL)
                        assert np.allclose(calc_phasedeg(data[check]), phase, rtol=RTOL)

    def test_strain(self):
        log = get_logger(level='warning')
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, f'{prefix}.bdf')
            op2_filename = os.path.join(basepath, f'{prefix}.op2')
            bdf.read_bdf(bdf_filename)
            op2 = read_op2(
                op2_filename,
                debug=False, log=log,
                exclude_results=['element_forces', 'stress'],
            )
            op2_new = data_in_material_coord(bdf, op2)

            for freq in freqs:
                for vecname in strain_vectors:
                    vector = getattr(op2_new.op2_results.strain, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, f'{vecname}_center_freq_{freq:1.1f}.txt')
                    else:
                        name = os.path.join(basepath, f'{vecname}_corner_freq_{freq:1.1f}.txt')
                    if not os.path.isfile(name):
                        raise AssertionError(f'Not found reference result {name}')
                    ref_strain = np.loadtxt(name)
                    mag = ref_strain[:, 1::2]
                    phase = ref_strain[:, 2::2]
                    if freq == 1.0:
                        data = vector.data[0]
                    elif freq == 9.5:
                        data = vector.data[17]
                    eids = get_eids_from_op2_vector(vector)
                    check = eids != 0
                    assert np.allclose(np.abs(data[check]), mag, rtol=RTOL)
                    phase[np.isclose(mag, 0)] = 0
                    assert np.allclose(calc_phasedeg(data[check]), phase, rtol=RTOL)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
