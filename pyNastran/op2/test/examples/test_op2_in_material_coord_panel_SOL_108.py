from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.utils.log import get_logger
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import read_op2
from pyNastran.op2.data_in_material_coord import (
    data_in_material_coord,
    get_eids_from_op2_vector, force_vectors, stress_vectors,
    strain_vectors)
pkg_path = pyNastran.__path__[0]


CASES = [
    ['test_panel_SOL_108', 'panel_SOL_108_center', [1.0, 9.5]],
    ['test_panel_SOL_108', 'panel_SOL_108_corner', [1.0, 9.5]],
]

RTOL = 0.001

def calc_phasedeg(vec):
    return np.round(np.rad2deg(np.arctan2(vec.imag, vec.real)) % 360., 4)


class TestMaterialCoordComplex(unittest.TestCase):
    def test_force(self):
        log = get_logger(level='warning')
        is_failed = False
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2 = read_op2(
                os.path.join(basepath, prefix + '.op2'),
                debug=False, log=log,
                exclude_results=['stress', 'strain'],
            )
            try:
                op2_new = data_in_material_coord(bdf, op2)
            except ValueError as e:
                op2.log.error('failed rotating %r' % prefix)
                is_failed = True
                #continue
                raise

            for freq in freqs:
                for vecname in force_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, '%s_center_freq_%1.1f.txt' % (vecname, freq))
                    else:
                        name = os.path.join(basepath, '%s_corner_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
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
                    assert np.allclose(calc_phasedeg(data[:, :]), phase, rtol=RTOL)

        if is_failed:
            raise ValueError('see previous message')

    def test_stress(self):
        log = get_logger(level='warning')
        is_failed = False
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2 = read_op2(
                os.path.join(basepath, prefix + '.op2'),
                debug=False, log=log,
                exclude_results=['element_forces', 'strain'],
            )
            try:
                op2_new = data_in_material_coord(bdf, op2)
            except ValueError as e:
                op2.log.error('failed rotating %r' % prefix)
                is_failed = True
                #continue
                raise

            for freq in freqs:
                for vecname in stress_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, '%s_center_freq_%1.1f.txt' % (vecname, freq))
                    else:
                        name = os.path.join(basepath, '%s_corner_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
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
        if is_failed:
            raise ValueError('see previous message')

    def test_strain(self):
        log = get_logger(level='warning')
        is_failed = False
        for folder, prefix, freqs in CASES:
            bdf = BDF(debug=False, log=log)
            basepath = os.path.join(pkg_path, 'op2', 'test', 'examples', folder)
            bdf.read_bdf(os.path.join(basepath, prefix + '.bdf'))
            op2 = read_op2(
                os.path.join(basepath, prefix + '.op2'),
                debug=False, log=log,
                exclude_results=['element_forces', 'stress'],
            )

            try:
                op2_new = data_in_material_coord(bdf, op2)
            except ValueError as e:
                op2.log.error('failed rotating %r' % prefix)
                is_failed = True
                #continue
                raise

            for freq in freqs:
                for vecname in strain_vectors:
                    vector = getattr(op2_new, vecname).get(1)
                    if vector is None:
                        continue
                    if 'center' in prefix:
                        name = os.path.join(basepath, '%s_center_freq_%1.1f.txt' % (vecname, freq))
                    else:
                        name = os.path.join(basepath, '%s_corner_freq_%1.1f.txt' % (vecname, freq))
                    if not os.path.isfile(name):
                        raise AssertionError('Not found reference result {0}'.format(name))
                    ref_strain = np.loadtxt(name)
                    mag = ref_strain[:, 1::2]
                    phase = ref_strain[:, 2::2]
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
        if is_failed:
            raise ValueError('see previous message')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
