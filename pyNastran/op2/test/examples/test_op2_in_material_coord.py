import os
from pathlib import Path
import unittest
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.test.test_op2 import IS_PANDAS
from pyNastran.op2.data_in_material_coord import (
    data_in_material_coord,
    get_eids_from_op2_vector, force_vectors, stress_vectors,
    strain_vectors, get_eid_to_theta_rad, get_eid_to_theta_rad2)

PKG_PATH = Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'op2' / 'test' / 'examples' / 'coord_transform'


CASES = [
    ['test_flat_plate_metallic', 'flat_plate_metallic', 3],
    ['test_dummy_wing_metallic', 'dummy_wing_metallic', 1],
    ['test_flat_plate_composite', 'flat_plate_composite', 1],
]

RTOL = 0.01
ATOL = 0.01


class TestMaterialCoordReal(unittest.TestCase):
    def test_ctria6_mcid(self):
        log = get_logger(level='error')
        model = BDF(debug=False, log=log)

        pid = 10
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3

        nids = [1, 2, 3, 4]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])  # element system is along 1-2, so +y
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_pshell(pid, mid1=mid, t=0.1)
        model.add_mat1(mid, E, G, nu)

        eid_to_theta_deg_expected = {
            1 : 90.,
            2 : 0.,
            10 : 0.,
            11 : 90.,
        }
        model.add_ctria3(1, pid, [1, 2, 3], theta_mcid=0)
        model.add_ctria3(2, pid, [1, 2, 3], theta_mcid=0.0)
        _check_theta(model, eid_to_theta_deg_expected)
        model.elements = {}

        model.add_cquad4(10, pid, nids, theta_mcid=0.0)
        model.add_cquad4(11, pid, nids, theta_mcid=0)
        _check_theta(model, eid_to_theta_deg_expected)
        #bdf.card_count = {'GRID': 4, 'CQUAD4': 2,}
        #op2_new = data_in_material_coord(bdf, op2)

        x = 1

    def test_ctria6_mcid_real(self):
        log = get_logger(level='error')
        bdf = BDF(debug=False, log=log)
        op2 = OP2(debug=False, log=log, mode='msc')
        basepath = TEST_PATH / 'coord_sys_new'
        bdf_filename = basepath / 'coord_sys_new.bdf'
        op2_filename = basepath / 'coord_sys_new.op2'
        f06_filename = basepath / 'coord_sys_new.test.f06'
        bdf.read_bdf(bdf_filename)
        op2.read_op2(op2_filename)
        log.level = 'debug'
        op2_new = data_in_material_coord(bdf, op2)
        op2_new.write_f06(f06_filename)

        if IS_PANDAS:
            op2_new.build_dataframe()
        ctria6_strain = op2_new.op2_results.strain.ctria6_strain[1]
        cols = ['ElementID', 'NodeID', 'exx', 'eyy', 'exy']
        if IS_PANDAS:
            df = ctria6_strain.dataframe.reset_index()[cols]
            #print(df.to_string(float_format='%e'))
        x = 1

    def test_force(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, f'{prefix}.bdf')
            op2_filename = os.path.join(basepath, f'{prefix}.op2')
            bdf.read_bdf(bdf_filename)
            op2.read_op2(op2_filename)
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in force_vectors:
                vector = getattr(op2_new.op2_results.force, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, f'{vecname}_subcase_{subcase:02d}.txt')
                if not os.path.isfile(name):
                    raise AssertionError(f'Not found reference result {name}\n{print_bad_path(name)}')
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                #check = eids != 0
                #if np.any(~check):
                    #print(check)
                #rel =
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, 0::5, :], ref_result[0::5], rtol=RTOL, atol=ATOL)
                else:
                    assert data.shape[0] == 1, data.shape
                    datai = data[0, :, :]
                    diff = datai - ref_result
                    rel = diff / ref_result
                    msg = f'diff:\n{diff}\n'
                    msg += f'rel:\n{rel}\n'
                    assert np.allclose(datai, ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')

    def test_stress(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, prefix + '.bdf')
            op2_filename = os.path.join(basepath, prefix + '.op2')
            bdf.read_bdf(bdf_filename)
            op2.read_op2(op2_filename)
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in stress_vectors:
                vector = getattr(op2_new.op2_results.stress, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, f'{vecname}_subcase_{subcase:02d}.txt')
                if not os.path.isfile(name):
                    raise AssertionError(f'Not found reference result {name}\n{print_bad_path(name)}')
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = (eids != 0)
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, 0::10, :], ref_result[0::10], rtol=RTOL, atol=ATOL)
                    assert np.allclose(data[:, 1::10, :], ref_result[1::10], rtol=RTOL, atol=ATOL)
                else:
                    assert np.allclose(data, ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')

    def test_strain(self):
        log = get_logger(level='warning')
        for folder, prefix, subcase in CASES:
            bdf = BDF(debug=False, log=log)
            op2 = OP2(debug=False, log=log)
            basepath = os.path.join(TEST_PATH, folder)
            bdf_filename = os.path.join(basepath, prefix + '.bdf')
            op2_filename = os.path.join(basepath, prefix + '.op2')
            bdf.read_bdf(bdf_filename)
            op2.read_op2(op2_filename)
            op2_new = data_in_material_coord(bdf, op2)
            for vecname in strain_vectors:
                vector = getattr(op2_new.op2_results.strain, vecname).get(subcase)
                if vector is None:
                    continue
                name = os.path.join(basepath, f'{vecname}_subcase_{subcase:02d}.txt')
                if not os.path.isfile(name):
                    raise AssertionError(f'Not found reference result {name}\n{print_bad_path(name)}')
                ref_result = np.loadtxt(name)
                data = vector.data
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                if 'cquad8' in vecname:
                    assert np.allclose(data[:, 0::10, :], ref_result[0::10], rtol=RTOL, atol=ATOL)
                    assert np.allclose(data[:, 1::10, :], ref_result[1::10], rtol=RTOL, atol=ATOL)
                else:
                    assert np.allclose(data, ref_result, rtol=RTOL, atol=ATOL)
            #print('OK')

def _check_theta(model: BDF, eid_to_theta_deg_expected: dict[int, float]):
    model.cross_reference()

    model.log.level = 'warning'
    eid_to_theta_rad1 = get_eid_to_theta_rad(model, debug=True)
    eid_to_theta_rad2 = get_eid_to_theta_rad2(model, debug=True)
    model.log.level = 'info'

    for eid, theta in eid_to_theta_rad1.items():
        theta_deg_expected = eid_to_theta_deg_expected[eid]
        theta_deg = np.degrees(theta)
        assert np.allclose(theta_deg, theta_deg_expected), f'eid={eid} theta_deg={theta_deg} expected={theta_deg_expected}'
    model.uncross_reference()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
