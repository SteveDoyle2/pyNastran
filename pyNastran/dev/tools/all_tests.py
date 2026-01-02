import numpy as np
import unittest
from pyNastran.dev.tools.inrel.test_inrel import TestInrel
from pyNastran.dev.tools.test_pressure_map import TestPressureMap
from pyNastran.dev.tools.huth_stiffness import get_huth

class TestHuth(unittest.TestCase):
    def test_huth(self):
        diameter = 0.1875  # inch
        t1 = 0.040  # inch
        t2 = 0.063  # inch
        E1 = 10500000  # psi
        E2 = 10500000  # psi
        E3 = 16000000  # psi

        # ----------------------------------------------------------
        # metallic_metallic_riveted
        # cis = 'single_shear'
        k_huth = get_huth(diameter, t1, t2, E1, E2, E3,
                          shear_type='single',
                          connection_type='rivet',
                          mat1_type='metallic',
                          mat2_type='metallic')
        # print(f'k_huth = {k_huth:g}')  # 147,424  # psi
        assert np.allclose(k_huth, 147424), k_huth

        cis = 'double_shear'
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='rivet',
            mat1_type='metallic',
            mat2_type='metallic')
        # print(f'k_huth = {k_huth:g}')  # 365,986  # psi
        assert np.allclose(k_huth, 365896.0556), k_huth
        # --------------------------
        #  metallic_metallic_bolted
        cis = 'single_shear'
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='metallic',
            mat2_type='metallic')
        # print(f'k_huth = {k_huth:g}')  # 152,589  # psi
        assert np.allclose(k_huth, 152589), k_huth

        cis = 'double_shear'
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='metallic',
            mat2_type='metallic')
        # print(f'k_huth = {k_huth:g}')  # 378,714  # psi
        assert np.allclose(k_huth, 378714), k_huth

        # -------------------------------------------
        # metalic_composite_bolted
        cis = 'single_shear'
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='composite',
            mat2_type='composite')
        # print(f'k_huth = {k_huth:g}')
        assert np.allclose(k_huth, 108992), k_huth

        cis = 'double_shear'  # n=2
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='composite',
            mat2_type='composite')
        # print(f'k_huth = {k_huth:g}')
        assert np.allclose(k_huth, 270510), k_huth
        # -------------------------------------------
        # metalic_composite_bolted
        cis = 'single_shear'
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='metallic',
            mat2_type='composite')
        # print(f'k_huth = {k_huth:g}')
        assert np.allclose(k_huth, 127157), k_huth

        cis = 'double_shear'  # n=2
        k_huth = get_huth(
            diameter, t1, t2, E1, E2, E3, cis,
            connection_type='bolt',
            mat1_type='metallic',
            mat2_type='composite')

        # print(f'k_huth = {k_huth:g}')
        assert np.allclose(k_huth, 315595), k_huth


if __name__ == '__main__':  # pragma: no covr
    unittest.main()
