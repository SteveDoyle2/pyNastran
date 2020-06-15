"""tests S/HABP"""
import os
import unittest
import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.converters.shabp.shabp import read_shabp
from pyNastran.converters.shabp.shabp_results import ShabpOut

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'shabp', 'models')

class TestShabp(unittest.TestCase):
    """tests S/HABP"""
    def test_shabp_1(self):
        """tests nose.mk5"""
        log = SimpleLogger(level='info', encoding='utf-8')
        shabp_filename = os.path.join(MODEL_PATH, 'nose', 'noseX_working.mk5')
        model = read_shabp(shabp_filename, log=log)

        npatches_expected = 1
        assert len(model.X) == npatches_expected, f'npatches_expected={npatches_expected} len(model.X)={len(model.X)}'
        assert len(model.Y) == npatches_expected, f'npatches_expected={npatches_expected} len(model.X)={len(model.Y)}'
        assert len(model.Z) == npatches_expected, f'npatches_expected={npatches_expected} len(model.X)={len(model.Z)}'


        #print(f'component_name_to_patch = {model.component_name_to_patch}')
        #print(f'patch_to_component_num = {model.patch_to_component_num}')
        #print(f'component_to_params = {model.component_to_params}')
        #print(f'component_num_to_name = {model.component_num_to_name}')
        #print(f'component_name_to_num = {model.component_name_to_num}')
        assert model.component_name_to_patch == {'ellipse': [1]}, model.component_name_to_patch
        assert model.patch_to_component_num == {0: 1}, model.patch_to_component_num
        assert model.component_to_params == {0: [5, 5, 1, 0, 0, 0.0, 1.0, 0.0, 1.0, 3.0, 3.0]}, model.component_to_params
        assert model.component_num_to_name == {0: 'ellipse'}, model.component_num_to_name
        assert model.component_name_to_num == {'ellipse': 0}, model.component_name_to_num

        areas = model.get_area_by_patch()
        assert np.allclose(areas, [266.47640991]), areas

        areas = model.get_area_by_component()
        assert np.allclose(areas['ellipse'], 266.47640991), areas

        areas, lengths = model.get_area_xlength_by_component()
        assert np.allclose(areas['ellipse'], 266.47640991), areas
        assert np.allclose(lengths['ellipse'], 20.0), lengths

        #self.title = ''
        #self.header = ''
        #self.shabp_cases = {}

    def test_shabp_2(self):
        """tests the flap"""
        log = SimpleLogger(level='info', encoding='utf-8')
        shabp_infilename = os.path.join(MODEL_PATH, 'flap', 'flap_inviscid.mk5')
        shabp_outfilename = os.path.join(MODEL_PATH, 'flap', 'SHABP.OUT')

        #test.model.load_shabp_geometry(shabp_infilename)
        #test.on_load_geometry(shabp_infilename, geometry_format='shabp', raise_error=True)
        model = read_shabp(shabp_infilename, read_special_routines=True,
                           log=log, debug=None)

        #print(f'component_name_to_patch = {model.component_name_to_patch}')
        #print(f'patch_to_component_num = {model.patch_to_component_num}')
        #print(f'component_to_params = {model.component_to_params}')
        #print(f'component_num_to_name = {model.component_num_to_name}')
        #print(f'component_name_to_num = {model.component_name_to_num}')

        assert model.component_name_to_patch == {'COMP1': [1], 'FLAP': [2]}, model.component_name_to_patch
        assert model.patch_to_component_num == {0: 1, 1: 2}, model.patch_to_component_num
        assert model.component_to_params == {0: [3, 1, 1, 0, 0, 0.0, 1.0, 0.0, 1.0, 3.0, 3.0], 1: [3, 1, 1, 0, 0, 0.0, 1.0, 0.0, 1.0, 3.0, 3.0]}, model.component_to_params
        assert model.component_num_to_name == {0: 'COMP1', 1: 'FLAP'}, model.component_num_to_name
        assert model.component_name_to_num == {'COMP1': 0, 'FLAP': 1}, model.component_name_to_num

        areas = model.get_area_by_patch()
        assert np.allclose(areas, [50., 124.6875]), areas

        areas = model.get_area_by_component()
        assert np.allclose(areas['COMP1'], 124.6875), areas
        assert np.allclose(areas['FLAP'], 50.0), areas

        areas, lengths = model.get_area_xlength_by_component()
        assert np.allclose(areas['COMP1'], 124.6875), areas
        assert np.allclose(areas['FLAP'], 50.0), areas
        assert np.allclose(lengths['COMP1'], 24.9375), lengths
        assert np.allclose(lengths['FLAP'], 34.9375), lengths

        out_model = ShabpOut(model, log=log)
        Cpd, unused_deltad = out_model.read_shabp_out(shabp_outfilename)

        #test.on_load_results(shabp_outfilename)


if __name__ == '__main__':   # pragma: no cover
    unittest.main()
