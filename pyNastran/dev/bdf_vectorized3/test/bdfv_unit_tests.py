from collections import defaultdict
import pathlib
import unittest
import numpy as np
import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.dev.bdf_vectorized3.bdf import read_bdf as read_bdf_vectorized

PKG_PATH = pathlib.Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'bdf' / 'test'
MODEL_PATH = PKG_PATH / '..' / 'models'

class TestBdfVectorized3(unittest.TestCase):
    def test_bwb(self):
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        bdf_filename_out = MODEL_PATH / 'bwb' / 'bwb_saero_vectorized.bdf'
        modelv = read_bdf_vectorized(
            bdf_filename, validate=True, xref=True, punch=False,
            save_file_structure=False, skip_cards=None, read_cards=None,
            encoding=None, log=None, debug=False, mode='msc')
        model = read_bdf(
            bdf_filename, validate=True, xref=True, punch=False,
            save_file_structure=False, skip_cards=None, read_cards=None,
            encoding=None, log=None, debug=None, mode='msc')

        area, mass = get_compare_model_vectorized(model, modelv)

        area_ctria3 = modelv.ctria3.area().sum()
        area_cquad4 = modelv.cquad4.area().sum()

        assert np.allclose(area['CTRIA3'], 7842.639526013481), area['CTRIA3']
        assert np.allclose(area['CQUAD4'], 2248660.627265879), area['CQUAD4']

        assert np.allclose(mass['CTRIA3'], 318.61913895889063), mass['CTRIA3']
        assert np.allclose(mass['CQUAD4'], 98307.632265389750), mass['CQUAD4']


        assert np.allclose(area_ctria3, 7842.639526013481), area_ctria3
        assert np.allclose(area_cquad4, 2248660.627265879), area_cquad4

        mass_ctria3 = modelv.ctria3.mass().sum()
        mass_cquad4 = modelv.cquad4.mass().sum()
        assert np.allclose(mass_ctria3, 318.61913895889063), mass_ctria3
        assert np.allclose(mass_cquad4, 98307.632265389750), mass_cquad4
        modelv.write_bdf(bdf_filename_out)

def get_compare_model_vectorized(model, modelv):
    area = defaultdict(list)
    mass = defaultdict(list)
    for eid, elem in model.elements.items():
        etype = elem.type
        if etype in {'CTRIA3', 'CQUAD4'}:
            area[etype] += elem.Area()
            mass[etype] += elem.Mass()
        else:
            continue
            #asdf
    return area, mass

if __name__ == '__main__':
    unittest.main()
