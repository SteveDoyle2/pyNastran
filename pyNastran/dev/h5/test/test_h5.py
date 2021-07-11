import os
import unittest

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.dev.h5.read_h5 import pyNastranH5

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

class TestH5(unittest.TestCase):
    def test_static_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_modes_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_time_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'time_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_modes_complex_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_freq_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_freq_elements2_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    def test_time_thermal_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\static_elements.h5',
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\modes_elements.h5',
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\time_elements.h5',
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\modes_complex_elements.h5',
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\freq_elements.h5',
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\freq_elements2.h5',
        ###r'C:\NASA\m4\formats\git\pyNastran\models\elements\loadstep_elements.h5',  # no nonlinear examples
        #r'C:\NASA\m4\formats\git\pyNastran\models\elements\time_thermal_elements.h5',

if __name__ == '__main__':
    unittest.main()
