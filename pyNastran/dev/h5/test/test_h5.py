import os
import unittest
import vtk

import pyNastran
from pyNastran.utils import print_bad_path

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
_VTK_VERSION = vtk.vtkVersion.GetVTKVersion()
VTK_VERSION_MAJOR = [int(val) for val in _VTK_VERSION.split('.')][0]
if VTK_VERSION_MAJOR >= 9:
    from pyNastran.dev.h5.read_h5 import pyNastranH5

class TestH5(unittest.TestCase):
    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_static_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_modes_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_time_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'time_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_modes_complex_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_freq_elements_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
    def test_freq_elements2_h5(self):
        h5_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.h5')
        assert os.path.exists(h5_filename), print_bad_path(h5_filename)
        model = pyNastranH5(add_aero=True, add_constraints=True,
                            add_results=True, subcases=None)
        model.read_h5_nastran(h5_filename, subcases=None)

    @unittest.skipIf(VTK_VERSION_MAJOR < 9, "skipping because test requires VTK >= 9")
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

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
