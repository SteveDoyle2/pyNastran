import os
import unittest
import warnings

import numpy as np
from cpylog import get_logger
warnings.simplefilter('always')
np.seterr(all='raise')


from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.tetgen.tetgen_to_usm3d import tetgen_to_usm3d
from pyNastran.converters.tetgen.tetgen import read_tetgen
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
import pyNastran

PKG_PATH = pyNastran.__path__[0]
model_path = os.path.join(PKG_PATH, 'converters', 'tetgen')
STL_PATH = os.path.join(PKG_PATH, 'converters', 'stl')


class TetgenGUI(TetgenIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = TetgenIO(self)
        self.build_fmts(['tetgen'], stop_on_failure=True)


class TestTetgenGUI(unittest.TestCase):
    def test_tetgen_to_usm3d_01(self):
        #geometry_filename = os.path.join(model_path, 'tetgen_test.1.smesh')
        base = os.path.join(model_path, 'tetgen_test.1')
        tetgen_to_usm3d(base)


    def test_tetgen_geometry_01(self):
        log = get_logger(level='warning')
        geometry_filename = os.path.join(model_path, 'gear.smesh')

        test = TetgenGUI()
        test.log = log
        test.model.load_tetgen_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='tetgen', raise_error=True)

    def test_tetgen_geometry_02(self):
        log = get_logger(level='warning')
        regenerate_files = False
        if regenerate_files: # pragma: no cover
            from pyNastran.converters.stl.stl import read_stl
            model = read_stl(os.path.join(STL_PATH, 'sphere.stl'))
            #model.flip_normals()
            model.write_stl('tetgen_test.stl')
            del model
            os.system('tetgen.exe -pqcvVqY tetgen_test.stl')

        test = TetgenGUI()
        test.log = log
        base = os.path.join(model_path, 'tetgen_test.1')
        test.model.load_tetgen_geometry(base + '.smesh')
        test.model.load_tetgen_geometry(base + '.ele')

        model2 = read_tetgen(base, dimension_flag=2, log=None, debug=False)
        model2.write_nastran(base + '.bdf')

        model3 = read_tetgen(base, dimension_flag=3, log=None, debug=False)
        model3.write_nastran(base + '.bdf')

        #base = 'tetgen_test_flipped.1'
        #m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
        #m.write_nastran(base + '.bdf')
        os.remove(base + '.bdf')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

