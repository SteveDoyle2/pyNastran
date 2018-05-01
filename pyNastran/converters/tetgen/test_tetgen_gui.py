import os
import unittest
import warnings

import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')


from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
from pyNastran.utils.log import get_logger
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tetgen')
STL_PATH = os.path.join(pkg_path, 'converters', 'stl')


class TetgenGUI(TetgenIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = TetgenIO(self)
        self.build_fmts(['tetgen'], stop_on_failure=True)


class TestTetgenGUI(unittest.TestCase):

    def test_tetgen_geometry_01(self):
        log = get_logger(level='warning')
        geometry_filename = os.path.join(model_path, 'gear.smesh')

        test = TetgenGUI()
        test.log = log
        test.model.load_tetgen_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='tetgen')

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
        #base = 'tetgen_test_flipped.1'
        #m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
        #m.write_nastran(base + '.bdf')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

