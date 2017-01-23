import os

import warnings
import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')


from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.tetgen.tetgen_io import TetgenIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tetgen')
STL_PATH = os.path.join(pkg_path, 'converters', 'stl')

import unittest

class TetgenGUI(TetgenIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        TetgenIO.__init__(self)


class TestTetgenGUI(unittest.TestCase):

    def test_tetgen_geometry_01(self):
        geometry_filename = os.path.join(model_path, 'gear.smesh')
        dirname = None

        test = TetgenGUI()
        test.load_tetgen_geometry(geometry_filename, dirname)

    def test_tetgen_geometry_02(self):
        regenerate_files = False
        if regenerate_files: # pragma: no cover
            from pyNastran.converters.stl.stl import STL
            m1 = STL()
            m1.read_stl(os.path.join(STL_PATH, 'sphere.stl'))
            #m1.flip_normals()
            m1.write_stl('tetgen_test.stl')
            del m1
            os.system('tetgen.exe -pqcvVqY tetgen_test.stl')

        test = TetgenGUI()
        base = os.path.join(model_path, 'tetgen_test.1')
        test.load_tetgen_geometry(base + '.smesh', dirname=None)
        test.load_tetgen_geometry(base + '.ele', dirname=None)
        #base = 'tetgen_test_flipped.1'
        #m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
        #m.write_nastran(base + '.bdf')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

