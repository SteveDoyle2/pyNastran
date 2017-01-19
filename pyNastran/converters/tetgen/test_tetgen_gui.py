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

import unittest

class TetgenGUI(TetgenIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        TetgenIO.__init__(self)


class TestTetgenGUI(unittest.TestCase):

    def test_tetgen_geometry(self):
        return
        geometry_filename = os.path.join(model_path, 'sphere.stl')
        dirname = None

        test = STL_GUI()
        #test.load_nastran_geometry(geometry_filename, None)
        test.load_stl_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

