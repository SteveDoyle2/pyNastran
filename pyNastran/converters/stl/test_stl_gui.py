import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.stl.stl_io import STL_IO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'stl')

import unittest

class STL_GUI(STL_IO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        STL_IO.__init__(self)


class STL_GUITest(unittest.TestCase):

    def test_stl_geometry(self):
        geometry_filename = os.path.join(model_path, 'sphere.stl')
        dirname = None

        test = STL_GUI()
        #test.load_nastran_geometry(geometry_filename, None)
        test.load_stl_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

