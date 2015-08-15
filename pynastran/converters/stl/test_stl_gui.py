import os

from pyNastran.gui.testing_methods import add_dummy_gui_functions
from pyNastran.converters.stl.stl_io import STL_IO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'stl')

import unittest

class STLGUI(unittest.TestCase):

    def test_stl_geometry(self):
        geometry_filename = os.path.join(model_path, 'sphere.stl')
        dirname = None

        test = STL_IO()
        test.is_nodal = True
        test.is_centroidal = True

        add_dummy_gui_functions(test)

        #test.load_nastran_geometry(geometry_filename, None)
        test.load_stl_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

