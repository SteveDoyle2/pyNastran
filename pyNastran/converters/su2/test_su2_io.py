import os

import warnings
import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')


from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.su2.su2_io import SU2_IO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'su2')

import unittest

class SU2_GUI(STL_IO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        SU2_IO.__init__(self)


class SU2_GUITest(unittest.TestCase):

    def test_su2_geometry(self):
        geometry_filename = os.path.join(model_path, 'mesh_naca0012_inv.su2')
        dirname = None

        test = SU2_GUI()
        test.load_su2_geometry(geometry_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

