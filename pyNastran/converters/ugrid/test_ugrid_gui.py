from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.converters.ugrid.ugrid_io import UGRID_IO
from pyNastran.converters.ugrid.surf_io import SurfIO
from pyNastran.gui.testing_methods import GUIMethods
#from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid
#from pyNastran.converters.ugrid.ugrid3d_to_nastran import ugrid3d_to_nastran
#from pyNastran.converters.ugrid.ugrid3d_to_tecplot import ugrid3d_to_tecplot_filename, ugrid_to_tecplot

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(pkg_path, '..', 'models')

class UGRID_GUI(UGRID_IO, SurfIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        UGRID_IO.__init__(self)
        SurfIO.__init__(self)

class TestUgridGUI(unittest.TestCase):
    def test_ugrid_gui_01(self):
        ugrid_filename = os.path.join(pkg_path, 'converters', 'ugrid', 'models',
                                      'two_blade_wake_sym_extended.surf')
        test = UGRID_GUI()
        test.load_surf_geometry(ugrid_filename, dirname=None)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
