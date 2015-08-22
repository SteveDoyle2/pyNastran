import os

from pyNastran.gui.testing_methods import GUIMethods
from pyNastran.converters.shabp.shabp_io import ShabpIO
import pyNastran

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'shabp')

import unittest

class ShabpGUI(ShabpIO, GUIMethods):
    def __init__(self):
        GUIMethods.__init__(self)
        ShabpIO.__init__(self)


class TestShabpGUI(unittest.TestCase):

    def test_shabp_geometry(self):
        return
        #geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        dirname = None

        test = ShabpGUI()


        test.load_shabp_geometry('models/NAC6.INP', '')
        #test.load_nastran_geometry(geometry_filename, None)
        #test.load_shabp_geometry(geometry_filename, dirname)

    def test_shabp_results(self):
        pass
        #geometry_filename = os.path.join(model_path, 'M100.inp')
        #agps_filename = os.path.join(model_path, 'agps')
        #out_filename = os.path.join(model_path, 'panair.out')
        #dirname = None

        #test = ShabpGUI()

        #test.load_panair_geometry(geometry_filename, dirname)
        #test.load_panair_results(agps_filename, dirname)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

