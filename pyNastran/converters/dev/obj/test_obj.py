import os
import unittest

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.dev.obj.obj_io import ObjIO
from pyNastran.converters.dev.obj.obj import OBJ
from pyNastran.utils.log import get_logger


pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, 'converters', 'dev', 'obj') #, 'models')


class ObjGUI(ObjIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        ObjIO.__init__(self, self)


class TestObjGUI(unittest.TestCase):

    def test_obj_geometry_01(self):
        """tests the ascii three plugs model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(model_path, 'shuttle.obj')

        test = ObjGUI()
        test.log = log
        test.load_obj_geometry(geometry_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

