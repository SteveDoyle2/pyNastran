import os
import unittest
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.dev.obj.obj_io import ObjIO
from pyNastran.converters.dev.obj.obj import OBJ, read_obj


PKG_PATH = pyNastran.__path__[0]
model_path = os.path.join(PKG_PATH, 'converters', 'dev', 'obj') #, 'models')


class ObjGUI(ObjIO, FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = ObjIO(self)
        self.build_fmts(['obj'], stop_on_failure=True)


class TestObjGUI(unittest.TestCase):

    def test_obj_geometry_01(self):
        """tests the ascii shuttle model"""
        log = get_logger(level='warning', encoding='utf-8')
        geometry_filename = os.path.join(model_path, 'shuttle.obj')

        test = ObjGUI()
        test.log = log
        #test.load_obj_geometry(geometry_filename)
        test.on_load_geometry(geometry_filename, geometry_format='obj', raise_error=True)

    def test_obj_geometry_02(self):
        """tests the ascii shuttle model"""
        log = get_logger(level='warning', encoding='utf-8')
        obj_filename = os.path.join(model_path, 'shuttle.obj')
        obj_filename2 = os.path.join(model_path, 'shuttle2.obj')

        model = read_obj(obj_filename, log=log)
        model.write_obj(obj_filename2)
        model2 = read_obj(obj_filename2, log=log)
        os.remove(obj_filename2)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

