import os
import unittest

import pyNastran
from pyNastran.converters.dev.vrml.vrml import read_vrml
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'dev', 'vrml')


class TestVrml(unittest.TestCase):
    def test_vrml_1(self):
        vrml_filename = os.path.join(MODEL_PATH, 'pyramid_sphere.wrl')
        log = None
        read_vrml(vrml_filename, log=log)

if __name__ == '__main__':
    unittest.main()
